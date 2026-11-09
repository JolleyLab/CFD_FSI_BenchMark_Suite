#!/usr/bin/env python3
"""Build one benchmark geometry at one or more resolutions from a case spec (JSON).

For each resolution:
  cases/<name>/<res>/
    mesh/                       svMultiPhysics mesh-complete (volume VTU, exterior, walls, mesh-surfaces/*.vtp)
    solver.xml                  rigid-wall CFD case (svMultiPhysics XML), BCs bound by face name
    <cap>.flow                  inflow waveform(s)
    mfem/<name>_<res>.msh       Gmsh 2.2 binary for MFEM (+ _attributes.json: boundary attribute -> face name)
    quality.json / quality.md   element quality report
    case.json                   provenance, element counts, edge length, timings, BC summary

Usage: build_case.py spec.json [--res coarse,medium,fine] [--out cases] [--calibrate-only]

Credits: meshing by SlicerExtension-VMTK CFD Mesh Generator (VMTK, Netgen); mesh-complete, solver.xml
and .flow files by SlicerSimVascular svmeshcomplete / svmpsetup; RCR outlets after Westerhof et al. 2009
(doi:10.1007/s11517-008-0359-2); flow split by area after Murray 1926 (doi:10.1073/pnas.12.3.207). See REFERENCES.md.
"""
import argparse, json, os, sys, time, shutil, math, logging, resource
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import meshlib, bclib, quality, to_mfem
from svmpsetup.case import write_case
from svromsetup.boundary_conditions import Inflow, RCR, Resistance

log = logging.getLogger("build_case")


def wait_for_memory(target_elements, mb_per_element=0.9e-3, reserve_mb=600, poll_s=30):
    """Block until MemAvailable can hold a mesh of this size (peak RSS ~0.85 kB per tet seen here),
    so two queues never run two large Netgen jobs into the OOM killer."""
    need = target_elements * mb_per_element + reserve_mb
    while True:
        avail = next(int(l.split()[1]) for l in open("/proc/meminfo") if l.startswith("MemAvailable")) / 1024
        if avail >= need:
            return
        log.info("waiting for memory: %.0f MB available, %.0f MB wanted", avail, need)
        time.sleep(poll_s)


def edge_length_for(target_elements, calib_h, calib_n):
    return calib_h * (calib_n / target_elements) ** (1.0 / 3.0)


def cap_areas(surface, faces, face_id_array):
    summ = meshlib.face_summary(surface, face_id_array)
    return {faces[str(fid)]: v["area"] for fid, v in summ.items() if faces[str(fid)].startswith("cap_")}, summ


def make_conditions(spec, areas):
    bc = spec["bc"]
    inflow_spec = bc["inflow"]
    conds = {}
    inflow_caps = inflow_spec["cap"] if isinstance(inflow_spec["cap"], list) else [inflow_spec["cap"]]
    period = inflow_spec.get("period_s", 1.0)
    for i, cap in enumerate(inflow_caps):
        q = inflow_spec["mean_flow_cm3s"][i] if isinstance(inflow_spec["mean_flow_cm3s"], list) else inflow_spec["mean_flow_cm3s"]
        kind = inflow_spec.get("waveform", "aortic")
        if kind == "aortic":
            conds[cap] = bclib.aortic_waveform(q, period)
        elif kind == "venous":
            conds[cap] = bclib.venous_waveform(q, period, pulsatility=inflow_spec.get("pulsatility", 0.15),
                                               phase=inflow_spec.get("phase", [0.0] * len(inflow_caps))[i] if isinstance(inflow_spec.get("phase"), list) else 0.0)
        elif kind == "steady":
            conds[cap] = Inflow.steady(q, period)
        elif kind == "file":
            from svromsetup.boundary_conditions import read_flow_file
            wf = read_flow_file(inflow_spec["file"] if not isinstance(inflow_spec["file"], list) else inflow_spec["file"][i])
            if wf.mean < 0:          # svSolver files are negative into the model
                wf = wf.flipped()
            if q:                    # optional rescale to a stated mean
                wf = Inflow(wf.time, tuple(float(x) * q / wf.mean for x in wf.flow), source=wf.source)
            conds[cap] = wf
        else:
            raise ValueError(kind)
    out_spec = bc["outlets"]
    outlet_areas = {k: v for k, v in areas.items() if k not in inflow_caps}
    total_inflow = sum(conds[c].mean for c in inflow_caps)
    if out_spec["type"] == "RCR":
        rcr, frac = bclib.rcr_from_targets(outlet_areas, total_inflow, out_spec["mean_pressure_mmHg"],
                                           flow_split=out_spec.get("flow_split"),
                                           rp_fraction=out_spec.get("Rp_fraction", 0.06),
                                           total_compliance=out_spec.get("C_total_cm5_per_dyn", 1.0e-3),
                                           distal_pressure_mmHg=out_spec.get("distal_pressure_mmHg", 0.0),
                                           area_exponent=out_spec.get("area_exponent", 1.0))
        conds.update(rcr)
    elif out_spec["type"] == "Resistance":
        res, frac = bclib.resistance_from_targets(outlet_areas, total_inflow,
                                                  pvr_wood_units=out_spec.get("pvr_wood_units"),
                                                  mean_pressure_mmHg=out_spec.get("mean_pressure_mmHg"),
                                                  distal_pressure_mmHg=out_spec.get("distal_pressure_mmHg", 0.0),
                                                  flow_split=out_spec.get("flow_split"),
                                                  area_exponent=out_spec.get("area_exponent", 1.0))
        conds.update(res)
    elif out_spec["type"] == "from_project":
        frac = {}
        for cap, v in out_spec["values_by_cap"].items():
            t = (v.get("type") or "RCR").lower(); vals = v["values"]
            if t == "rcr" and len(vals) >= 3:
                conds[cap] = RCR(Rp=vals[0], C=vals[1], Rd=vals[2], Pd=v.get("distal_pressure", 0.0))
            elif t == "resistance" and vals:
                conds[cap] = Resistance(R=vals[0])
            else:
                raise ValueError(f"{cap}: unsupported project BC {v}")
        missing = [c for c in outlet_areas if c not in conds]
        if missing:
            raise ValueError(f"project gives no BC for outlets {missing}")
    else:
        raise ValueError(out_spec["type"])
    return conds, frac


def conditions_json(conds, frac):
    out = {}
    for k, c in conds.items():
        if isinstance(c, Inflow):
            out[k] = {"type": "inflow", "mean_cm3s": c.mean, "period_s": c.period, "waveform": c.source,
                      "time": list(c.time), "flow_into_domain": list(c.flow)}
        elif isinstance(c, RCR):
            out[k] = {"type": "RCR", "Rp": c.Rp, "C": c.C, "Rd": c.Rd, "Pd": c.Pd, "flow_fraction": frac.get(k),
                      "units": "dyn s/cm^5, cm^5/dyn, dyn/cm^2"}
        elif isinstance(c, Resistance):
            out[k] = {"type": "Resistance", "R": c.R, "flow_fraction": frac.get(k), "units": "dyn s/cm^5"}
    return out


def solver_values(spec, period):
    s = spec.get("solver", {})
    dt = s.get("time_step_s", 1.0e-3)
    return {
        "mesh.scale_factor": 1.0,                      # meshes are already in cm
        "general.time_step_size": dt,
        "general.number_of_time_steps": int(round(s.get("cycles", 1) * period / dt)),
        "general.increment_in_saving_vtk_files": s.get("save_every", 50),
        "general.increment_in_saving_restart_files": s.get("save_every", 50),
        "fluid.density": s.get("density", 1.06),
        "fluid.viscosity": s.get("viscosity", 0.04),
        "fluid.max_iterations": s.get("newton_max_iterations", 5),
        "fluid.tolerance": s.get("newton_tolerance", 1e-4),
        "ls.max_iterations": s.get("ls_max_iterations", 10),
        "ls.ns_gm_max_iterations": s.get("ns_gm_max_iterations", 200),
        "ls.ns_cg_max_iterations": s.get("ns_cg_max_iterations", 500),
        "ls.tolerance": s.get("ls_tolerance", 0.4),
        "ls.ns_gm_tolerance": s.get("ns_gm_tolerance", 0.01),
        "ls.ns_cg_tolerance": s.get("ns_cg_tolerance", 0.2),
        "bc.rcr_initial_pressure": s.get("rcr_initial_pressure", 0.0),
    }


def build(spec, res_names, out_root, calibrate_only=False):
    name = spec["name"]
    surface = meshlib.read_polydata(spec["surface"])
    fid_arr = spec.get("face_id_array", "ModelFaceID")
    cap_arr = spec.get("cap_id_array", "CapID")
    faces = spec["faces"]
    areas, summ = cap_areas(surface, faces, fid_arr)
    conds, frac = make_conditions(spec, areas)
    period = max(c.period for c in conds.values() if isinstance(c, Inflow))
    values = solver_values(spec, period)
    names = {int(k): v for k, v in faces.items()}
    log.info("%s: caps %s", name, {k: round(v, 3) for k, v in areas.items()})

    calib = spec.get("calibration")
    if calib is None or calibrate_only:
        h0 = spec.get("calibration_h", 0.2)
        log.info("calibrating at h=%g", h0)
        mesh, dt = meshlib.generate(surface, h0, cell_entity_ids=fid_arr, cap_ids=cap_arr,
                                    mesher=spec.get("mesher", "netgen"), skip_capping=spec.get("surface_is_capped", True))
        n = meshlib.count_cells(mesh)["tet4"]
        calib = {"h": h0, "n_tets": n, "seconds": dt}
        spec["calibration"] = calib
        log.info("calibration: h=%g -> %d tets in %.0fs", h0, n, dt)
        if calibrate_only:
            return spec

    # (h, n_tets) observations: element count does not scale as a clean h^-3 with Netgen's
    # grading, so the edge length for a target is taken from the two nearest observations.
    obs = [(calib["h"], calib["n_tets"])] + [tuple(o) for o in spec.get("observations", [])]

    def h_for(target):
        """Edge length for a target count from a log-log least-squares fit of every (h, n)
        observation, exponent bounded to [2.5, 5.5] (n ~ h^-3 in 3D, Netgen grading bends it).
        Anchored at the observation nearest the target, so the fit only sets the slope."""
        hs = np.array([o[0] for o in obs], float); ns = np.array([o[1] for o in obs], float)
        expo = 3.0
        if len(set(np.round(hs, 5))) >= 2:
            A = np.vstack([np.log(hs), np.ones_like(hs)]).T
            slope = np.linalg.lstsq(A, np.log(ns), rcond=None)[0][0]
            expo = float(min(max(-slope, 2.5), 5.5))
        i = int(np.argmin(np.abs(np.log(ns / target))))
        return float(hs[i] * (ns[i] / target) ** (1.0 / expo)), expo

    results = {}
    for res in res_names:
        target = spec["resolutions"][res]
        out_dir = os.path.join(out_root, name, res)
        if os.path.exists(os.path.join(out_dir, "case.json")):
            log.info("%s/%s exists, skipping", name, res)
            continue
        os.makedirs(out_dir, exist_ok=True)
        lock = None
        if target >= 4.0e6:   # one large job at a time across queues
            import fcntl
            lock = open("/tmp/cvbench_large.lock", "w"); fcntl.flock(lock, fcntl.LOCK_EX)
        wait_for_memory(target)
        t0 = time.time()
        best = None
        for attempt in range(spec.get("max_attempts", 2)):
            h, expo = h_for(target)
            if any(abs(h / o[0] - 1.0) < 0.01 for o in obs):
                h *= 0.97  # do not repeat an edge length already tried (Netgen's size field is quantized)
            log.info("== %s/%s: target %.3g elements -> h = %.4f %s (exponent %.2f, attempt %d)", name, res, target, h,
                     spec.get("units", "cm"), expo, attempt + 1)
            mesh_try, t_try = meshlib.generate(surface, h, cell_entity_ids=fid_arr, cap_ids=cap_arr,
                                               mesher=spec.get("mesher", "netgen"), skip_capping=spec.get("surface_is_capped", True))
            n = meshlib.count_cells(mesh_try)["tet4"]
            obs.append((h, n))
            spec.setdefault("observations", []).append([h, n])
            if best is None or abs(math.log(n / target)) < abs(math.log(best[2] / target)):
                # park the best attempt on disk, not in memory: two 5 M meshes do not fit in 7 GB
                tmp = os.path.join(out_dir, "_best_attempt.vtu")
                meshlib.write_dataset(mesh_try, tmp)
                best = (tmp, h, n, t_try)
            del mesh_try
            import gc; gc.collect()
            if abs(n / target - 1.0) <= spec.get("tolerance", 0.25):
                break
            log.info("   %d tets is outside +/-%.0f%% of target %d; re-estimating h", n, 100 * spec.get("tolerance", 0.25), target)
        tmp, h, n, t_mesh = best
        r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(tmp); r.Update(); mesh = r.GetOutput()
        os.remove(tmp)
        log.info("   keeping h = %.4f with %d tets (closest to target %d)", h, n, target)
        counts = meshlib.count_cells(mesh)
        log.info("   meshed: %s in %.0fs", counts, t_mesh)
        # svMultiPhysics case (mesh-complete + solver.xml + .flow) via SimVascular Mesh Prep / MultiPhysics packages
        t1 = time.time()
        case = write_case(out_dir, mesh, names, conds, values, face_id_array_name=fid_arr)
        t_case = time.time() - t1
        # MFEM mesh
        t2 = time.time()
        mfem_dir = os.path.join(out_dir, "mfem"); os.makedirs(mfem_dir, exist_ok=True)
        attrs = to_mfem.write_mfem_gmsh(mesh, os.path.join(mfem_dir, f"{name}_{res}.msh"), face_id_array=fid_arr,
                                   names_by_face_id=names)
        t_mfem = time.time() - t2
        # quality
        t3 = time.time()
        rep = quality.quality_report(mesh, spec.get("units", "cm"))
        json.dump(rep, open(os.path.join(out_dir, "quality.json"), "w"), indent=2)
        open(os.path.join(out_dir, "quality.md"), "w").write(quality.report_markdown(rep, f"{name} / {res}") + "\n")
        t_q = time.time() - t3
        meta = {
            "name": name, "resolution": res, "target_elements": target, "edge_length": h,
            "units": spec.get("units", "cm"), "mesher": spec.get("mesher", "netgen"),
            "surface_source": spec.get("source"), "license": spec.get("license"),
            "elements": counts, "faces": names, "cap_areas_cm2": areas,
            "boundary_conditions": conditions_json(conds, frac), "solver_values": values,
            "timings_s": {"mesh": t_mesh, "mesh_complete_and_xml": t_case, "mfem": t_mfem, "quality": t_q,
                          "total": time.time() - t0},
            "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
            "mfem_attributes": attrs,
            "pipeline": "SlicerExtension-VMTK CfdMeshGenerator (MeshingPipeline, Netgen) -> SlicerSimVascular svmeshcomplete/svmpsetup -> to_mfem/to_gmsh",
        }
        json.dump(meta, open(os.path.join(out_dir, "case.json"), "w"), indent=2)
        results[res] = meta
        log.info("   done %s/%s in %.0fs (peak RSS %.0f MB)", name, res, meta["timings_s"]["total"], meta["peak_rss_mb"])
        if lock:
            lock.close()
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--res", default="coarse,medium,fine")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cases"))
    ap.add_argument("--calibrate-only", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    spec = json.load(open(a.spec))
    r = build(spec, a.res.split(","), a.out, a.calibrate_only)
    json.dump(spec, open(a.spec, "w"), indent=2)   # keeps calibration + (h, n) observations
