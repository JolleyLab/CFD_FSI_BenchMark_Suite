#!/usr/bin/env python3
"""Fetch Vascular Model Repository (VMR) SimVascular projects and turn each into a benchmark
case spec for build_case.py, reusing the project's own face names and boundary conditions.

    vmr_fetch.py 0031_H_ABAO_AAA 0054_H_PULMFON_TAT 0065_H_PULMFON_SVD 0007_H_AO_H
    vmr_fetch.py --from-zip downloads/0031_H_ABAO_AAA.zip      # if downloaded by hand

Project zips are https://www.vascularmodel.com/svprojects/<NAME>.zip (100-300 MB each; the
catalogue is dataset/dataset-svprojects.csv in github.com/SimVascular/vascularmodel). What is used:
  Models/<model>.vtp   SimVascular PolyData model: closed surface with ModelFaceID (+ CapID)
  Models/<model>.mdl   XML: <face id=.. name=.. type="wall|cap">  -> face names
  Simulations/<job>.sjb XML: per-cap "BC Type" / "Values" / flow file -> inflow + RCR/resistance
  flow-files/*.flow, Simulations/<job>/rcrt.dat                -> fallbacks for the above
Licence: VMR README-COPYRIGHT (research and development use, keep the notice, acknowledge).

Credits: data from the Vascular Model Repository (Wilson et al. 2013, doi:10.1115/1.4025983), VMR licence.
See REFERENCES.md for the per-model dataset DOIs and source studies.
"""
import argparse, glob, json, os, re, sys, zipfile, urllib.request, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

VMR_URL = "https://www.vascularmodel.com/svprojects/{name}.zip"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

DEFAULT_RESOLUTIONS = {"coarse": 1500000, "medium": 3000000, "fine": 5000000}


def download(name, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, f"{name}.zip")
    if not os.path.exists(dest):
        url = VMR_URL.format(name=name)
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, dest)
    return dest


def extract(zip_path, dest_dir):
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(dest_dir)
    # the project is the folder holding a .svproj file
    for root, dirs, files in os.walk(dest_dir):
        if any(f.endswith(".svproj") for f in files):
            return root
    return dest_dir


def _parse_sv_xml(path):
    """SimVascular files hold two root elements (<format/> then the document); wrap them."""
    text = open(path, encoding="utf-8", errors="replace").read()
    text = re.sub(r"^\s*<\?xml[^>]*\?>", "", text)
    return ET.fromstring("<sv_root>" + text + "</sv_root>")


def read_mdl(path):
    """{face id: (name, type)} from a SimVascular .mdl."""
    faces = {}
    tree = _parse_sv_xml(path)
    for face in tree.iter("face"):
        fid = int(face.get("id")); faces[fid] = (face.get("name"), face.get("type", ""))
    return faces


def read_sjb(path):
    """Per-cap boundary condition properties and basic props from a SimVascular .sjb job."""
    tree = _parse_sv_xml(path)
    job = {"basic": {}, "caps": {}, "wall": {}, "solver": {}}
    for sec, key in (("basic_props", "basic"), ("wall_props", "wall"), ("solver_props", "solver")):
        for el in tree.iter(sec):
            for prop in el.iter("prop"):
                job[key][prop.get("key")] = prop.get("value")
    for el in tree.iter("cap_props"):
        for cap in el.iter("cap"):
            job["caps"][cap.get("name")] = {p.get("key"): p.get("value") for p in cap.iter("prop")}
    return job


def read_legacy_rcrt(path):
    """svSolver rcrt.dat without face names: list of (Rp, C, Rd, Pd) in outlet order."""
    toks = [l.strip() for l in open(path).read().splitlines() if l.strip()]
    n = int(toks[0]); out = []; i = 1
    while i < len(toks):
        npts = int(toks[i]); rp, c, rd = float(toks[i + 1]), float(toks[i + 2]), float(toks[i + 3])
        pd = float(toks[i + 4].split()[1]) if npts else 0.0
        out.append((rp, c, rd, pd)); i += 4 + npts
    return out


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")


def spec_from_project(project_dir, model_id, out_surfaces, resolutions=None):
    mdls = glob.glob(os.path.join(project_dir, "Models", "*.mdl"))
    vtps = glob.glob(os.path.join(project_dir, "Models", "*.vtp"))
    if not vtps:
        raise FileNotFoundError(f"no Models/*.vtp in {project_dir}")
    # prefer the vtp sharing a stem with the mdl
    vtp = vtps[0]
    faces_mdl = {}
    if mdls:
        faces_mdl = read_mdl(mdls[0])
        stem = os.path.splitext(os.path.basename(mdls[0]))[0]
        for v in vtps:
            if os.path.splitext(os.path.basename(v))[0] == stem:
                vtp = v
    import meshlib
    surface = meshlib.read_polydata(vtp)
    cd = surface.GetCellData()
    if cd.GetArray("ModelFaceID") is None:
        raise ValueError(f"{vtp}: no ModelFaceID cell array")
    summ = meshlib.face_summary(surface, "ModelFaceID")
    faces = {}
    for fid, info in summ.items():
        nm, typ = faces_mdl.get(fid, (f"face_{fid}", ""))
        nm = sanitize(nm)
        if typ == "cap" or nm.lower().startswith("cap") or (not typ and cd.GetArray("CapID") is not None
                                                            and _cap_by_capid(surface, fid)):
            faces[str(fid)] = nm if nm.startswith("cap_") else f"cap_{nm}"
        else:
            faces[str(fid)] = nm if (nm == "wall" or nm.startswith("wall_")) else f"wall_{nm}"
    # boundary conditions from the job file
    sjbs = glob.glob(os.path.join(project_dir, "Simulations", "*.sjb"))
    bc = {"inflow": None, "outlets": {"type": "RCR", "values_by_cap": {}}}
    notes = []
    if sjbs:
        job = read_sjb(sjbs[0])
        notes.append(f"BCs read from {os.path.basename(sjbs[0])}")
        for cap, props in job["caps"].items():
            bctype = (props.get("BC Type") or "").lower()
            sc = sanitize(cap)
            cname = next((v for v in faces.values() if v in (f"cap_{sc}", sc)), None) or \
                next((v for v in faces.values() if v.endswith("_" + sc)), f"cap_{sc}")
            if "velocit" in bctype or "inflow" in bctype or props.get("Original File"):
                flow_file = props.get("Original File") or ""
                cands = glob.glob(os.path.join(project_dir, "**", flow_file), recursive=True) if flow_file else []
                bc["inflow"] = {"cap": cname, "waveform": "file", "file": cands[0] if cands else None,
                                "period_s": float(props.get("Period", 1.0) or 1.0), "sjb_props": props}
            elif bctype in ("rcr", "resistance", "coronary"):
                vals = [float(x) for x in (props.get("Values") or "").split()]
                bc["outlets"]["values_by_cap"][cname] = {"type": props.get("BC Type"), "values": vals,
                                                        "distal_pressure": float(props.get("Pressure", 0) or 0)}
        if job["basic"]:
            bc["fluid"] = {"density": job["basic"].get("Fluid Density"), "viscosity": job["basic"].get("Fluid Viscosity")}
    rcrts = glob.glob(os.path.join(project_dir, "Simulations", "**", "rcrt.dat"), recursive=True)
    if rcrts:
        bc["outlets"]["rcrt_dat"] = rcrts[0]
        try:
            bc["outlets"]["rcrt_legacy_blocks"] = read_legacy_rcrt(rcrts[0])
        except Exception as e:  # noqa
            notes.append(f"rcrt.dat not parsed: {e}")
    flows = glob.glob(os.path.join(project_dir, "flow-files", "*.flow")) + glob.glob(os.path.join(project_dir, "Simulations", "**", "*.flow"), recursive=True)
    if bc["inflow"] is None and flows:
        bc["inflow"] = {"cap": None, "waveform": "file", "file": flows[0], "period_s": None}
        notes.append("inflow cap not identified from .sjb; assign bc.inflow.cap by hand (largest CapID==1 face)")
    os.makedirs(out_surfaces, exist_ok=True)
    dest_vtp = os.path.join(out_surfaces, f"{model_id}_model.vtp")
    import shutil; shutil.copy(vtp, dest_vtp)
    spec = {
        "name": model_id,
        "description": f"VMR model {model_id} (SimVascular project), rigid-wall CFD at benchmark resolutions",
        "surface": os.path.relpath(dest_vtp, ROOT),
        "source": {"repository": "Vascular Model Repository", "model": model_id,
                   "project_zip": VMR_URL.format(name=model_id), "model_vtp": os.path.relpath(vtp, project_dir),
                   "mdl": os.path.relpath(mdls[0], project_dir) if mdls else None},
        "license": "VMR README-COPYRIGHT (research and development use; keep notice; acknowledge NLM R01LM013120 / NHLBI HHSN268201100035C and/or cite Wilson et al. 2013 doi:10.1115/1.4025983)",
        "units": "cm", "face_id_array": "ModelFaceID",
        "cap_id_array": "CapID" if cd.GetArray("CapID") is not None else None,
        "surface_is_capped": True,
        "faces": faces,
        "vmr_boundary_conditions": bc,
        "bc": _bc_for_build(bc, faces, summ),
        "notes": notes,
        "solver": {"time_step_s": 1.0e-3, "cycles": 1, "save_every": 50, "density": 1.06, "viscosity": 0.04},
        "mesher": "netgen",
        "calibration_h": float(0.045 * np.sqrt(sum(v["area"] for v in summ.values()) / len(summ))),
        "resolutions": resolutions or DEFAULT_RESOLUTIONS,
    }
    return spec


def _cap_by_capid(surface, fid):
    from vtk.util.numpy_support import vtk_to_numpy
    f = vtk_to_numpy(surface.GetCellData().GetArray("ModelFaceID")); c = vtk_to_numpy(surface.GetCellData().GetArray("CapID"))
    return bool((c[f == fid] > 0).all())


def _bc_for_build(vmr_bc, faces, summ):
    """Translate the project's conditions into build_case's spec form where possible; the raw
    values are kept beside it in `vmr_boundary_conditions` either way."""
    caps = [v for v in faces.values() if v.startswith("cap_")]
    inflow = vmr_bc.get("inflow") or {}
    inflow_cap = inflow.get("cap") or (caps[0] if caps else None)
    if inflow_cap not in caps:   # .sjb calls it "inflow" but the model face carries a longer name
        inflow_cap = next((c for c in caps if c.lower().startswith("cap_inflow")), inflow_cap)
    out = {"inflow": {"cap": inflow_cap, "waveform": "file" if inflow.get("file") else "aortic",
                      "file": inflow.get("file"), "period_s": inflow.get("period_s") or 1.0, "mean_flow_cm3s": None},
           "outlets": {"type": "RCR", "mean_pressure_mmHg": 90.0, "Rp_fraction": 0.06, "C_total_cm5_per_dyn": 1.0e-3}}
    vals = vmr_bc.get("outlets", {}).get("values_by_cap", {})
    if vals:
        out["outlets"] = {"type": "from_project", "values_by_cap": vals}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="*")
    ap.add_argument("--from-zip", nargs="*", default=[])
    ap.add_argument("--downloads", default=os.path.join(ROOT, "downloads"))
    ap.add_argument("--specs", default=os.path.join(ROOT, "specs"))
    ap.add_argument("--surfaces", default=os.path.join(ROOT, "surfaces"))
    a = ap.parse_args()
    zips = list(a.from_zip) + [download(m, a.downloads) for m in a.models]
    for z in zips:
        mid = os.path.splitext(os.path.basename(z))[0]
        proj = extract(z, os.path.join(a.downloads, mid))
        spec = spec_from_project(proj, mid, a.surfaces)
        path = os.path.join(a.specs, f"{mid}.json")
        json.dump(spec, open(path, "w"), indent=2)
        print(f"wrote {path}: faces {spec['faces']}; notes {spec['notes']}")
        print("  review bc.inflow.cap / bc.outlets, then: python tools/build_case.py", path)
