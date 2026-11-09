#!/usr/bin/env python3
"""Structured hexahedral O-grid meshes of a straight pipe: fluid only (CFD) and fluid + wall
(ALE FSI), at three resolutions, written as svMultiPhysics mesh-complete folders, solver.xml and
MFEM meshes (fluid attribute 1, solid attribute 2). CGS units.

Geometry: radius R = 1.0 cm, length L = 11 cm (same as the tet `pipe` case), wall thickness
t = 0.1 cm (R/10) in n_t layers. Resolution is set by the cell size h; cells are near-cubic.

Credits: mesh-complete writer from SlicerSimVascular svmeshcomplete; FSI input modelled on
svMultiPhysics tests/cases/fsi/pipe_RCR_3d. See REFERENCES.md.
"""
import json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

import bclib, quality, to_mfem
from svmeshcomplete.mesh_complete import write_mesh_complete
from svmeshcomplete.face_table import FaceTable, Face
from svmpsetup.case import write_case
from svromsetup.boundary_conditions import Inflow, RCR

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def ogrid_section(R, n_core, n_rad, core_frac=0.5, R_outer=None, n_t=0):
    """Points and quads of an O-grid disc (and an optional annulus of n_t layers out to R_outer).
    Returns (points[N,2], quads[M,4], ring_index_of_outer_boundary_points, annulus_quads)."""
    # core: square of half-width a = core_frac*R mapped slightly toward a circle
    a = core_frac * R
    pts = {}
    def pid(key, xy):
        if key not in pts:
            pts[key] = (len(pts), xy)
        return pts[key][0]
    quads = []
    # core square grid
    u = np.linspace(-1, 1, n_core + 1)
    def core_pt(i, j):
        x, y = a * u[i], a * u[j]
        # bulge the square so the outer ring blocks stay regular
        s = 0.15
        return (x * (1 - s * (y / a) ** 2) if False else x, y) if False else (x + s * x * (1 - (y / a) ** 2) * 0, y)
    core_ids = np.zeros((n_core + 1, n_core + 1), int)
    for i in range(n_core + 1):
        for j in range(n_core + 1):
            core_ids[i, j] = pid(("c", i, j), (a * u[i], a * u[j]))
    for i in range(n_core):
        for j in range(n_core):
            quads.append((core_ids[i, j], core_ids[i + 1, j], core_ids[i + 1, j + 1], core_ids[i, j + 1]))
    # four outer blocks: from each square edge radially to the circle
    # parametrize boundary of square by t in [0, 4): edge k, local s in [0,1]
    def square_pt(k, s):
        if k == 0: return (a, -a + 2 * a * s)
        if k == 1: return (a - 2 * a * s, a)
        if k == 2: return (-a, a - 2 * a * s)
        return (-a + 2 * a * s, -a)
    def square_id(k, i):
        if k == 0: return core_ids[n_core, i]
        if k == 1: return core_ids[n_core - i, n_core]
        if k == 2: return core_ids[0, n_core - i]
        return core_ids[i, 0]
    ring_ids = {}   # (k,i,r) -> id ; r = 0 square ... n_rad circle
    outer_boundary = []
    for k in range(4):
        for i in range(n_core + 1):
            s = i / n_core
            sx, sy = square_pt(k, s)
            ang = math.atan2(sy, sx)
            cx, cy = R * math.cos(ang), R * math.sin(ang)
            for r in range(n_rad + 1):
                key = ("r", k, i, r)
                # share edges between blocks
                if i == n_core:
                    key = ("r", (k + 1) % 4, 0, r)
                if r == 0:
                    ring_ids[(k, i, r)] = square_id(k, i)
                    continue
                w = r / n_rad
                # blend linearly between square point and circle point (radial lines)
                x = (1 - w) * sx + w * cx
                y = (1 - w) * sy + w * cy
                ring_ids[(k, i, r)] = pid(key, (x, y))
    for k in range(4):
        for i in range(n_core):
            for r in range(n_rad):
                quads.append((ring_ids[(k, i, r)], ring_ids[(k, i, r + 1)], ring_ids[(k, i + 1, r + 1)], ring_ids[(k, i + 1, r)]))
    outer = sorted({ring_ids[(k, i, n_rad)] for k in range(4) for i in range(n_core + 1)})
    ann_quads = []
    if n_t and R_outer:
        # annulus layers from R to R_outer, following the same angular sampling as the circle ring
        ann_ids = {}
        for k in range(4):
            for i in range(n_core + 1):
                base = ring_ids[(k, i, n_rad)]
                ang = math.atan2(pts_xy(pts, base)[1], pts_xy(pts, base)[0])
                for l in range(1, n_t + 1):
                    rr = R + (R_outer - R) * l / n_t
                    key = ("a", (k + 1) % 4 if i == n_core else k, 0 if i == n_core else i, l)
                    ann_ids[(k, i, l)] = pid(key, (rr * math.cos(ang), rr * math.sin(ang)))
                ann_ids[(k, i, 0)] = base
        for k in range(4):
            for i in range(n_core):
                for l in range(n_t):
                    ann_quads.append((ann_ids[(k, i, l)], ann_ids[(k, i, l + 1)], ann_ids[(k, i + 1, l + 1)], ann_ids[(k, i + 1, l)]))
    P = np.zeros((len(pts), 2))
    for key, (idx, xy) in pts.items():
        P[idx] = xy
    return P, np.array(quads), outer, np.array(ann_quads) if len(ann_quads) else np.zeros((0, 4), int)


def pts_xy(pts, idx):
    for key, (i, xy) in pts.items():
        if i == idx:
            return xy
    raise KeyError(idx)


def fix_orientation(P, quads):
    """Counter-clockwise quads so extruded hexes have positive Jacobian."""
    q = quads.copy()
    for n, (a, b, c, d) in enumerate(q):
        area = 0.0
        poly = [P[a], P[b], P[c], P[d]]
        for i in range(4):
            x0, y0 = poly[i]; x1, y1 = poly[(i + 1) % 4]
            area += x0 * y1 - x1 * y0
        if area < 0:
            q[n] = (a, d, c, b)
    return q


def extrude(P2, quads, z, region_id):
    """Extrude a 2D quad mesh along z -> hex8 (VTK ordering: bottom quad then top quad)."""
    nz = len(z); n2 = len(P2)
    pts = np.zeros((n2 * nz, 3))
    pts[:, :2] = np.tile(P2, (nz, 1)); pts[:, 2] = np.repeat(z, n2)
    q = np.asarray(quads, np.int64)
    k = np.arange(nz - 1)[:, None, None]
    bottom = q[None, :, :] + k * n2
    top = bottom + n2
    hexes = np.concatenate([bottom, top], axis=2).reshape(-1, 8)
    return pts, hexes, n2


def boundary_quads(P2, quads, z, n2, outer_ids, inner_ids=None):
    """Boundary quads with face ids: 1 = outer lateral surface, 2 = z0 cap, 3 = z1 cap,
    4 = inner lateral surface (annulus only)."""
    nz = len(z)
    out = []
    # caps
    for (a, b, c, d) in quads:
        out.append((2, (a, d, c, b)))                                           # z = z0, normal -z
        o = (nz - 1) * n2
        out.append((3, (o + a, o + b, o + c, o + d)))                           # z = z1, normal +z
    # lateral: consecutive outer boundary points ordered by angle
    def lateral(ids, fid, outward=True):
        ang = np.arctan2(P2[ids, 1], P2[ids, 0]); order = [ids[i] for i in np.argsort(ang)]
        for k in range(nz - 1):
            o0, o1 = k * n2, (k + 1) * n2
            for i in range(len(order)):
                a, b = order[i], order[(i + 1) % len(order)]
                q = (o0 + a, o0 + b, o1 + b, o1 + a) if outward else (o0 + b, o0 + a, o1 + a, o1 + b)
                out.append((fid, q))
    lateral(outer_ids, 1, True)
    if inner_ids is not None:
        lateral(inner_ids, 4, False)
    return out


def to_ug(points, hexes, bnd, region_id, face_array="ModelFaceID"):
    """Hexes first, then boundary quads, as one vtkUnstructuredGrid (vectorized)."""
    hexes = np.asarray(hexes, np.int64)
    quads = np.array([q for _, q in bnd], np.int64).reshape(-1, 4)
    fids_b = np.array([f for f, _ in bnd], np.int32)
    n_h, n_q = len(hexes), len(quads)
    conn = np.concatenate([hexes.ravel(), quads.ravel()])
    offsets = np.concatenate([np.arange(0, 8 * n_h + 1, 8), 8 * n_h + np.arange(4, 4 * n_q + 1, 4)]) if n_q else np.arange(0, 8 * n_h + 1, 8)
    ca = vtk.vtkCellArray()
    ca.SetData(numpy_to_vtk(offsets.astype(np.int64), deep=True, array_type=vtk.VTK_ID_TYPE),
               numpy_to_vtk(conn.astype(np.int64), deep=True, array_type=vtk.VTK_ID_TYPE))
    types = np.concatenate([np.full(n_h, vtk.VTK_HEXAHEDRON, np.uint8), np.full(n_q, vtk.VTK_QUAD, np.uint8)])
    ug = vtk.vtkUnstructuredGrid()
    vp = vtk.vtkPoints(); vp.SetData(numpy_to_vtk(np.ascontiguousarray(points), deep=True)); ug.SetPoints(vp)
    ug.SetCells(numpy_to_vtk(types, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR), ca)
    a = numpy_to_vtk(np.concatenate([np.zeros(n_h, np.int32), fids_b]), deep=True); a.SetName(face_array); ug.GetCellData().AddArray(a)
    r = numpy_to_vtk(np.full(n_h + n_q, region_id, np.int32), deep=True); r.SetName("ModelRegionID"); ug.GetCellData().AddArray(r)
    return ug


def build_pipe(h, R=1.0, L=11.0, t=0.1, n_t=3, with_wall=True):
    n_core = max(2, int(round(R / h)))
    n_rad = max(1, int(round(0.5 * R / h)))
    nz = max(2, int(round(L / h)) + 1)
    z = np.linspace(0.5, 0.5 + L, nz)   # same axial placement as the tet pipe (clipped at 0.5 and 11.5)
    P2, quads, outer, ann = ogrid_section(R, n_core, n_rad, R_outer=R + t if with_wall else None, n_t=n_t if with_wall else 0)
    quads = fix_orientation(P2, quads)
    pts, hexes, n2 = extrude(P2, quads, z, 1)
    bnd = boundary_quads(P2, quads, z, n2, outer)
    fluid = to_ug(pts, hexes, bnd, 1)
    solid = None
    if with_wall:
        ann = fix_orientation(P2, ann)
        # the annulus uses the same 2D point array; its own boundary: outer ring at R+t, inner ring at R
        ann_outer = sorted({v for q in ann for v in q if abs(np.hypot(*P2[v]) - (R + t)) < 1e-9})
        pts_s, hex_s, _ = extrude(P2, ann, z, 2)
        bnd_s = boundary_quads(P2, ann, z, n2, ann_outer, inner_ids=outer)
        solid = to_ug(pts_s, hex_s, bnd_s, 2)
        solid = compact(solid)
    fluid = compact(fluid)
    return fluid, solid, dict(n_core=n_core, n_rad=n_rad, nz=nz, h=h)


def compact(ug):
    """Drop unused points (the 2D array is shared between disc and annulus)."""
    cl = vtk.vtkStaticCleanUnstructuredGrid() if hasattr(vtk, "vtkStaticCleanUnstructuredGrid") else None
    if cl is None:
        return ug
    cl.SetInputData(ug); cl.RemoveUnusedPointsOn(); cl.ToleranceIsAbsoluteOn(); cl.SetAbsoluteTolerance(0.0); cl.Update()
    return cl.GetOutput()


FSI_XML = """<?xml version="1.0" encoding="UTF-8" ?>
<!-- cvbench pipe_fsi_hex: ALE FSI, fluid hex O-grid + hex wall, svMultiPhysics -->
<svMultiPhysicsFile version="0.1">
<GeneralSimulationParameters>
  <Continue_previous_simulation> 0 </Continue_previous_simulation>
  <Number_of_spatial_dimensions> 3 </Number_of_spatial_dimensions>
  <Number_of_time_steps> {nsteps} </Number_of_time_steps>
  <Time_step_size> {dt} </Time_step_size>
  <Spectral_radius_of_infinite_time_step> 0.50 </Spectral_radius_of_infinite_time_step>
  <Searched_file_name_to_trigger_stop> STOP_SIM </Searched_file_name_to_trigger_stop>
  <Save_results_to_VTK_format> true </Save_results_to_VTK_format>
  <Name_prefix_of_saved_VTK_files> result </Name_prefix_of_saved_VTK_files>
  <Increment_in_saving_VTK_files> 50 </Increment_in_saving_VTK_files>
  <Start_saving_after_time_step> 1 </Start_saving_after_time_step>
  <Increment_in_saving_restart_files> 50 </Increment_in_saving_restart_files>
  <Convert_BIN_to_VTK_format> 0 </Convert_BIN_to_VTK_format>
  <Verbose> 1 </Verbose> <Warning> 0 </Warning> <Debug> 0 </Debug>
</GeneralSimulationParameters>
<Add_mesh name="lumen" >
  <Mesh_file_path> mesh/fluid/mesh-complete.mesh.vtu </Mesh_file_path>
  <Add_face name="lumen_inlet">  <Face_file_path> mesh/fluid/mesh-surfaces/cap_inlet.vtp </Face_file_path> </Add_face>
  <Add_face name="lumen_outlet"> <Face_file_path> mesh/fluid/mesh-surfaces/cap_outlet.vtp </Face_file_path> </Add_face>
  <Add_face name="lumen_wall">   <Face_file_path> mesh/fluid/mesh-surfaces/wall.vtp </Face_file_path> </Add_face>
  <Domain> 0 </Domain>
</Add_mesh>
<Add_mesh name="wall" >
  <Mesh_file_path> mesh/solid/mesh-complete.mesh.vtu </Mesh_file_path>
  <Add_face name="wall_inlet">  <Face_file_path> mesh/solid/mesh-surfaces/cap_inlet_ring.vtp </Face_file_path> </Add_face>
  <Add_face name="wall_outlet"> <Face_file_path> mesh/solid/mesh-surfaces/cap_outlet_ring.vtp </Face_file_path> </Add_face>
  <Add_face name="wall_inner">  <Face_file_path> mesh/solid/mesh-surfaces/wall_interface.vtp </Face_file_path> </Add_face>
  <Add_face name="wall_outer">  <Face_file_path> mesh/solid/mesh-surfaces/wall_outer.vtp </Face_file_path> </Add_face>
  <Domain> 1 </Domain>
</Add_mesh>
<Add_projection name="wall_inner" > <Project_from_face> lumen_wall </Project_from_face> </Add_projection>
<Add_equation type="FSI" >
   <Coupled> true </Coupled> <Min_iterations> 1 </Min_iterations> <Max_iterations> 7 </Max_iterations> <Tolerance> 1e-6 </Tolerance>
   <Domain id="0" >
      <Equation> fluid </Equation> <Density> 1.06 </Density>
      <Viscosity model="Constant" > <Value> 0.04 </Value> </Viscosity>
      <Backflow_stabilization_coefficient> 0.2 </Backflow_stabilization_coefficient>
   </Domain>
   <Domain id="1" >
      <Equation> struct </Equation>
      <Constitutive_model type="neoHookean"> </Constitutive_model>
      <Dilational_penalty_model> M94 </Dilational_penalty_model>
      <Density> 1.0 </Density> <Elasticity_modulus> {E} </Elasticity_modulus> <Poisson_ratio> 0.49 </Poisson_ratio>
   </Domain>
   <Output type="Spatial" > <Velocity> true </Velocity> <Pressure> true </Pressure> <Displacement> true </Displacement> <WSS> true </WSS> </Output>
   <LS type="GMRES" >
      <Linear_algebra type="fsils" > <Preconditioner> fsils </Preconditioner> </Linear_algebra>
      <Tolerance> 1e-6 </Tolerance> <Max_iterations> 100 </Max_iterations> <Krylov_space_dimension> 50 </Krylov_space_dimension>
   </LS>
   <Add_BC name="lumen_inlet" >
      <Type> Dir </Type> <Time_dependence> Unsteady </Time_dependence>
      <Temporal_values_file_path> cap_inlet.flow </Temporal_values_file_path>
      <Profile> Parabolic </Profile> <Impose_flux> true </Impose_flux>
   </Add_BC>
   <Add_BC name="lumen_outlet" >
      <Type> Neu </Type> <Time_dependence> RCR </Time_dependence>
      <RCR_values> <Capacitance> {C} </Capacitance> <Distal_resistance> {Rd} </Distal_resistance>
        <Proximal_resistance> {Rp} </Proximal_resistance> <Distal_pressure> 0 </Distal_pressure> <Initial_pressure> 0 </Initial_pressure> </RCR_values>
   </Add_BC>
   <Add_BC name="wall_inlet" > <Type> Dir </Type> <Value> 0.0 </Value> <Impose_on_state_variable_integral> true </Impose_on_state_variable_integral>
      <Zero_out_perimeter> false </Zero_out_perimeter> <Effective_direction> (0, 0, 1) </Effective_direction> </Add_BC>
   <Add_BC name="wall_outlet" > <Type> Dir </Type> <Value> 0.0 </Value> <Impose_on_state_variable_integral> true </Impose_on_state_variable_integral>
      <Zero_out_perimeter> false </Zero_out_perimeter> <Effective_direction> (0, 0, 1) </Effective_direction> </Add_BC>
   <Add_BC name="wall_outer" > <Type> Robin </Type> <Stiffness> 1.0e4 </Stiffness> <Damping> 1.0e2 </Damping> <Apply_along_normal_direction> true </Apply_along_normal_direction> </Add_BC>
</Add_equation>
<Add_equation type="mesh" >
   <Coupled> true </Coupled> <Min_iterations> 1 </Min_iterations> <Max_iterations> 7 </Max_iterations> <Tolerance> 1e-6 </Tolerance>
   <Poisson_ratio> 0.3 </Poisson_ratio>
   <LS type="CG" > <Linear_algebra type="fsils" > <Preconditioner> fsils </Preconditioner> </Linear_algebra> <Tolerance> 1e-6 </Tolerance> </LS>
   <Add_BC name="lumen_inlet" > <Type> Dir </Type> <Value> 0.0 </Value> </Add_BC>
   <Add_BC name="lumen_outlet" > <Type> Dir </Type> <Value> 0.0 </Value> </Add_BC>
</Add_equation>
</svMultiPhysicsFile>
"""


def write_flow(path, inflow):
    with open(path, "w") as f:
        f.write(f"{len(inflow.time)} 16\n")
        for t_, q in zip(inflow.time, inflow.flow):
            f.write(f"{t_!r} {-q!r}\n")


def build_all(out_root, resolutions, E=1.0e7):
    inflow = bclib.aortic_waveform(30.0, 1.0)
    rcr, frac = bclib.rcr_from_targets({"cap_outlet": math.pi}, 30.0, 90.0, rp_fraction=0.09, total_compliance=1.5e-4)
    rcr = rcr["cap_outlet"]
    summary = {}
    for res, h in resolutions.items():
        fluid, solid, info = build_pipe(h)
        nhex = sum(1 for i in range(fluid.GetNumberOfCells()) if fluid.GetCellType(i) == vtk.VTK_HEXAHEDRON)
        print(res, info, "fluid hex", nhex, flush=True)
        # ---- CFD-only case (hex): same BCs as the tet pipe, via svmpsetup
        cfd_dir = os.path.join(out_root, "pipe_hex", res); os.makedirs(cfd_dir, exist_ok=True)
        names = {1: "wall", 2: "cap_inlet", 3: "cap_outlet"}
        values = {"mesh.scale_factor": 1.0, "general.time_step_size": 1e-3, "general.number_of_time_steps": 1000,
                  "general.increment_in_saving_vtk_files": 50, "general.increment_in_saving_restart_files": 50}
        write_case(cfd_dir, fluid, names, {"cap_inlet": inflow, "cap_outlet": rcr}, values, face_id_array_name="ModelFaceID")
        os.makedirs(os.path.join(cfd_dir, "mfem"), exist_ok=True)
        to_mfem.write_mfem_gmsh(fluid, os.path.join(cfd_dir, "mfem", f"pipe_hex_{res}.msh"), names_by_face_id=names)
        rep = quality.quality_report(fluid)
        json.dump(rep, open(os.path.join(cfd_dir, "quality.json"), "w"), indent=2)
        open(os.path.join(cfd_dir, "quality.md"), "w").write(quality.report_markdown(rep, f"pipe_hex / {res}") + "\n")
        json.dump({"name": "pipe_hex", "resolution": res, "grid": info, "elements": {"hex8": nhex, "points": fluid.GetNumberOfPoints()},
                   "faces": names, "boundary_conditions": {"cap_inlet": "aortic waveform, mean 30 cm3/s", "cap_outlet": rcr.__dict__},
                   "units": "cm"}, open(os.path.join(cfd_dir, "case.json"), "w"), indent=2)
        # ---- FSI case (hex fluid + hex wall)
        fsi_dir = os.path.join(out_root, "pipe_fsi_hex", res); os.makedirs(fsi_dir, exist_ok=True)
        write_mesh_complete(fluid, FaceTable([Face(1, "wall"), Face(2, "cap_inlet"), Face(3, "cap_outlet")]),
                            os.path.join(fsi_dir, "mesh", "fluid"), face_id_array_name="ModelFaceID")
        write_mesh_complete(solid, FaceTable([Face(1, "wall_outer"), Face(2, "cap_inlet_ring"), Face(3, "cap_outlet_ring"), Face(4, "wall_interface")]),
                            os.path.join(fsi_dir, "mesh", "solid"), face_id_array_name="ModelFaceID")
        write_flow(os.path.join(fsi_dir, "cap_inlet.flow"), inflow)
        open(os.path.join(fsi_dir, "solver.xml"), "w").write(FSI_XML.format(nsteps=1000, dt=1e-3, E=E, C=rcr.C, Rd=rcr.Rd, Rp=rcr.Rp))
        # combined MFEM mesh: merge fluid+solid (shared interface nodes), attributes 1/2, boundary ids offset for solid
        app = vtk.vtkAppendFilter(); app.MergePointsOn(); app.SetToleranceIsAbsolute(True); app.SetTolerance(1e-9)
        thr = vtk.vtkThreshold(); thr.SetInputData(solid)
        thr.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "ModelFaceID")
        thr.SetLowerThreshold(4); thr.SetUpperThreshold(4); thr.SetInvert(True); thr.Update()   # drop solid interface quads (duplicate of fluid wall)
        solid_off = vtk.vtkUnstructuredGrid(); solid_off.DeepCopy(thr.GetOutput())
        fa = solid_off.GetCellData().GetArray("ModelFaceID")
        for i in range(fa.GetNumberOfTuples()):
            v = int(fa.GetTuple1(i)); fa.SetTuple1(i, v + 10 if v > 0 else 0)
        app.AddInputData(fluid); app.AddInputData(solid_off); app.Update()
        both = app.GetOutput()
        os.makedirs(os.path.join(fsi_dir, "mfem"), exist_ok=True)
        to_mfem.write_mfem(both, os.path.join(fsi_dir, "mfem", f"pipe_fsi_hex_{res}.mesh"),
                           names_by_face_id={1: "fluid_wall_interface", 2: "fluid_cap_inlet", 3: "fluid_cap_outlet",
                                             11: "solid_wall_outer", 12: "solid_cap_inlet_ring", 13: "solid_cap_outlet_ring",
                                             "note": "attribute 1 (fluid wall) is the fluid-solid interface; solid side not duplicated"})
        rep_s = quality.quality_report(solid)
        json.dump({"fluid": rep, "solid": rep_s}, open(os.path.join(fsi_dir, "quality.json"), "w"), indent=2)
        open(os.path.join(fsi_dir, "quality.md"), "w").write(quality.report_markdown(rep, f"pipe_fsi_hex fluid / {res}") + "\n\n" +
                                                              quality.report_markdown(rep_s, f"pipe_fsi_hex solid / {res}") + "\n")
        nhs = sum(1 for i in range(solid.GetNumberOfCells()) if solid.GetCellType(i) == vtk.VTK_HEXAHEDRON)
        json.dump({"name": "pipe_fsi_hex", "resolution": res, "grid": info, "wall_thickness_cm": 0.1, "wall_layers": 3,
                   "elements": {"fluid_hex8": nhex, "solid_hex8": nhs}, "solid": {"model": "neoHookean", "E_dyn_cm2": E, "nu": 0.49, "density": 1.0},
                   "note": "fluid and solid share interface nodes (conforming); MFEM mesh has attributes 1=fluid, 2=solid; solid boundary attributes offset by 10",
                   "units": "cm"}, open(os.path.join(fsi_dir, "case.json"), "w"), indent=2)
        summary[res] = {"fluid_hex": nhex, "solid_hex": nhs, **info}
    return summary


if __name__ == "__main__":
    # cell sizes giving ~1.5M / 3M / 5M fluid hexes for R=1, L=11 (V = 34.6 cm^3; N ~ V/h^3)
    res = {"coarse": 0.0285, "medium": 0.0226, "fine": 0.0191}
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        res = {"test": 0.1}
    print(json.dumps(build_all(os.path.join(ROOT, "cases"), res), indent=1))
