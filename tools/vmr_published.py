#!/usr/bin/env python3
"""Add a VMR project's as-published mesh (Simulations/<job>/mesh-complete, the mesh its boundary
conditions were tuned on) as the `published` level of a case: repackaged through the same
svmpsetup/MFEM/quality path as the generated levels, with the project's own BCs.

    vmr_published.py specs/0054_H_PULMFON_TAT.json downloads/0054_H_PULMFON_TAT/0054_H_PULMFON_TAT

Credits: the published meshes are the VMR projects' own (Wilson et al. 2013, doi:10.1115/1.4025983);
repackaged with SlicerSimVascular svmpsetup. See REFERENCES.md.
"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import meshlib, quality, to_mfem, build_case
from svmpsetup.case import write_case

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def read(path):
    r = (vtk.vtkXMLPolyDataReader if path.endswith(".vtp") else vtk.vtkXMLUnstructuredGridReader)(); r.SetFileName(path); r.Update(); return r.GetOutput()


def combine(volume, exterior):
    """One unstructured grid holding the volume cells (face id 0) and the exterior cells
    (ModelFaceID) on the volume's points, joined through GlobalNodeID (float32/float64 point
    coordinates differ between the two files, so merging by position is not reliable)."""
    vg = vtk_to_numpy(volume.GetPointData().GetArray("GlobalNodeID")).astype(np.int64)
    eg = vtk_to_numpy(exterior.GetPointData().GetArray("GlobalNodeID")).astype(np.int64)
    lookup = -np.ones(vg.max() + 1, np.int64); lookup[vg] = np.arange(len(vg))
    emap = lookup[eg]
    if (emap < 0).any():
        raise ValueError("exterior GlobalNodeID not found in volume")
    tri = vtk.vtkTriangleFilter(); tri.SetInputData(exterior); tri.PassLinesOff(); tri.PassVertsOff(); tri.Update(); ext = tri.GetOutput()
    econn = vtk_to_numpy(ext.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    efid = vtk_to_numpy(ext.GetCellData().GetArray("ModelFaceID")).astype(np.int32)
    vtypes = np.array([volume.GetCellType(i) for i in range(volume.GetNumberOfCells())], np.uint8)
    voffs = vtk_to_numpy(volume.GetCells().GetOffsetsArray()); vconn = vtk_to_numpy(volume.GetCells().GetConnectivityArray())
    conn = np.concatenate([vconn, emap[econn].ravel()])
    offs = np.concatenate([voffs, voffs[-1] + np.arange(3, 3 * len(econn) + 1, 3)])
    types = np.concatenate([vtypes, np.full(len(econn), vtk.VTK_TRIANGLE, np.uint8)])
    ug = vtk.vtkUnstructuredGrid(); pts = vtk.vtkPoints(); pts.DeepCopy(volume.GetPoints()); ug.SetPoints(pts)
    ca = vtk.vtkCellArray()
    ca.SetData(numpy_to_vtk(offs.astype(np.int64), deep=True, array_type=vtk.VTK_ID_TYPE), numpy_to_vtk(conn.astype(np.int64), deep=True, array_type=vtk.VTK_ID_TYPE))
    ug.SetCells(numpy_to_vtk(types, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR), ca)
    fid = np.concatenate([np.zeros(len(vtypes), np.int32), efid])
    a = numpy_to_vtk(fid, deep=True); a.SetName("ModelFaceID"); ug.GetCellData().AddArray(a)
    return ug


def main(spec_path, project_dir, out_root=os.path.join(ROOT, "cases")):
    spec = json.load(open(spec_path))
    name = spec["name"]
    mc = glob.glob(os.path.join(project_dir, "Simulations", "*", "mesh-complete"))
    if not mc:
        raise FileNotFoundError("no Simulations/*/mesh-complete in project")
    mc = mc[0]
    volume = read(os.path.join(mc, "mesh-complete.mesh.vtu"))
    exterior = read(os.path.join(mc, "mesh-complete.exterior.vtp"))
    mesh = combine(volume, exterior)
    fids = np.unique(vtk_to_numpy(exterior.GetCellData().GetArray("ModelFaceID")))
    names = {int(k): v for k, v in spec["faces"].items()}
    missing = [int(f) for f in fids if int(f) not in names]
    if missing:
        raise ValueError(f"published mesh has face ids {missing} not in the model spec")
    # conditions from the project, areas from the published exterior
    summ = meshlib.face_summary(exterior, "ModelFaceID")
    areas = {names[k]: v["area"] for k, v in summ.items() if names[k].startswith("cap_")}
    conds, frac = build_case.make_conditions(spec, areas)
    period = max(c.period for c in conds.values() if hasattr(c, "period"))
    values = build_case.solver_values(spec, period)
    out_dir = os.path.join(out_root, name, "published"); os.makedirs(out_dir, exist_ok=True)
    case = write_case(out_dir, mesh, names, conds, values, face_id_array_name="ModelFaceID")
    os.makedirs(os.path.join(out_dir, "mfem"), exist_ok=True)
    attrs = to_mfem.write_mfem_gmsh(mesh, os.path.join(out_dir, "mfem", f"{name}_published.msh"), names_by_face_id=names)
    rep = quality.quality_report(mesh)
    json.dump(rep, open(os.path.join(out_dir, "quality.json"), "w"), indent=2)
    open(os.path.join(out_dir, "quality.md"), "w").write(quality.report_markdown(rep, f"{name} / published") + "\n")
    counts = meshlib.count_cells(mesh)
    meta = {"name": name, "resolution": "published", "elements": counts, "faces": names, "cap_areas_cm2": areas,
            "edge_length": None, "units": "cm", "mesher": "SimVascular (as published in the VMR project, TetGen/MeshSim)",
            "surface_source": spec.get("source"), "license": spec.get("license"),
            "boundary_conditions": build_case.conditions_json(conds, frac), "solver_values": values, "mfem_attributes": attrs,
            "note": f"as-published mesh from {os.path.relpath(mc, project_dir)}; the mesh the project's RCR/inflow values were tuned on"}
    json.dump(meta, open(os.path.join(out_dir, "case.json"), "w"), indent=2)
    print(name, "published:", counts)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
