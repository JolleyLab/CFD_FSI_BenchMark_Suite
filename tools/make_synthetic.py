#!/usr/bin/env python3
"""Idealized (synthetic) vascular surfaces, built from centerline tubes with a smooth radius
law, as a signed-distance field -> marching cubes -> VMTK smoothing/remeshing -> planar end
clips -> capped and face-labelled the way a SimVascular model is (ModelFaceID, CapID), so they
run through exactly the same build as the patient-specific cases.

Geometries (all in cm):
  aaa_idealized   infrarenal aorta (d 2.0) with a fusiform aneurysm (max d 5.5, length 9) and an
                  aorto-iliac bifurcation (iliacs d 1.2). 1 inlet, 2 outlets.
  tcpc_idealized  extracardiac total cavopulmonary connection: SVC (d 1.4) from above, IVC conduit
                  (d 1.8) from below offset 0.8 cm, into LPA/RPA (d 1.2/1.3). 2 inlets, 2 outlets.
  pipe            straight tube d 2.0 cm, length 12 cm (verification baseline).

Credits: VMTK surface remeshing and capping (Antiga et al. 2008, doi:10.1007/s11517-008-0420-1). See REFERENCES.md.
"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import meshlib
from meshlib import MeshingPipeline  # noqa (aliases vmtk modules)
from vmtk import vtkvmtkComputationalGeometryPython as vcg  # noqa: F401
from vmtk import vtkvmtkMiscPython as vmisc  # noqa: F401
from vmtk import vtkvmtkDifferentialGeometryPython as vdg


# ----------------------------------------------------------------------------- distance field
def _seg_dist(p, a, b):
    ab = b - a
    t = np.clip(((p - a) @ ab) / (ab @ ab), 0.0, 1.0)
    proj = a + t[:, None] * ab
    return np.linalg.norm(p - proj, axis=1), t


def tube_sdf(points, polyline, radii):
    """Signed distance to a tube around a polyline with per-vertex radius (linear along segments)."""
    best = np.full(len(points), np.inf)
    for i in range(len(polyline) - 1):
        d, t = _seg_dist(points, polyline[i], polyline[i + 1])
        r = radii[i] * (1 - t) + radii[i + 1] * t
        best = np.minimum(best, d - r)
    return best


def smooth_union(a, b, k=0.4):
    """Polynomial smooth-min, blends branches so junctions are rounded (k = blend radius, cm)."""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0.0, 1.0)
    return b * (1 - h) + a * h - k * h * (1 - h)


def field_to_surface(sdf_fn, bounds, spacing=0.05, smooth_iterations=30):
    xs = np.arange(bounds[0], bounds[1] + spacing, spacing)
    ys = np.arange(bounds[2], bounds[3] + spacing, spacing)
    zs = np.arange(bounds[4], bounds[5] + spacing, spacing)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    F = sdf_fn(P).reshape(X.shape)
    img = vtk.vtkImageData(); img.SetDimensions(*X.shape); img.SetOrigin(xs[0], ys[0], zs[0]); img.SetSpacing(spacing, spacing, spacing)
    arr = numpy_to_vtk(np.ascontiguousarray(F.transpose(2, 1, 0)).ravel(), deep=True); arr.SetName("sdf")
    img.GetPointData().SetScalars(arr)
    mc = vtk.vtkFlyingEdges3D(); mc.SetInputData(img); mc.SetValue(0, 0.0); mc.ComputeNormalsOff(); mc.Update()
    surf = mc.GetOutput()
    conn = vtk.vtkPolyDataConnectivityFilter(); conn.SetInputData(surf); conn.SetExtractionModeToLargestRegion(); conn.Update()
    sm = vtk.vtkWindowedSincPolyDataFilter(); sm.SetInputData(conn.GetOutput()); sm.SetNumberOfIterations(smooth_iterations)
    sm.SetPassBand(0.05); sm.BoundarySmoothingOff(); sm.NonManifoldSmoothingOn(); sm.NormalizeCoordinatesOn(); sm.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputData(sm.GetOutput()); cl.Update()
    return cl.GetOutput()


def clip_with_plane(surface, origin, normal):
    """Remove the part of the surface beyond the plane, `normal` pointing outward (toward the
    discarded end)."""
    plane = vtk.vtkPlane(); plane.SetOrigin(*origin); plane.SetNormal(*[-float(c) for c in normal])
    c = vtk.vtkClipPolyData(); c.SetInputData(surface); c.SetClipFunction(plane); c.InsideOutOff(); c.Update()
    conn = vtk.vtkPolyDataConnectivityFilter(); conn.SetInputData(c.GetOutput()); conn.SetExtractionModeToLargestRegion(); conn.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputData(conn.GetOutput()); cl.Update()
    return cl.GetOutput()


def remesh(surface, edge_length, reduction=0.7, iterations=5):
    """Quadric-decimate the marching-cubes triangles (which are sliver-ridden) and remesh them
    with VMTK's surface remesher to near-equilateral triangles of the given edge length."""
    dec = vtk.vtkQuadricDecimation(); dec.SetInputData(surface); dec.SetTargetReduction(reduction); dec.Update()
    rm = vdg.vtkvmtkPolyDataSurfaceRemeshing(); rm.SetInputData(dec.GetOutput())
    rm.SetElementSizeModeToTargetArea(); rm.SetTargetArea(edge_length ** 2 * math.sqrt(3) / 4)
    rm.SetNumberOfIterations(iterations); rm.SetCollapseAngleThreshold(0.2); rm.SetTriangleSplitFactor(5.0)
    rm.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputData(rm.GetOutput()); cl.Update()
    return cl.GetOutput()


def cap_and_label(open_surface, inlet_centroids, tol=0.6):
    """Cap an open surface with the pipeline's capper, then write SimVascular-style arrays:
    ModelFaceID (1 = wall, caps numbered upward) and CapID (-1 wall, 1 inlet, 2 outlet)."""
    p = MeshingPipeline()
    s = vtk.vtkPolyData(); s.DeepCopy(open_surface)
    ids = vtk.vtkIntArray(); ids.SetName("CellEntityIds"); ids.SetNumberOfTuples(s.GetNumberOfCells()); ids.Fill(1)
    s.GetCellData().AddArray(ids)
    capped = p.capSurface(s, "CellEntityIds", "simple")
    if isinstance(capped, tuple):
        capped = capped[0]
    tri = vtk.vtkTriangleFilter(); tri.SetInputData(capped); tri.PassLinesOff(); tri.PassVertsOff(); tri.Update()
    capped = tri.GetOutput()
    ents = vtk_to_numpy(capped.GetCellData().GetArray("CellEntityIds")).copy()
    capid = np.full(len(ents), -1, np.int32)
    summ = meshlib.face_summary(capped, "CellEntityIds")
    inlet_centroids = [np.array(c) for c in inlet_centroids]
    cap_info = {}
    for fid, info in summ.items():
        if fid == 1:
            continue
        c = np.array(info["centroid"])
        is_inlet = any(np.linalg.norm(c - ic) < tol for ic in inlet_centroids)
        capid[ents == fid] = 1 if is_inlet else 2
        cap_info[fid] = dict(info, inlet=bool(is_inlet))
    mf = numpy_to_vtk(ents.astype(np.int32), deep=True); mf.SetName("ModelFaceID")
    ci = numpy_to_vtk(capid, deep=True); ci.SetName("CapID")
    capped.GetCellData().RemoveArray("CellEntityIds")
    capped.GetCellData().AddArray(mf); capped.GetCellData().AddArray(ci)
    return capped, cap_info


# ----------------------------------------------------------------------------- geometries
def aaa_idealized():
    """Infrarenal AAA with bifurcation. z axis = cranio-caudal (inlet at top, z = 0)."""
    # aorta centerline, slight anterior bow (y) through the sac
    z = np.linspace(0, -14, 57)
    y = 0.6 * np.exp(-((z + 7.5) / 3.0) ** 2)          # sac bulges anteriorly
    x = np.zeros_like(z)
    aorta = np.stack([x, y, z], axis=1)
    r = 1.0 + 1.75 * np.exp(-((z + 7.5) / 2.6) ** 2)    # d 2.0 cm neck -> 5.5 cm max, back to 2.0
    r = np.where(z < -12.5, 1.0 - 0.15 * (-12.5 - z), r)  # taper toward bifurcation (d 1.7 at -14)
    def sdf(P):
        f = tube_sdf(P, aorta, r)
        for sx in (-1, 1):
            il = np.array([[0, 0, -13.6], [sx * 0.9, 0.1, -15.5], [sx * 2.1, 0.3, -18.5], [sx * 2.8, 0.4, -21.0]])
            ri = np.array([0.8, 0.7, 0.62, 0.6])
            f = smooth_union(f, tube_sdf(P, il, ri), k=0.5)
        return f
    surf = field_to_surface(sdf, (-5.5, 5.5, -3.5, 4.5, -22.0, 1.0), spacing=0.06)
    surf = clip_with_plane(surf, (0, 0, -0.6), (0, 0, 1))                  # inlet plane (discard above)
    surf = clip_with_plane(surf, (-2.6, 0.37, -20.3), (-0.7, 0.1, -1.85))    # left iliac outlet
    surf = clip_with_plane(surf, (2.6, 0.37, -20.3), (0.7, 0.1, -1.85))      # right iliac outlet
    return surf, [(0, 0, -0.6)]


def tcpc_idealized():
    """Extracardiac TCPC. Pulmonary arteries run along x (patient left-right), cavae along z."""
    pa = np.array([[-5.0, 0.0, 0.0], [-2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    rpa = np.array([0.60, 0.65, 0.75, 0.70, 0.65])   # LPA d 1.2 / RPA d 1.3 with a wider confluence
    svc = np.array([[0.0, 0.0, 6.0], [0.0, 0.0, 2.0], [0.0, 0.0, 0.0]])
    rsvc = np.array([0.70, 0.72, 0.85])            # d 1.4, flared at the anastomosis
    ivc = np.array([[0.8, 0.3, -7.0], [0.8, 0.3, -3.0], [0.8, 0.2, -1.2], [0.8, 0.1, 0.0]])
    rivc = np.array([0.90, 0.90, 0.95, 1.05])      # conduit d 1.8, offset 0.8 cm toward RPA, flared
    def sdf(P):
        f = tube_sdf(P, pa, rpa)
        f = smooth_union(f, tube_sdf(P, svc, rsvc), k=0.6)
        f = smooth_union(f, tube_sdf(P, ivc, rivc), k=0.7)
        return f
    surf = field_to_surface(sdf, (-6.0, 6.0, -2.0, 2.5, -8.0, 7.0), spacing=0.05)
    surf = clip_with_plane(surf, (0, 0, 5.4), (0, 0, 1))        # SVC inlet
    surf = clip_with_plane(surf, (0.8, 0.3, -6.4), (0, 0, -1))  # IVC inlet
    surf = clip_with_plane(surf, (-4.4, 0, 0), (-1, 0, 0))      # LPA outlet
    surf = clip_with_plane(surf, (4.4, 0, 0), (1, 0, 0))        # RPA outlet
    return surf, [(0, 0, 5.4), (0.8, 0.3, -6.4)]


def pipe_analytic(R=1.0, z0=0.5, z1=11.5, h=0.08):
    """Closed, labelled cylinder built directly (structured triangles, fan caps): ModelFaceID
    1 = wall, 2 = outlet cap (z1), 3 = inlet cap (z0); CapID -1 / 2 / 1."""
    n_ang = int(round(2 * math.pi * R / h)); n_z = int(round((z1 - z0) / (h * math.sqrt(3) / 2))) + 1
    ang = np.linspace(0, 2 * math.pi, n_ang, endpoint=False); z = np.linspace(z0, z1, n_z)
    pts = []; tris = []; fid = []
    for k, zz in enumerate(z):
        off = (0.5 * (k % 2)) * (ang[1] - ang[0])
        for a in ang + off:
            pts.append((R * math.cos(a), R * math.sin(a), zz))
    def idx(k, i): return k * n_ang + (i % n_ang)
    for k in range(n_z - 1):
        for i in range(n_ang):
            if k % 2 == 0:
                tris += [(idx(k, i), idx(k, i + 1), idx(k + 1, i)), (idx(k, i + 1), idx(k + 1, i + 1), idx(k + 1, i))]
            else:
                tris += [(idx(k, i), idx(k, i + 1), idx(k + 1, i + 1)), (idx(k, i), idx(k + 1, i + 1), idx(k + 1, i))]
            fid += [1, 1]
    # caps: rim points + concentric interior rings, triangulated with Delaunay2D (convex disc)
    for k, f, orient in ((0, 3, -1), (n_z - 1, 2, 1)):
        rim = [idx(k, i) for i in range(n_ang)]
        disc_ids = list(rim)
        n_r = max(1, int(round(R / h)))
        for j in range(1, n_r):
            rr = R * (1 - j / n_r); m = max(6, int(round(2 * math.pi * rr / h)))
            off = (0.5 * (j % 2)) * (2 * math.pi / m)
            for a in np.linspace(0, 2 * math.pi, m, endpoint=False) + off:
                disc_ids.append(len(pts)); pts.append((rr * math.cos(a), rr * math.sin(a), z[k]))
        disc_ids.append(len(pts)); pts.append((0.0, 0.0, z[k]))
        dp = vtk.vtkPolyData(); vp = vtk.vtkPoints()
        for i in disc_ids:
            vp.InsertNextPoint(pts[i][0], pts[i][1], 0.0)
        dp.SetPoints(vp)
        dl = vtk.vtkDelaunay2D(); dl.SetInputData(dp); dl.SetTolerance(1e-6); dl.Update()
        dc = vtk_to_numpy(dl.GetOutput().GetPolys().GetConnectivityArray()).reshape(-1, 3)
        for t in dc:
            a, b, c = (disc_ids[t[0]], disc_ids[t[1]], disc_ids[t[2]])
            # orient: Delaunay2D gives CCW in xy (normal +z); inlet cap normal must be -z
            tri = (a, b, c) if orient > 0 else (a, c, b)
            tris.append(tri); fid.append(f)
    pd = vtk.vtkPolyData(); vp = vtk.vtkPoints(); vp.SetData(numpy_to_vtk(np.array(pts), deep=True)); pd.SetPoints(vp)
    ca = vtk.vtkCellArray()
    T = np.array(tris, np.int64)
    ca.SetData(numpy_to_vtk(np.arange(0, 3 * len(T) + 1, 3).astype(np.int64), deep=True, array_type=vtk.VTK_ID_TYPE), numpy_to_vtk(T.ravel(), deep=True, array_type=vtk.VTK_ID_TYPE))
    pd.SetPolys(ca)
    mf = numpy_to_vtk(np.array(fid, np.int32), deep=True); mf.SetName("ModelFaceID"); pd.GetCellData().AddArray(mf)
    ci = numpy_to_vtk(np.array([-1 if f == 1 else (1 if f == 3 else 2) for f in fid], np.int32), deep=True); ci.SetName("CapID"); pd.GetCellData().AddArray(ci)
    cl = vtk.vtkCleanPolyData(); cl.SetInputData(pd); cl.Update()
    return cl.GetOutput()


def pipe():
    line = np.array([[0, 0, 0.0], [0, 0, 12.0]]); rr = np.array([1.0, 1.0])
    surf = field_to_surface(lambda P: tube_sdf(P, line, rr), (-1.5, 1.5, -1.5, 1.5, -0.5, 12.5), spacing=0.05)
    surf = clip_with_plane(surf, (0, 0, 0.5), (0, 0, -1))
    surf = clip_with_plane(surf, (0, 0, 11.5), (0, 0, 1))
    return surf, [(0, 0, 0.5)]


GEOMETRIES = {"aaa_idealized": aaa_idealized, "tcpc_idealized": tcpc_idealized, "pipe": pipe}


def remesh_capped(capped, edge_length):
    """Remesh the capped, labelled surface (caps included) keeping ModelFaceID, re-deriving CapID."""
    fid = vtk_to_numpy(capped.GetCellData().GetArray("ModelFaceID")); cid = vtk_to_numpy(capped.GetCellData().GetArray("CapID"))
    capid_of = {int(f): int(cid[fid == f][0]) for f in np.unique(fid)}
    s = vtk.vtkPolyData(); s.DeepCopy(capped)
    ids = vtk.vtkIntArray(); ids.SetName("CellEntityIds"); ids.SetNumberOfTuples(s.GetNumberOfCells())
    for i, f in enumerate(fid):
        ids.SetValue(i, int(f))
    s.GetCellData().AddArray(ids)
    rm = vdg.vtkvmtkPolyDataSurfaceRemeshing(); rm.SetInputData(s); rm.SetCellEntityIdsArrayName("CellEntityIds")
    rm.SetElementSizeModeToTargetArea(); rm.SetTargetArea(edge_length ** 2 * math.sqrt(3) / 4)
    rm.SetNumberOfIterations(8); rm.SetCollapseAngleThreshold(0.2); rm.SetTriangleSplitFactor(5.0); rm.SetPreserveBoundaryEdges(0)
    rm.Update()
    out = rm.GetOutput()
    ents = vtk_to_numpy(out.GetCellData().GetArray("CellEntityIds")).astype(np.int32)
    mf = numpy_to_vtk(ents, deep=True); mf.SetName("ModelFaceID")
    ci = numpy_to_vtk(np.array([capid_of[int(e)] for e in ents], np.int32), deep=True); ci.SetName("CapID")
    out.GetCellData().RemoveArray("CellEntityIds"); out.GetCellData().AddArray(mf); out.GetCellData().AddArray(ci)
    return out


def build(name, out_dir, remesh_h=0.08):
    if name == "pipe":
        capped = remesh_capped(pipe_analytic(h=remesh_h), remesh_h)
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, f"{name}_model.vtp"); meshlib.write_dataset(capped, path)
        summ = meshlib.face_summary(capped, "ModelFaceID")
        caps = {k: dict(v, inlet=(k == 3)) for k, v in summ.items() if k != 1}
        info = {"surface": path, "n_triangles": capped.GetNumberOfCells(), "caps": {int(k): v for k, v in caps.items()}}
        json.dump(info, open(os.path.join(out_dir, f"{name}_faces.json"), "w"), indent=2)
        return info
    surf, inlets = GEOMETRIES[name]()
    surf = remesh(surf, remesh_h)
    capped, caps = cap_and_label(surf, inlets)
    capped = remesh_capped(capped, remesh_h)
    caps = {k: dict(meshlib.face_summary(capped, "ModelFaceID")[k], inlet=v["inlet"]) for k, v in caps.items()}
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}_model.vtp")
    meshlib.write_dataset(capped, path)
    info = {"surface": path, "n_triangles": capped.GetNumberOfCells(), "caps": {int(k): v for k, v in caps.items()}}
    json.dump(info, open(os.path.join(out_dir, f"{name}_faces.json"), "w"), indent=2)
    return info


if __name__ == "__main__":
    out = sys.argv[2] if len(sys.argv) > 2 else "surfaces"
    for n in (sys.argv[1].split(",") if len(sys.argv) > 1 else GEOMETRIES):
        info = build(n, out)
        print(n, json.dumps({k: (round(v["eq_diameter"], 2), [round(c, 2) for c in v["centroid"]], v["inlet"]) for k, v in info["caps"].items()}))
