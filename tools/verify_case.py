#!/usr/bin/env python3
"""Independent checks on a built case directory (does not use the code that built it):

  * mesh-complete: volume VTU holds one element type, no boundary cells; GlobalNodeID is 1..N;
    every mesh-surfaces/*.vtp's GlobalNodeID points at coincident volume nodes; exterior is
    closed (no open edges) and its face areas equal the sum of the per-face files;
  * solver.xml references files that exist, and every Add_face has a BC or is a wall;
  * MFEM input (mfem/*.msh, Gmsh 2.2; a legacy mfem/*.mesh is checked too if present): header,
    counts, physical names and boundary attribute set equal *_attributes.json, all elements
    positively oriented, every boundary face is a face of exactly one element.
  * quality gate: no inverted elements; isotropic (generated) meshes FAIL below a 2 deg minimum or
    above a 178 deg maximum dihedral angle and WARN above 0.5 % of elements with scaled Jacobian
    < 0.2; boundary-layer meshes (VMR published, flagged in case.json) are only checked for
    inversion, since their thin wall elements are intentional;
  * volume check: sum of element volumes equals the volume enclosed by the exterior surface
    (catches holes and overlaps) to 0.1 %;
  * MFEM load test (if PyMFEM is installed): the mesh loads in MFEM, MFEM's own orientation checks
    need no fixes, and the volume MFEM integrates matches the VTK volume.
Prints PASS/WARN/FAIL lines; exit code 1 on any failure.

Independent of the build code by design: reads only the written files.
"""
import glob, json, os, re, sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

fails = 0


def warn(cond, msg):
    print(("PASS " if cond else "WARN ") + msg)


def check(cond, msg):
    global fails
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails += 1


def read(path):
    r = (vtk.vtkXMLPolyDataReader if path.endswith(".vtp") else vtk.vtkXMLUnstructuredGridReader)(); r.SetFileName(path); r.Update(); return r.GetOutput()


def area(pd):
    """Surface area of a polydata of triangles and/or quads."""
    tri = vtk.vtkTriangleFilter(); tri.SetInputData(pd); tri.Update()
    mp = vtk.vtkMassProperties(); mp.SetInputData(tri.GetOutput()); mp.Update(); return mp.GetSurfaceArea()


def check_mesh_complete(mesh_dir):
    vol = read(os.path.join(mesh_dir, "mesh-complete.mesh.vtu"))
    types = np.unique([vol.GetCellType(i) for i in range(vol.GetNumberOfCells())])
    check(len(types) == 1 and types[0] in (vtk.VTK_TETRA, vtk.VTK_HEXAHEDRON), f"{mesh_dir}: single volume element type {types}")
    gid = vtk_to_numpy(vol.GetPointData().GetArray("GlobalNodeID"))
    check(gid.min() == 1 and gid.max() == vol.GetNumberOfPoints() and len(np.unique(gid)) == len(gid), f"{mesh_dir}: GlobalNodeID is 1..N ({vol.GetNumberOfPoints():,} nodes)")
    P = vtk_to_numpy(vol.GetPoints().GetData())
    ext_path = os.path.join(mesh_dir, "mesh-complete.exterior.vtp")
    if not os.path.exists(ext_path):
        print(f"NOTE {mesh_dir}: no mesh-complete.exterior.vtp (allowed; solver reads mesh-surfaces only)")
        for f in sorted(glob.glob(os.path.join(mesh_dir, "mesh-surfaces", "*.vtp"))):
            face = read(f); fg = vtk_to_numpy(face.GetPointData().GetArray("GlobalNodeID"))
            d = np.linalg.norm(vtk_to_numpy(face.GetPoints().GetData()) - P[fg - 1], axis=1).max()
            check(d < 1e-6, f"  {os.path.basename(f)}: GlobalNodeID maps onto volume nodes (max offset {d:.1e})")
        return vol
    ext = read(ext_path)
    fe = vtk.vtkFeatureEdges(); fe.SetInputData(ext); fe.BoundaryEdgesOn(); fe.FeatureEdgesOff(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOn(); fe.Update()
    check(fe.GetOutput().GetNumberOfCells() == 0, f"{mesh_dir}: exterior closed and manifold ({ext.GetNumberOfCells():,} boundary cells)")
    ext_area = area(ext)
    total = 0.0
    for f in sorted(glob.glob(os.path.join(mesh_dir, "mesh-surfaces", "*.vtp"))):
        face = read(f)
        fg = vtk_to_numpy(face.GetPointData().GetArray("GlobalNodeID"))
        d = np.linalg.norm(vtk_to_numpy(face.GetPoints().GetData()) - P[fg - 1], axis=1).max()
        check(d < 1e-9, f"  {os.path.basename(f)}: GlobalNodeID maps onto volume nodes (max offset {d:.1e}), {face.GetNumberOfCells():,} cells")
        if os.path.basename(f).startswith("cap_"):
            # resolution across the cap: equivalent diameter over mean edge length
            ed = vtk.vtkExtractEdges(); ed.SetInputData(face); ed.Update()
            e = ed.GetOutput(); EP = vtk_to_numpy(e.GetPoints().GetData()); EC = vtk_to_numpy(e.GetLines().GetConnectivityArray()).reshape(-1, 2)
            h = float(np.linalg.norm(EP[EC[:, 0]] - EP[EC[:, 1]], axis=1).mean())
            across = 2.0 * np.sqrt(area(face) / np.pi) / h
            is_inflow = "inflow" in os.path.basename(f) or "inlet" in os.path.basename(f)
            warn(across >= (10 if is_inflow else 5), f"  {os.path.basename(f)}: {across:.1f} elements across the cap (warn below {10 if is_inflow else 5})")
        total += area(face)
    check(abs(total - ext_area) / ext_area < 1e-6, f"{mesh_dir}: face areas sum to exterior area ({total:.4f} vs {ext_area:.4f} cm^2)")
    return vol


def check_solver_xml(case_dir):
    xml = open(os.path.join(case_dir, "solver.xml")).read()
    for ref in re.findall(r"<(?:Mesh_file_path|Face_file_path|Temporal_values_file_path)>\s*([^<\s]+)\s*<", xml):
        check(os.path.exists(os.path.join(case_dir, ref)), f"solver.xml references {ref}")
    faces = re.findall(r"<Add_face name=\"([^\"]+)\"", xml)
    bcs = re.findall(r"<Add_BC name=\"([^\"]+)\"", xml)
    for f in faces:
        check(f in bcs or f.startswith("wall") or f.endswith("wall"), f"face {f} has a BC (or is a wall)")


def read_mfem_ascii(path):
    """MFEM mesh v1.0 -> (el, bd, V); rows are [attribute, geometry, vertex ids...] (0-based)."""
    with open(path) as f:
        head = f.readline().strip()
        check(head == "MFEM mesh v1.0", f"{os.path.basename(path)}: header")
        text = f.read()
    def section(name):
        m = re.search(rf"\n{name}\n(\d+)\n", text); return int(m.group(1)), m.end()
    ne, pe = section("elements"); nb, pb = section("boundary"); nv, pv = section("vertices")
    el = np.loadtxt(io_lines(text, pe, ne), dtype=np.int64)
    bd = np.loadtxt(io_lines(text, pb, nb), dtype=np.int64)
    verts_dim = int(text[pv:].split("\n", 1)[0]); check(verts_dim == 3, "vertices dimension 3")
    V = np.loadtxt(io_lines(text, text.find("\n", pv) + 1, nv), dtype=float)
    check(len(el) == ne and len(bd) == nb and len(V) == nv, f"counts: {ne:,} elements, {nb:,} boundary, {nv:,} vertices")
    return el, bd, V, None


GMSH_TO_MFEM = {4: (4, 4, 3), 5: (5, 8, 3), 2: (2, 3, 2), 3: (3, 4, 2), 6: (6, 6, 3)}   # gmsh type -> (MFEM geom, nodes, dim)


def read_gmsh22(path):
    """Gmsh 2.2 (ASCII or binary, written by any tool) -> (el, bd, V, physical names), rows as above.
    The physical tag (first tag) is the MFEM attribute, as MFEM's own Gmsh reader takes it."""
    data = open(path, "rb").read()
    def block(name):
        a = data.index(f"${name}\n".encode()) + len(name) + 2
        return a, data.index(f"$End{name}".encode(), a)
    a, b = block("MeshFormat")
    ver, ftype, dsize = data[a:b].split(b"\n")[0].split()
    check(ver == b"2.2" and dsize == b"8", f"{os.path.basename(path)}: Gmsh 2.2 header ({'binary' if ftype == b'1' else 'ASCII'})")
    binary = ftype == b"1"
    names = {}
    if b"$PhysicalNames" in data:
        a, b = block("PhysicalNames")
        for ln in data[a:b].decode().strip().split("\n")[1:]:
            d, t, nm = ln.split(None, 2); names[(int(d), int(t))] = nm.strip('"')
    a, b = block("Nodes")
    nl = data.index(b"\n", a); nv = int(data[a:nl]); p = nl + 1
    if binary:
        rec = np.frombuffer(data, dtype=[("id", "<i4"), ("x", "<f8", (3,))], count=nv, offset=p)
        ids, V = rec["id"].astype(np.int64), rec["x"].astype(float)
    else:
        arr = np.loadtxt(io_lines(data[p:b].decode(), 0, nv)); ids, V = arr[:, 0].astype(np.int64), arr[:, 1:4]
    check(np.array_equal(ids, np.arange(1, nv + 1)), "Gmsh node ids are 1..N")
    a, b = block("Elements")
    nl = data.index(b"\n", a); ne_tot = int(data[a:nl]); p = nl + 1
    rows = {3: [], 2: []}
    if binary:
        seen = 0
        while seen < ne_tot:
            gt, n, nt = np.frombuffer(data, "<i4", 3, p); p += 12
            g, nn, dim = GMSH_TO_MFEM[int(gt)]
            blk = np.frombuffer(data, "<i4", n * (1 + nt + nn), p).reshape(n, 1 + nt + nn); p += blk.nbytes
            r = np.empty((n, 2 + nn), np.int64); r[:, 0] = blk[:, 1]; r[:, 1] = g; r[:, 2:] = blk[:, 1 + nt:] - 1
            rows[dim].append(r); seen += n
    else:
        for ln in data[p:b].decode().strip().split("\n"):
            v = list(map(int, ln.split())); g, nn, dim = GMSH_TO_MFEM[v[1]]; nt = v[2]
            rows[dim].append(np.array([[v[3], g] + [x - 1 for x in v[3 + nt:3 + nt + nn]]], np.int64))
    el = np.concatenate(rows[3]); bd = np.concatenate(rows[2])
    check(len(el) + len(bd) == ne_tot, f"counts: {len(el):,} elements, {len(bd):,} boundary, {nv:,} vertices")
    return el, bd, V, names


def _face_keys(F, nvtx):
    """Sorted vertex ids of each face -> structured (k1, k2) keys; overflow-safe for any mesh size
    (k1 = v0*N + v1, k2 = v2*N + v3, each < N^2)."""
    F = np.sort(F.astype(np.int64), axis=1)
    k = np.empty(len(F), dtype=[("k1", "<i8"), ("k2", "<i8")])
    k["k1"] = F[:, 0] * nvtx + F[:, 1]
    k["k2"] = F[:, 2] * nvtx + (F[:, 3] if F.shape[1] > 3 else 0)
    return k


def face_owners(n, faces, bnd, nvtx, el_name, face_name):
    """Every boundary face is a face of 1 element (exterior) or 2 (interface); the exterior count
    equals the number of element faces owned once."""
    uk, cnt = np.unique(np.concatenate([_face_keys(n[:, list(f)], nvtx) for f in faces]), return_counts=True)
    bk = _face_keys(bnd, nvtx)
    pos = np.minimum(np.searchsorted(uk, bk), len(uk) - 1)
    owners = np.where(uk[pos] == bk, cnt[pos], 0)
    interior = int((owners == 2).sum())
    check(bool(np.isin(owners, (1, 2)).all()), f"every boundary {face_name} belongs to 1 (exterior) or 2 (interface) {el_name}; {interior:,} interface faces")
    single = int((cnt == 1).sum())
    check(single == len(bnd) - interior, f"exterior boundary count equals single-owner {el_name[:-1] if el_name == "tets" else "hex"} faces ({single:,} vs {len(bnd) - interior:,})")


def check_mfem(path, attr_json):
    """Structure of the mesh MFEM will read (.msh Gmsh 2.2, or legacy .mesh MFEM v1.0)."""
    attrs = json.load(open(attr_json))
    el, bd, V, names = read_gmsh22(path) if path.endswith(".msh") else read_mfem_ascii(path)
    ne, nb = len(el), len(bd)
    if names is not None:
        want = {(3, int(k)): v for k, v in attrs["element_attributes"].items()}
        want.update({(2, int(k)): v for k, v in attrs["boundary_attributes"].items()})
        check(names == want, f"Gmsh physical names match attributes.json ({len(want)} groups)")
    geoms = np.unique(el[:, 1]); check(all(g in (4, 5) for g in geoms), f"element geometries {geoms}")
    battrs = set(int(a) for a in np.unique(bd[:, 0])); check(battrs == set(int(k) for k in attrs["boundary_attributes"]), f"boundary attributes {sorted(battrs)} match attributes.json")
    if (el[:, 1] == 4).all():
        n = el[:, 2:6]; a_, b_, c_, d_ = V[n[:, 0]], V[n[:, 1]], V[n[:, 2]], V[n[:, 3]]
        vol = np.einsum("ij,ij->i", np.cross(b_ - a_, c_ - a_), d_ - a_) / 6.0
        check((vol > 0).all(), f"all {ne:,} tets positive volume (min {vol.min():.2e}, total {vol.sum():.4f})")
        del a_, b_, c_, d_, vol
        face_owners(n, ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)), bd[:, 2:5], len(V), "tets", "tri")
    if (el[:, 1] == 5).all():
        n = el[:, 2:10]
        v0 = V[n[:, 0]]; j = np.einsum("ij,ij->i", np.cross(V[n[:, 1]] - v0, V[n[:, 3]] - v0), V[n[:, 4]] - v0)
        check((j > 0).all(), f"all {ne:,} hexes positively oriented at node 0 (min corner volume {j.min():.2e})")
        del v0, j
        face_owners(n, ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)), bd[:, 2:6], len(V), "hexes", "quad")
    check(V.min() > -1e3 and V.max() < 1e3, "vertex coordinates in a sane cm range")


def io_lines(text, start, n):
    import io
    lines = text[start:].split("\n", n)
    return io.StringIO("\n".join(lines[:n]))


def check_quality_and_volume(case_dir, vol):
    qj = os.path.join(case_dir, "quality.json")
    cj = os.path.join(case_dir, "case.json")
    if not os.path.exists(qj):
        return None
    q = json.load(open(qj)); q = q.get("fluid", q)
    c = json.load(open(cj)) if os.path.exists(cj) else {}
    boundary_layer = c.get("resolution") == "published" or "published" in str(c.get("mesher", ""))
    for et, e in q["element_types"].items():
        check(e.get("negative_volume_count", 0) == 0, f"quality: no inverted {et} ({e['count']:,} elements)")
        m = e["metrics"]
        if boundary_layer:
            print(f"NOTE quality: boundary-layer mesh; anisotropy expected (max aspect ratio {m['aspect_ratio']['max']:.0f})")
            continue
        if et == "tet4":
            check(m["min_dihedral_angle_deg"]["min"] >= 2.0, f"quality: min dihedral {m['min_dihedral_angle_deg']['min']:.1f} deg >= 2")
            check(m["max_dihedral_angle_deg"]["max"] <= 178.0, f"quality: max dihedral {m['max_dihedral_angle_deg']['max']:.1f} deg <= 178")
        frac = e.get("fraction_scaled_jacobian_lt_0p2", 0.0)
        warn(frac <= 0.005, f"quality: {100*frac:.3f} % of {et} with scaled Jacobian < 0.2 (warn above 0.5 %)")
    # volume: elements vs enclosed by exterior
    mq = vtk.vtkMeshQuality(); mq.SetInputData(vol); mq.SaveCellQualityOn()
    mq.SetTetQualityMeasureToVolume(); mq.SetHexQualityMeasureToVolume(); mq.Update()
    v_el = float(vtk_to_numpy(mq.GetOutput().GetCellData().GetArray("Quality")).sum())
    ext = os.path.join(case_dir, "mesh", "mesh-complete.exterior.vtp")
    if os.path.exists(ext):
        tri = vtk.vtkTriangleFilter(); tri.SetInputData(read(ext)); tri.Update()
        mp = vtk.vtkMassProperties(); mp.SetInputData(tri.GetOutput()); mp.Update()
        v_surf = mp.GetVolume()
        check(abs(v_el - v_surf) / v_surf < 1e-3, f"volume: elements {v_el:.4f} vs enclosed {v_surf:.4f} cm^3 ({100*abs(v_el-v_surf)/v_surf:.4f} %)")
    return v_el


def check_mfem_load(path, v_ref):
    try:
        import mfem.ser as mfem
    except ImportError:
        print("NOTE MFEM load test skipped (PyMFEM not installed)")
        return
    m = mfem.Mesh(path, 1, 1)
    check(m.GetNE() > 0, f"MFEM loads {os.path.basename(path)}: {m.GetNE():,} elements, {m.GetNBE():,} boundary, attributes {m.attributes.ToList()}")
    check(m.CheckElementOrientation(False) == 0, "MFEM CheckElementOrientation: no elements needed fixing")
    check(m.CheckBdrElementOrientation(False) == 0, "MFEM CheckBdrElementOrientation: no boundary elements needed fixing")
    fes = mfem.FiniteElementSpace(m, mfem.H1_FECollection(1, m.Dimension()))
    lf = mfem.LinearForm(fes); lf.AddDomainIntegrator(mfem.DomainLFIntegrator(mfem.ConstantCoefficient(1.0))); lf.Assemble()
    v = lf.Sum()
    if v_ref:
        check(abs(v - v_ref) / v_ref < 1e-6, f"MFEM volume {v:.4f} matches VTK volume {v_ref:.4f} cm^3")


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    for case_dir in [a for a in sys.argv[1:] if not a.startswith("--")]:
        print(f"== {case_dir}")
        vol_total = None
        if os.path.isdir(os.path.join(case_dir, "mesh", "mesh-surfaces")):
            vmesh = check_mesh_complete(os.path.join(case_dir, "mesh"))
            vol_total = check_quality_and_volume(case_dir, vmesh)
        for sub in ("fluid", "solid"):
            if os.path.isdir(os.path.join(case_dir, "mesh", sub)):
                check_mesh_complete(os.path.join(case_dir, "mesh", sub))
        if os.path.exists(os.path.join(case_dir, "solver.xml")):
            check_solver_xml(case_dir)
        if os.path.isdir(os.path.join(case_dir, "mfem")):
            check(bool(glob.glob(os.path.join(case_dir, "mfem", "*.msh")) or glob.glob(os.path.join(case_dir, "mfem", "*.mesh"))),
                  "MFEM input present (mfem/*.msh)")
        for g in glob.glob(os.path.join(case_dir, "mfem", "*.msh")):
            check_mfem(g, g[:-4] + "_attributes.json")
            check_mfem_load(g, vol_total)
        for m in glob.glob(os.path.join(case_dir, "mfem", "*.mesh")):
            check_mfem(m, m[:-5] + "_attributes.json")
            check_mfem_load(m, vol_total)
        for g in [x for x in glob.glob(os.path.join(case_dir, "mfem", "*.msh")) if os.path.exists(x[:-4] + ".mesh")]:
            try:
                import mfem.ser as mfem
                a = mfem.Mesh(g[:-4] + ".mesh", 1, 1); b = mfem.Mesh(g, 1, 1)
                same = (a.GetNE(), a.GetNBE(), a.GetNV(), a.attributes.ToList(), a.bdr_attributes.ToList()) == \
                       (b.GetNE(), b.GetNBE(), b.GetNV(), b.attributes.ToList(), b.bdr_attributes.ToList())
                check(same and b.CheckElementOrientation(False) == 0, f"Gmsh {os.path.basename(g)}: MFEM reads it as the same mesh as the .mesh file")
            except ImportError:
                print("NOTE Gmsh check skipped (PyMFEM not installed)")
        if "--no-physics" not in sys.argv and os.path.exists(os.path.join(case_dir, "solver.xml")):
            import check_physics
            check_physics.run(case_dir, check, warn, lambda msg: print("NOTE " + msg))
    print("FAILURES:", fails)
    sys.exit(1 if fails else 0)
