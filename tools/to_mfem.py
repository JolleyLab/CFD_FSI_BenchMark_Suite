"""Write an MFEM "mesh v1.0" file from a face-labelled VTK unstructured grid.

Volume elements keep their region id as the MFEM element attribute (fluid = 1, solid = 2 in the
FSI cases); boundary cells take the face id (ModelFaceID / CellEntityIds) as the boundary
attribute, so `attributes.json` beside the mesh maps every boundary attribute to the face name
used by the svMultiPhysics case (`mesh/mesh-surfaces/<name>.vtp`), and boundary conditions can
be bound by name on either side.

MFEM geometry codes: 2 = TRIANGLE, 3 = SQUARE, 4 = TETRAHEDRON, 5 = CUBE, 6 = PRISM.
Node ordering of VTK tet/hex/wedge matches MFEM's for these linear elements.

Credits: writes the MFEM mesh v1.0 format (Anderson et al. 2021, doi:10.1016/j.camwa.2020.06.009).
See REFERENCES.md.
"""
import json
import os
import sys

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

GEOM = {vtk.VTK_TRIANGLE: 2, vtk.VTK_QUAD: 3, vtk.VTK_TETRA: 4, vtk.VTK_HEXAHEDRON: 5, vtk.VTK_WEDGE: 6}
NODES = {2: 3, 3: 4, 4: 4, 5: 8, 6: 6}


def cell_arrays(mesh):
    types = np.array([mesh.GetCellType(i) for i in range(mesh.GetNumberOfCells())], dtype=np.int32)
    offsets = vtk_to_numpy(mesh.GetCells().GetOffsetsArray())
    conn = vtk_to_numpy(mesh.GetCells().GetConnectivityArray())
    return types, offsets, conn


def write_mfem(mesh, path, *, face_id_array="ModelFaceID", region_id_array="ModelRegionID",
               default_region=1, names_by_face_id=None, scale=1.0):
    types, offsets, conn = cell_arrays(mesh)
    conn = conn.copy()
    cd = mesh.GetCellData()
    face_ids = vtk_to_numpy(cd.GetArray(face_id_array)).astype(np.int64) if cd.GetArray(face_id_array) else np.zeros(len(types), np.int64)
    if cd.GetArray(region_id_array) is not None:
        region = vtk_to_numpy(cd.GetArray(region_id_array)).astype(np.int64)
    else:
        region = np.full(len(types), default_region, np.int64)
    region = np.where(region <= 0, default_region, region)

    is_vol = np.isin(types, [vtk.VTK_TETRA, vtk.VTK_HEXAHEDRON, vtk.VTK_WEDGE])
    is_bnd = np.isin(types, [vtk.VTK_TRIANGLE, vtk.VTK_QUAD])
    vol_idx = np.flatnonzero(is_vol)
    bnd_idx = np.flatnonzero(is_bnd)
    if bnd_idx.size == 0:
        raise ValueError("mesh carries no boundary cells; MFEM needs them for boundary attributes")
    bad = bnd_idx[face_ids[bnd_idx] <= 0]
    if bad.size:
        raise ValueError(f"{bad.size} boundary cells have face id <= 0; MFEM attributes must be >= 1")

    # Compact vertices to those used by volume elements (boundary cells share them).
    used = np.zeros(mesh.GetNumberOfPoints(), bool)
    for i in vol_idx:
        used[conn[offsets[i]:offsets[i + 1]]] = True
    new_id = -np.ones(mesh.GetNumberOfPoints(), np.int64)
    new_id[used] = np.arange(used.sum())
    pts = vtk_to_numpy(mesh.GetPoints().GetData())[used] * scale

    # positive tets and outward boundary faces (MFEM fixes both on load, but write them right)
    P = vtk_to_numpy(mesh.GetPoints().GetData())
    tet_idx = vol_idx[types[vol_idx] == vtk.VTK_TETRA]
    if tet_idx.size:
        T = np.stack([conn[offsets[tet_idx] + k] for k in range(4)], axis=1)
        a, b, c, d = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]], P[T[:, 3]]
        neg = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) < 0
        if neg.any():
            i1, i2 = offsets[tet_idx[neg]] + 1, offsets[tet_idx[neg]] + 2
            conn[i1], conn[i2] = conn[i2].copy(), conn[i1].copy()
            T = np.stack([conn[offsets[tet_idx] + k] for k in range(4)], axis=1)
        # owner lookup for boundary triangles via sorted-face keys
        nv = len(P); faces = np.concatenate([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]])
        faces.sort(axis=1)
        keys = (faces[:, 0] * nv + faces[:, 1]) * nv + faces[:, 2]
        owners = np.tile(np.arange(len(T)), 4)
        order = np.argsort(keys, kind="stable"); keys_s = keys[order]; owners_s = owners[order]
        tri_b = bnd_idx[types[bnd_idx] == vtk.VTK_TRIANGLE]
        if tri_b.size:
            B = np.stack([conn[offsets[tri_b] + k] for k in range(3)], axis=1)
            bs = np.sort(B, axis=1); bkeys = (bs[:, 0] * nv + bs[:, 1]) * nv + bs[:, 2]
            pos = np.searchsorted(keys_s, bkeys)
            found = (pos < len(keys_s)) & (keys_s[np.minimum(pos, len(keys_s) - 1)] == bkeys)
            own = T[owners_s[pos[found]]]
            cent = P[own].mean(axis=1)
            Bf = B[found]
            normal = np.cross(P[Bf[:, 1]] - P[Bf[:, 0]], P[Bf[:, 2]] - P[Bf[:, 0]])
            flip = np.einsum("ij,ij->i", normal, P[Bf].mean(axis=1) - cent) < 0
            fl = tri_b[found][flip]
            conn[offsets[fl] + 1], conn[offsets[fl] + 2] = conn[offsets[fl] + 2].copy(), conn[offsets[fl] + 1].copy()
    # hexes from the structured generator are already positive with outward quads

    with open(path, "w") as f:
        f.write("MFEM mesh v1.0\n\ndimension\n3\n\n")
        f.write(f"elements\n{len(vol_idx)}\n")
        for i in vol_idx:
            g = GEOM[types[i]]
            ids = new_id[conn[offsets[i]:offsets[i + 1]]]
            f.write(f"{region[i]} {g} " + " ".join(map(str, ids)) + "\n")
        f.write(f"\nboundary\n{len(bnd_idx)}\n")
        for i in bnd_idx:
            g = GEOM[types[i]]
            ids = new_id[conn[offsets[i]:offsets[i + 1]]]
            if (ids < 0).any():
                raise ValueError("boundary cell references a vertex not used by any volume element")
            f.write(f"{face_ids[i]} {g} " + " ".join(map(str, ids)) + "\n")
        f.write(f"\nvertices\n{len(pts)}\n3\n")
        np.savetxt(f, pts, fmt="%.10g")

    attrs = {
        "element_attributes": {int(r): ("fluid" if r == 1 else ("solid" if r == 2 else f"region_{r}"))
                               for r in np.unique(region[vol_idx])},
        "boundary_attributes": {int(fid): (names_by_face_id or {}).get(int(fid), f"face_{fid}")
                                for fid in np.unique(face_ids[bnd_idx])},
        "notes": (names_by_face_id or {}).get("note"),
        "length_unit": "cm" if scale == 1.0 else f"scaled_by_{scale}",
    }
    json.dump(attrs, open(str(path).rsplit(".", 1)[0] + "_attributes.json", "w"), indent=2)
    return attrs


def write_mfem_gmsh(mesh, msh_path, keep_mfem=False, **kw):
    """Write the release format for MFEM: Gmsh 2.2 binary at `msh_path` plus `<name>_attributes.json`.

    The MFEM v1.0 file is written first (it applies the orientation fixes above), converted to Gmsh
    by to_gmsh, then removed unless keep_mfem. MFEM reads the .msh as the identical mesh (same
    vertices, connectivity and attributes; checked by tools/compare_mfem_gmsh.py). To get the
    ASCII MFEM file back from a .msh, run tools/msh_to_mfem.py.
    """
    import to_gmsh
    mfem_path = str(msh_path)[:-4] + ".mesh"
    attrs = write_mfem(mesh, mfem_path, **kw)
    V, el, bd = to_gmsh.read_mfem(mfem_path)
    to_gmsh.write_gmsh22(str(msh_path), V, el, bd, attrs, binary=True)
    if not keep_mfem:
        os.remove(mfem_path)
    return attrs


def read_ug(path):
    r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(path); r.Update()
    return r.GetOutput()


if __name__ == "__main__":
    print(write_mfem(read_ug(sys.argv[1]), sys.argv[2]))
