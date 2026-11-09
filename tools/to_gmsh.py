#!/usr/bin/env python3
"""Write Gmsh .msh for a case from its MFEM mesh, so both files describe the identical mesh.

    to_gmsh.py cases/<case>/<level>      # writes mfem/<case>_<level>.msh (Gmsh 2.2 binary)
    to_gmsh.py --ascii cases/...         # Gmsh 2.2 ASCII

Physical groups carry the MFEM attributes: 3D group 1 = "fluid" (2 = "solid" in FSI meshes),
2D groups = boundary attributes, named as in mfem/*_attributes.json (cap_inflow, cap_<outlet>,
wall). Format 2.2 is the one every Gmsh-era reader takes (MFEM, Gmsh, FEniCS/meshio, deal.II);
MFEM takes element and boundary attributes from the physical tags. Binary is the default
(about 3x smaller than ASCII and much faster to read).

Credits: Gmsh file format (Geuzaine & Remacle 2009, doi:10.1002/nme.2579). See REFERENCES.md.
"""
import argparse, glob, json, os, sys
import numpy as np

GEOM = {4: ("tetra", 4), 5: ("hexahedron", 8), 2: ("triangle", 3), 3: ("quad", 4), 6: ("wedge", 6)}


def read_mfem(path):
    """Parse MFEM mesh v1.0 (linear elements) into points, element blocks and attributes."""
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    def seek(word):
        nonlocal i
        while lines[i].strip() != word:
            i += 1
        i += 1
    seek("elements"); ne = int(lines[i]); i += 1
    el = [list(map(int, lines[i + k].split())) for k in range(ne)]; i += ne
    seek("boundary"); nb = int(lines[i]); i += 1
    bd = [list(map(int, lines[i + k].split())) for k in range(nb)]; i += nb
    seek("vertices"); nv = int(lines[i]); i += 1; dim = int(lines[i]); i += 1
    V = np.loadtxt(lines[i:i + nv]).reshape(nv, dim)
    return V, el, bd


GMSH_TYPE = {4: (4, 4), 5: (5, 8), 2: (2, 3), 3: (3, 4), 6: (6, 6)}   # MFEM geom -> (gmsh elm-type, nodes)


def write_gmsh22(dst, V, el, bd, attrs, binary=True):
    """Gmsh 2.2 (ASCII or binary). Each element carries two tags: physical = MFEM attribute,
    elementary = the same; nodes are 1-based. Element blocks are grouped by type for binary."""
    import struct
    rows = [(r, 3) for r in el] + [(r, 2) for r in bd]
    with open(dst, "wb") as f:
        w = lambda t: f.write(t.encode())
        w("$MeshFormat\n2.2 %d 8\n" % (1 if binary else 0))
        if binary:
            f.write(struct.pack("<i", 1)); w("\n")
        w("$EndMeshFormat\n$PhysicalNames\n")
        names = [(3, int(k), v) for k, v in attrs["element_attributes"].items()] + \
                [(2, int(k), v) for k, v in attrs["boundary_attributes"].items()]
        w("%d\n" % len(names))
        for d, k, v in names:
            w('%d %d "%s"\n' % (d, k, v))
        w("$EndPhysicalNames\n$Nodes\n%d\n" % len(V))
        if binary:
            rec = np.zeros(len(V), dtype=[("id", "<i4"), ("x", "<f8", (3,))])
            rec["id"] = np.arange(1, len(V) + 1); rec["x"][:, :V.shape[1]] = V
            f.write(rec.tobytes()); w("\n")
        else:
            for i, p in enumerate(V):
                w("%d %.17g %.17g %.17g\n" % (i + 1, p[0], p[1], p[2] if len(p) > 2 else 0.0))
        w("$EndNodes\n$Elements\n%d\n" % len(rows))
        eid = 1
        if binary:
            by = {}
            for r, _ in rows:
                by.setdefault(r[1], []).append(r)
            for g, rs in by.items():
                gt, n = GMSH_TYPE[g]
                a = np.array(rs, dtype=np.int64)
                f.write(struct.pack("<iii", gt, len(a), 2))
                blk = np.empty((len(a), 3 + n), "<i4")
                blk[:, 0] = np.arange(eid, eid + len(a)); blk[:, 1] = a[:, 0]; blk[:, 2] = a[:, 0]
                blk[:, 3:] = a[:, 2:2 + n] + 1
                f.write(blk.tobytes()); eid += len(a)
            w("\n")
        else:
            for r, _ in rows:
                gt, n = GMSH_TYPE[r[1]]
                w("%d %d 2 %d %d %s\n" % (eid, gt, r[0], r[0], " ".join(str(x + 1) for x in r[2:2 + n]))); eid += 1
        w("$EndElements\n")


def write(case_dir, fmt="gmsh22", binary=True):
    out = []
    for mpath in glob.glob(os.path.join(case_dir, "mfem", "*.mesh")):
        attrs = json.load(open(mpath[:-5] + "_attributes.json"))
        V, el, bd = read_mfem(mpath)
        dst = mpath[:-5] + ".msh"
        write_gmsh22(dst, V, el, bd, attrs, binary=binary)
        out.append(dst)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--ascii", action="store_true")
    a = ap.parse_args()
    for c in a.cases:
        for p in write(c, binary=not a.ascii):
            print("wrote", p, f"{os.path.getsize(p)/1e6:.1f} MB")
