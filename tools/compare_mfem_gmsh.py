#!/usr/bin/env python3
"""Check that mfem/<name>.msh (Gmsh 2.2) carries exactly what mfem/<name>.mesh (MFEM v1.0) does.

Both files are loaded by MFEM itself (PyMFEM, Anderson et al. 2021, doi:10.1016/j.camwa.2020.06.009)
and compared array by array: vertex coordinates, element connectivity and attributes, boundary
element connectivity and attributes. Run before dropping the .mesh files from a release.

    python tools/compare_mfem_gmsh.py cases/*/*/
"""
import glob, os, sys
import numpy as np
import mfem.ser as mfem


def arrays(path):
    m = mfem.Mesh(path, 1, 1)
    v = np.array(m.GetVertexArray())
    e = np.array([m.GetElementVertices(i) for i in range(m.GetNE())], dtype=np.int64)
    ea = np.array(m.GetAttributeArray())
    b = np.array([m.GetBdrElementVertices(i) for i in range(m.GetNBE())], dtype=np.int64)
    ba = np.array(m.GetBdrAttributeArray())
    geo = (m.GetElementBaseGeometry(0), m.GetBdrElementBaseGeometry(0))
    ori = (m.CheckElementOrientation(False), m.CheckBdrElementOrientation(False))
    return dict(v=v, e=e, ea=ea, b=b, ba=ba, geo=geo, ori=ori)


def compare(mesh_path, msh_path):
    a, g = arrays(mesh_path), arrays(msh_path)
    out = {}
    out["vertices"] = a["v"].shape == g["v"].shape and float(np.abs(a["v"] - g["v"]).max()) == 0.0
    out["elements"] = np.array_equal(a["e"], g["e"]) and np.array_equal(a["ea"], g["ea"])
    out["boundary"] = np.array_equal(a["b"], g["b"]) and np.array_equal(a["ba"], g["ba"])
    out["geometry"] = a["geo"] == g["geo"]
    out["orientation"] = a["ori"] == g["ori"] == (0, 0)
    return out, a["e"].shape[0]


if __name__ == "__main__":
    bad = 0
    for d in sys.argv[1:]:
        for mp in sorted(glob.glob(os.path.join(d, "mfem", "*.mesh"))):
            gp = mp[:-5] + ".msh"
            if not os.path.exists(gp):
                print(f"MISSING {gp}"); bad += 1; continue
            res, ne = compare(mp, gp)
            ok = all(res.values())
            bad += not ok
            print(f"{'SAME' if ok else 'DIFF'} {mp}  ({ne:,} elements)  " + " ".join(f"{k}={'ok' if v else 'NO'}" for k, v in res.items()), flush=True)
    print(f"DIFFERENCES: {bad}")
    sys.exit(1 if bad else 0)
