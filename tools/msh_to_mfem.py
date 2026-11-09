#!/usr/bin/env python3
"""Write MFEM mesh v1.0 (ASCII) files from the released Gmsh .msh files, if a reader needs them.

The release ships mfem/<case>_<level>.msh (Gmsh 2.2 binary) only. MFEM reads it directly, e.g.
`mfem::Mesh mesh("0031_H_ABAO_AAA_fine.msh");`, with element attribute 1 = fluid and boundary
attribute = face id (names in *_attributes.json), so this conversion is optional. The output is
what MFEM's own Mesh::Print writes for the loaded mesh (PyMFEM, Anderson et al. 2021,
doi:10.1016/j.camwa.2020.06.009).

    python tools/msh_to_mfem.py cases/0031_H_ABAO_AAA/fine          # writes mfem/*.mesh beside the .msh
    python tools/msh_to_mfem.py --gzip cases/*/coarse                # writes mfem/*.mesh.gz (MFEM reads it)
"""
import argparse, glob, gzip, os, shutil


def convert(msh, gz=False):
    import mfem.ser as mfem
    m = mfem.Mesh(msh, 1, 1)
    out = msh[:-4] + ".mesh"
    m.Print(out, 16)
    if gz:
        with open(out, "rb") as a, gzip.open(out + ".gz", "wb") as b:
            shutil.copyfileobj(a, b)
        os.remove(out); out += ".gz"
    return out, m.GetNE()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cases", nargs="+", help="case level folders (cases/<case>/<level>) or .msh files")
    ap.add_argument("--gzip", action="store_true")
    a = ap.parse_args()
    for c in a.cases:
        for msh in ([c] if c.endswith(".msh") else sorted(glob.glob(os.path.join(c, "mfem", "*.msh")))):
            out, ne = convert(msh, a.gzip)
            print(f"wrote {out} ({ne:,} elements)")
