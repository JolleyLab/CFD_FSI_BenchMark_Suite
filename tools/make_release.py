#!/usr/bin/env python3
"""Assemble the release (maintainers): a git-ready tree (everything except the large meshes) and one archive
per case level holding its meshes.

    make_release.py [--out release]

    release/CFD_FSI_BenchMark_Suite/          git tree: README, docs, tools, specs, surfaces, and every
                                              case level's solver.xml, *.flow, run.sh, case.json,
                                              quality reports and MFEM attribute maps (all small)
    release/meshes/<case>__<level>__svmp.tar.gz   mesh/ (svMultiPhysics mesh-complete: VTU + VTP faces)
    release/meshes/<case>__<level>__mfem.tar.gz   mfem/*.msh (Gmsh 2.2 binary; MFEM reads it directly)
    release/meshes/SHA256SUMS                     checksums of every archive
    release/bundles/<case>.zip                    (--bundles) one zip per case: every level's folder,
                                                  meshes and case files, for sharing outside GitHub

The archives are uploaded as GitHub Release assets of one tagged suite version (default v1.0),
so meshes.json points at https://github.com/JolleyLab/CFD_FSI_BenchMark_Suite/releases/download/<tag>.
GitHub limits: each asset under 2 GiB, up to 1000 assets per release, no total size limit
(docs.github.com, "About releases"); the script fails if an archive is too large.
`tools/fetch_meshes.py` in the tree downloads and unpacks any subset of them into place.

Publishing from a checkout: unpack the meshes into cases/ (fetch_meshes.py --all, or the per-case
bundle zips), run `python tools/make_release.py --out <dir>`, upload <dir>/meshes/* to the
release with `gh release upload <tag> ...`, and commit the meshes.json and SHA256SUMS it writes
(checksums change whenever archives are rebuilt).
The MFEM ASCII .mesh files are not shipped: the .msh carries the identical mesh
(tools/compare_mfem_gmsh.py) and tools/msh_to_mfem.py writes .mesh files on demand.
ALE/FSI folders are left out of this release (deferred).
"""
import argparse, glob, hashlib, json, os, shutil, sys, tarfile

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
EXCLUDE_CASES = {"pipe_fsi_hex"}
TOP_FILES = ["LICENSE", "DATA_LICENSE.md", "CITATION.cff", "README.md", "CATALOG.md", "REFERENCES.md", "SUMMARY.md", "RELEASE_NOTES.md", "requirements.txt"]
REPO = "JolleyLab/CFD_FSI_BenchMark_Suite"
MAX_ASSET = 2 * 1024**3 - 1          # GitHub release asset limit (each file under 2 GiB)
MAX_ASSETS = 1000                    # per release
SMALL_PATTERNS = ["solver.xml", "*.flow", "run.sh", "case.json", "quality.md", "quality.json", "mfem/*_attributes.json"]


def case_notices(case):
    """Source licence notices that must travel with a case's data (see DATA_LICENSE.md)."""
    if case.startswith("0") or case == "aorta_0159_H_AO_H":
        return ["VMR-README-COPYRIGHT.txt"]                 # VMR README-COPYRIGHT
    if case in ("aortoiliac_svmp_iliac", "stanford_reference"):
        return ["svMultiPhysics-LICENSE.txt"]               # BSD 3-Clause
    return []


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def index_bundles(bundles_dir, tag, levels):
    """meshes.json for per-case zip assets: case -> archive, bytes, sha256, and its levels."""
    cases = {}
    for z in sorted(glob.glob(os.path.join(bundles_dir, "*.zip"))):
        case = os.path.basename(z)[:-4]
        cases[case] = {"archive": os.path.basename(z), "bytes": os.path.getsize(z), "sha256": sha256(z),
                       "levels": {k.split("/")[1]: {kk: v[kk] for kk in ("element_type", "elements", "points")}
                                  for k, v in levels.items() if k.split("/")[0] == case}}
    missing = {k.split("/")[0] for k in levels} - set(cases)
    if missing:
        sys.exit(f"bundles missing for {sorted(missing)}; run with --bundles first")
    return {"repository": REPO, "tag": tag, "base_url": f"https://github.com/{REPO}/releases/download/{tag}",
            "note": "one zip per case holding every level: mesh/ (svMultiPhysics mesh-complete) and "
                    "mfem/<case>_<level>.msh (Gmsh 2.2 for MFEM); tools/fetch_meshes.py extracts the levels and solver you ask for",
            "cases": cases}


def write_bundles(out, cases):
    """One zip per case holding cases/<case>/<level>/ exactly as the release tree plus fetched meshes."""
    import zipfile
    os.makedirs(out, exist_ok=True); sums = []
    for case in cases:
        z = os.path.join(out, f"{case}.zip")
        if not os.path.exists(z):
            with zipfile.ZipFile(z + ".tmp", "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
                for d in sorted(glob.glob(os.path.join(ROOT, "cases", case, "*"))):
                    if not os.path.isfile(os.path.join(d, "case.json")):
                        continue
                    level = os.path.basename(d)
                    for dirpath, _, files in os.walk(d):
                        for f in files:
                            src = os.path.join(dirpath, f); rel = os.path.relpath(src, d)
                            if f.endswith(".mesh") or f.endswith(".log") or rel.split(os.sep)[0] in ("fsi",) or "-procs" in rel:
                                continue
                            zf.write(src, f"cases/{case}/{level}/{rel}")
                    for n in ("LICENSE", "DATA_LICENSE.md"):
                        zf.write(os.path.join(ROOT, n), f"cases/{case}/{level}/{n}")
                    for n in case_notices(case):
                        zf.write(os.path.join(ROOT, "surfaces", n), f"cases/{case}/{level}/{n}")
            os.replace(z + ".tmp", z)
        sums.append(f"{sha256(z)}  {os.path.basename(z)}")
        print(f"bundle {case}: {os.path.getsize(z)/1e6:.0f} MB", flush=True)
    open(os.path.join(out, "SHA256SUMS"), "w").write("\n".join(sums) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "release"))
    ap.add_argument("--tag", default="v1.0", help="release tag the archives are attached to")
    ap.add_argument("--bundles", action="store_true", help="also write one zip per case (release/bundles)")
    ap.add_argument("--assets", choices=("bundles", "levels"), default="bundles",
                    help="what the release assets are: per-case zips (v1.0) or per-level, per-solver tar.gz")
    a = ap.parse_args()
    tree = os.path.join(a.out, "CFD_FSI_BenchMark_Suite"); meshes = os.path.join(a.out, "meshes")
    if os.path.exists(tree):
        shutil.rmtree(tree)
    os.makedirs(tree); os.makedirs(meshes, exist_ok=True)
    for f in TOP_FILES:
        shutil.copy(os.path.join(ROOT, f), tree)
    shutil.copytree(os.path.join(ROOT, ".github"), os.path.join(tree, ".github"))
    os.makedirs(os.path.join(tree, "tools"))
    for f in sorted(glob.glob(os.path.join(ROOT, "tools", "*.py"))):
        if os.path.basename(f) in ("trim_branches.py", "wall_extrude.py", "fsi_mfem.py", "package.py", "relabel.py"):
            continue   # deferred (FSI) or build-housekeeping only
        shutil.copy(f, os.path.join(tree, "tools"))
    shutil.copytree(os.path.join(ROOT, "specs"), os.path.join(tree, "specs"))
    shutil.copytree(os.path.join(ROOT, "tests"), os.path.join(tree, "tests"), ignore=shutil.ignore_patterns("__pycache__", "*.log", "*-procs"))
    shutil.copytree(os.path.join(ROOT, "surfaces"), os.path.join(tree, "surfaces"),
                    ignore=shutil.ignore_patterns("*_centerlines.vtp"))
    index = {}; sums = []
    for d in sorted(glob.glob(os.path.join(ROOT, "cases", "*", "*"))):
        case, level = d.split(os.sep)[-2:]
        if case in EXCLUDE_CASES or not os.path.isfile(os.path.join(d, "case.json")):
            continue
        dst = os.path.join(tree, "cases", case, level); os.makedirs(os.path.join(dst, "mfem"), exist_ok=True)
        for pat in SMALL_PATTERNS:
            for f in glob.glob(os.path.join(d, pat)):
                shutil.copy(f, os.path.join(dst, os.path.relpath(f, d)))
        for n in case_notices(case):   # source licence notices are kept in git with each case level
            shutil.copy(os.path.join(ROOT, "surfaces", n), dst)
        c = json.load(open(os.path.join(d, "case.json")))
        e = c["elements"]
        notices = [(os.path.join(ROOT, "LICENSE"), f"cases/{case}/{level}/LICENSE"),
                   (os.path.join(ROOT, "DATA_LICENSE.md"), f"cases/{case}/{level}/DATA_LICENSE.md")]
        notices += [(os.path.join(ROOT, "surfaces", n), f"cases/{case}/{level}/{n}") for n in case_notices(case)]
        parts = {"svmp": [(os.path.join(d, "mesh"), f"cases/{case}/{level}/mesh")] + notices,
                 "mfem": [(f, f"cases/{case}/{level}/mfem/{os.path.basename(f)}") for f in sorted(glob.glob(os.path.join(d, "mfem", "*.msh")))] + notices}
        if not parts["mfem"] or not os.path.isfile(os.path.join(d, "mesh", "mesh-complete.mesh.vtu")):
            sys.exit(f"{case}/{level}: meshes missing; unpack them first (fetch_meshes.py --all, or the bundle zips)")
        entry = {"element_type": "hex8" if e.get("hex8") else "tet4", "elements": e.get("hex8") or e.get("tet4"),
                 "points": e.get("points"), "parts": {}}
        for part, members in parts.items():
            arc = os.path.join(meshes, f"{case}__{level}__{part}.tar.gz")
            if not os.path.exists(arc):
                with tarfile.open(arc + ".tmp", "w:gz", compresslevel=6) as t:
                    for src, name in members:
                        t.add(src, arcname=name)
                os.replace(arc + ".tmp", arc)
            size = os.path.getsize(arc)
            if size > MAX_ASSET:
                sys.exit(f"{arc} is {size/2**30:.2f} GiB, over GitHub's per-asset limit; split it")
            entry["parts"][part] = {"archive": os.path.basename(arc), "bytes": size, "sha256": sha256(arc)}
            sums.append(f"{entry['parts'][part]['sha256']}  {os.path.basename(arc)}")
        index[f"{case}/{level}"] = entry
        print(f"{case}/{level}: " + ", ".join(f"{p} {v['bytes']/1e6:.0f} MB" for p, v in entry["parts"].items()), flush=True)
    open(os.path.join(meshes, "SHA256SUMS"), "w").write("\n".join(sums) + "\n")
    if a.bundles:
        write_bundles(os.path.join(a.out, "bundles"), sorted({k.split("/")[0] for k in index}))
    if len(sums) + 1 > MAX_ASSETS:
        sys.exit(f"{len(sums)} archives exceed GitHub's {MAX_ASSETS} assets per release")
    if a.assets == "bundles":
        idx = index_bundles(os.path.join(a.out, "bundles"), a.tag, index)
        shutil.copy(os.path.join(a.out, "bundles", "SHA256SUMS"), tree)
    else:
        shutil.copy(os.path.join(meshes, "SHA256SUMS"), tree)
        idx = {"repository": REPO, "tag": a.tag, "base_url": f"https://github.com/{REPO}/releases/download/{a.tag}",
               "note": "one archive per case level and solver: <case>__<level>__svmp.tar.gz and <case>__<level>__mfem.tar.gz",
               "levels": index}
    json.dump(idx, open(os.path.join(tree, "meshes.json"), "w"), indent=2)
    open(os.path.join(tree, ".gitignore"), "w").write("# meshes are fetched, not committed\ncases/*/*/mesh/\ncases/*/*/mfem/*.mesh\ncases/*/*/mfem/*.msh\ncases/*/*/LICENSE\ncases/*/*/DATA_LICENSE.md\n*.tar.gz\n*-procs/\nsolver.run.xml\nrun_*.log\n__pycache__/\n.downloads/\n")
    total = sum(p["bytes"] for v in index.values() for p in v["parts"].values())
    print(f"{len(index)} levels, {len(sums)} archives, {total/1e9:.2f} GB; tree at {tree}")
    up = f"{a.out}/bundles/*.zip {a.out}/bundles/SHA256SUMS" if a.assets == "bundles" else f"{meshes}/*.tar.gz {meshes}/SHA256SUMS"
    print(f"upload: gh release upload {a.tag} {up} --repo {REPO}")


if __name__ == "__main__":
    main()
