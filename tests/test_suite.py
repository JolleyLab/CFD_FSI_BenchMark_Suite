#!/usr/bin/env python3
"""Fast checks for CI on every push (no mesh downloads; about a minute).

    python tests/test_suite.py            # all tests; also runs under pytest

1. Release index (meshes.json, SHA256SUMS): one zip per case, named <case>.zip, with a SHA-256
   that matches SHA256SUMS and a size under GitHub's 2 GiB asset limit; every level listed.
2. Case tree: every level in meshes.json has its case files in git (solver.xml, flow files it
   references, case.json, quality reports, attribute map, executable run.sh), every case folder
   in git is in meshes.json, and each case carries its source's licence notice (DATA_LICENSE.md).
3. Tools: every script in tools/ compiles; fetch_meshes.py --list runs.
4. Tiny case (tests/cases/pipe_tiny/tiny, about 15 k tets, committed): passes the full verifier
   (needs vtk, numpy and PyMFEM), and tools/msh_to_mfem.py writes an MFEM file that MFEM reads as
   the same mesh.
"""
import glob, json, os, py_compile, re, shutil, stat, subprocess, sys, tempfile

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
TINY = os.path.join(ROOT, "tests", "cases", "pipe_tiny", "tiny")
VMR = re.compile(r"^\d{4}_|^aorta_0159")


def _index():
    return json.load(open(os.path.join(ROOT, "meshes.json")))


def test_release_index():
    idx = _index()
    assert idx["base_url"] == f"https://github.com/{idx['repository']}/releases/download/{idx['tag']}"
    sums = dict(reversed(l.split()) for l in open(os.path.join(ROOT, "SHA256SUMS")).read().split("\n") if l.strip())
    for case, v in idx["cases"].items():
        assert v["archive"] == f"{case}.zip", v["archive"]
        assert sums.get(v["archive"]) == v["sha256"], f"{v['archive']}: SHA256SUMS mismatch"
        assert 0 < v["bytes"] < 2 * 1024**3, f"{v['archive']}: over GitHub's 2 GiB asset limit"
        assert v["levels"], case
        for lv, info in v["levels"].items():
            assert info["element_type"] in ("tet4", "hex8") and info["elements"] > 0, f"{case}/{lv}"
    assert set(sums) == {v["archive"] for v in idx["cases"].values()}, "SHA256SUMS lists exactly the case zips"


def _levels():
    return [f"{c}/{lv}" for c, v in _index()["cases"].items() for lv in v["levels"]]


def test_case_tree():
    levels = _levels()
    in_git = {os.path.relpath(os.path.dirname(p), os.path.join(ROOT, "cases"))
              for p in glob.glob(os.path.join(ROOT, "cases", "*", "*", "case.json"))}
    assert in_git == set(levels), f"case folders vs meshes.json: {sorted(in_git ^ set(levels))}"
    for level in levels:
        d = os.path.join(ROOT, "cases", level)
        for f in ("solver.xml", "case.json", "quality.json", "quality.md", "run.sh"):
            assert os.path.isfile(os.path.join(d, f)), f"{level}: {f}"
        assert os.stat(os.path.join(d, "run.sh")).st_mode & stat.S_IXUSR, f"{level}: run.sh not executable"
        assert glob.glob(os.path.join(d, "mfem", "*_attributes.json")), f"{level}: attribute map"
        xml = open(os.path.join(d, "solver.xml")).read()
        for ref in re.findall(r">\s*([^<>\s]+\.(?:flow|vtu|vtp))\s*<", xml):
            if ref.endswith(".flow"):
                assert os.path.isfile(os.path.join(d, ref)), f"{level}: {ref} referenced but missing"
            else:
                assert ref.startswith("mesh/"), f"{level}: mesh path {ref} outside mesh/"
        c = json.load(open(os.path.join(d, "case.json")))
        assert c.get("boundary_conditions"), f"{level}: no BCs in case.json"
        case = level.split("/")[0]
        if VMR.match(case):
            assert os.path.isfile(os.path.join(d, "VMR-README-COPYRIGHT.txt")), f"{level}: VMR notice"
        if case in ("aortoiliac_svmp_iliac", "stanford_reference"):
            assert os.path.isfile(os.path.join(d, "svMultiPhysics-LICENSE.txt")), f"{level}: svMultiPhysics BSD-3 notice"
    for f in ("LICENSE", "DATA_LICENSE.md", "CITATION.cff"):
        assert os.path.isfile(os.path.join(ROOT, f)), f


def test_tools_compile_and_list():
    for f in glob.glob(os.path.join(ROOT, "tools", "*.py")):
        py_compile.compile(f, doraise=True)
    out = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "fetch_meshes.py"), "--list"],
                         capture_output=True, text=True, check=True).stdout
    assert len(out.strip().split("\n")) == len(_levels()) + len(_index()["cases"])


def test_tiny_case_verifies():
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify_case.py"), TINY], capture_output=True, text=True)
    assert r.returncode == 0 and "FAILURES: 0" in r.stdout, r.stdout[-3000:] + r.stderr[-2000:]
    assert "MFEM smoke solve: flux balance" in r.stdout, "smoke solve did not run"


def test_msh_to_mfem_round_trip():
    import numpy as np
    import mfem.ser as mfem
    tmp = tempfile.mkdtemp()
    try:
        shutil.copytree(os.path.join(TINY, "mfem"), os.path.join(tmp, "mfem"))
        subprocess.run([sys.executable, os.path.join(ROOT, "tools", "msh_to_mfem.py"), tmp], check=True, capture_output=True)
        a = mfem.Mesh(glob.glob(os.path.join(tmp, "mfem", "*.msh"))[0], 1, 1)
        b = mfem.Mesh(glob.glob(os.path.join(tmp, "mfem", "*.mesh"))[0], 1, 1)
        assert (a.GetNE(), a.GetNBE(), a.GetNV()) == (b.GetNE(), b.GetNBE(), b.GetNV())
        assert np.abs(np.array(a.GetVertexArray()) - np.array(b.GetVertexArray())).max() == 0.0
        assert np.array_equal(np.array(a.GetBdrAttributeArray()), np.array(b.GetBdrAttributeArray()))
    finally:
        shutil.rmtree(tmp)


if __name__ == "__main__":
    failed = 0
    for name, fn in [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]:
        try:
            fn(); print(f"PASS {name}")
        except Exception as e:
            failed += 1; print(f"FAIL {name}: {e}")
    print(f"FAILURES: {failed}")
    sys.exit(1 if failed else 0)
