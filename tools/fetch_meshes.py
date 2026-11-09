#!/usr/bin/env python3
"""Download mesh files from the suite's GitHub Release into cases/, checking SHA-256.

The release has one zip per case (e.g. 0031_H_ABAO_AAA.zip) holding every level of that case.
This tool downloads the zip of each case you name and extracts only what you ask for:
  svmp  mesh/ (svMultiPhysics mesh-complete: mesh-complete.mesh.vtu + mesh-surfaces/*.vtp)
  mfem  mfem/<case>_<level>.msh (Gmsh 2.2 binary; MFEM reads it directly)
The small case files (solver.xml, *.flow, case.json, quality reports, mfem/*_attributes.json)
are already in the git checkout. Licence notices in the zip are extracted with the meshes.

    python tools/fetch_meshes.py --list                              # cases, levels, zip sizes
    python tools/fetch_meshes.py 0031_H_ABAO_AAA/fine --solver mfem  # one level, MFEM file only
    python tools/fetch_meshes.py 0031_H_ABAO_AAA                     # one case, all levels, both
    python tools/fetch_meshes.py --level coarse --solver svmp        # every case's coarse level
    python tools/fetch_meshes.py --all                               # everything (about 3.4 GB)

Downloaded zips are kept in .downloads/ (reused by later calls) unless --no-keep is given.
No login is needed for the public repository. If GITHUB_TOKEN (or GH_TOKEN) is set, downloads go
through the GitHub REST API with it instead, which also works for a private fork or mirror.
Standard library only.
"""
import argparse, hashlib, json, os, shutil, sys, urllib.error, urllib.request, zipfile

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
NOTICES = ("LICENSE", "DATA_LICENSE.md", "VMR-README-COPYRIGHT.txt", "svMultiPhysics-LICENSE.txt")


def _get(url, token=None, accept=None):
    req = urllib.request.Request(url, headers={"User-Agent": "cfd-benchmark-fetch"})
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    if accept:
        req.add_header("Accept", accept)
    return urllib.request.urlopen(req)


def asset_urls(idx, token):
    """Archive name -> download URL. Public: the release download URL. With a token: the API asset URL."""
    base = idx["base_url"].rstrip("/")
    if not token or not base.startswith("https://github.com/"):
        return lambda name: f"{base}/{name}", None
    api = f"https://api.github.com/repos/{idx['repository']}/releases/tags/{idx['tag']}"
    try:
        rel = json.load(_get(api, token, "application/vnd.github+json"))
    except urllib.error.HTTPError as e:
        sys.exit(f"release {idx['tag']} of {idx['repository']} not found (HTTP {e.code}); "
                 "check the tag in meshes.json and that the token can read the repository")
    assets = {a["name"]: a["url"] for a in rel["assets"]}
    return (lambda name: assets[name]), "application/octet-stream"


def download(url, dst, token, accept, size):
    with _get(url, token, accept) as r, open(dst + ".part", "wb") as f:
        got = 0
        for b in iter(lambda: r.read(1 << 20), b""):
            f.write(b); got += len(b)
            print(f"\r  {got/1e6:7.0f} / {size/1e6:.0f} MB", end="", flush=True)
    print()
    os.replace(dst + ".part", dst)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def wanted(name, case, levels, parts):
    """Is zip member `name` part of the selection? Members are cases/<case>/<level>/..."""
    p = name.split("/")
    if len(p) < 4 or p[0] != "cases" or p[1] != case or p[2] not in levels or ".." in p or name.endswith("/"):
        return False
    rest = p[3:]
    if rest[0] == "mesh":
        return "svmp" in parts
    if rest[0] == "mfem" and rest[-1].endswith(".msh"):
        return "mfem" in parts
    return len(rest) == 1 and rest[0] in NOTICES


def main():
    ap = argparse.ArgumentParser(description="Fetch benchmark meshes from the GitHub Release.")
    ap.add_argument("select", nargs="*", help="<case> or <case>/<level>")
    ap.add_argument("--level", help="only this level (coarse, medium, fine, published, ...)")
    ap.add_argument("--solver", choices=("svmp", "mfem", "both"), default="both")
    ap.add_argument("--all", action="store_true", help="every case and level")
    ap.add_argument("--list", action="store_true", help="list cases and levels, download nothing")
    ap.add_argument("--no-keep", action="store_true", help="delete each zip after extracting")
    a = ap.parse_args()
    idx = json.load(open(os.path.join(ROOT, "meshes.json")))
    cases = idx["cases"]
    if a.list:
        for c, v in cases.items():
            print(f"{c}  ({v['archive']}, {v['bytes']/1e6:.0f} MB)")
            for lv, info in v["levels"].items():
                print(f"    {lv:10s} {info['element_type']} {info['elements']:>10,}")
        return
    if not (a.select or a.level or a.all):
        ap.error("name a case or level, or use --level / --all (see --list)")
    plan = {}   # case -> levels
    for c, v in cases.items():
        for lv in v["levels"]:
            k = f"{c}/{lv}"
            if (a.all or not a.select or any(k == s or c == s for s in a.select)) and (not a.level or lv == a.level):
                plan.setdefault(c, set()).add(lv)
    if not plan:
        sys.exit("nothing matches; see --list")
    parts = ("svmp", "mfem") if a.solver == "both" else (a.solver,)
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    url_for, accept = asset_urls(idx, token)
    dl = os.path.join(ROOT, ".downloads"); os.makedirs(dl, exist_ok=True)
    for c, levels in plan.items():
        info = cases[c]; arc = os.path.join(dl, info["archive"])
        if os.path.exists(arc) and sha256(arc) != info["sha256"]:
            os.remove(arc)
        if not os.path.exists(arc):
            print(f"downloading {info['archive']} ({info['bytes']/1e6:.0f} MB)")
            try:
                download(url_for(info["archive"]), arc, token, accept, info["bytes"])
            except urllib.error.HTTPError as e:
                hint = " (if the repository is private, set GITHUB_TOKEN)" if e.code == 404 and not token else ""
                sys.exit(f"download failed for {info['archive']}: HTTP {e.code}{hint}")
            if sha256(arc) != info["sha256"]:
                os.remove(arc); sys.exit(f"checksum mismatch for {info['archive']}; removed the download, try again")
        n = 0
        with zipfile.ZipFile(arc) as z:
            for m in z.infolist():
                if wanted(m.filename, c, levels, parts):
                    dst = os.path.join(ROOT, *m.filename.split("/"))
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    with z.open(m) as src, open(dst, "wb") as out:
                        shutil.copyfileobj(src, out)
                    n += 1
        if a.no_keep:
            os.remove(arc)
        print(f"ok {c}: {', '.join(sorted(levels))} [{'+'.join(parts)}] ({n} files)")


if __name__ == "__main__":
    main()
