#!/usr/bin/env python3
"""Download and unpack mesh archives from the suite's GitHub Release into cases/, checking SHA-256.

Each case level has two archives, so you fetch only what your solver reads:
  svmp  mesh/ (svMultiPhysics mesh-complete: mesh-complete.mesh.vtu + mesh-surfaces/*.vtp)
  mfem  mfem/<case>_<level>.msh (Gmsh 2.2 binary; MFEM reads it directly)
The small files (solver.xml, *.flow, case.json, quality reports, mfem/*_attributes.json) are
already in the git checkout.

    python tools/fetch_meshes.py --list                              # what is available, with sizes
    python tools/fetch_meshes.py 0031_H_ABAO_AAA/fine --solver mfem  # one level, MFEM file only
    python tools/fetch_meshes.py 0031_H_ABAO_AAA                     # one case, all levels, both parts
    python tools/fetch_meshes.py --level coarse --solver svmp        # every case's coarse level
    python tools/fetch_meshes.py --all                               # everything (about 3.6 GB)

No login is needed for the public repository. If GITHUB_TOKEN (or GH_TOKEN) is set, downloads go
through the GitHub REST API with it instead, which also works for a private fork or mirror.
Standard library only.
"""
import argparse, hashlib, json, os, sys, tarfile, urllib.error, urllib.request

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


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


def main():
    ap = argparse.ArgumentParser(description="Fetch benchmark meshes from the GitHub Release.")
    ap.add_argument("select", nargs="*", help="<case> or <case>/<level>")
    ap.add_argument("--level", help="only this level (coarse, medium, fine, published)")
    ap.add_argument("--solver", choices=("svmp", "mfem", "both"), default="both")
    ap.add_argument("--all", action="store_true", help="every case and level")
    ap.add_argument("--list", action="store_true", help="list levels and archive sizes, download nothing")
    ap.add_argument("--keep-archives", action="store_true")
    a = ap.parse_args()
    idx = json.load(open(os.path.join(ROOT, "meshes.json")))
    levels = idx["levels"]
    if a.list:
        for k, v in levels.items():
            print(f"{k:42s} {v['element_type']} {v['elements']:>10,}  " +
                  "  ".join(f"{p} {q['bytes']/1e6:5.0f} MB" for p, q in v["parts"].items()))
        return
    if not (a.select or a.level or a.all):
        ap.error("name a case or level, or use --level / --all (see --list)")
    todo = [k for k in levels if (a.all or not a.select or any(k == s or k.split("/")[0] == s for s in a.select))
            and (not a.level or k.endswith("/" + a.level))]
    if not todo:
        sys.exit("nothing matches; see --list")
    parts = ("svmp", "mfem") if a.solver == "both" else (a.solver,)
    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    url_for, accept = asset_urls(idx, token)
    dl = os.path.join(ROOT, ".downloads"); os.makedirs(dl, exist_ok=True)
    for k in todo:
        for p in parts:
            info = levels[k]["parts"][p]; arc = os.path.join(dl, info["archive"])
            if not os.path.exists(arc):
                print(f"downloading {k} [{p}] ({info['bytes']/1e6:.0f} MB)")
                try:
                    download(url_for(info["archive"]), arc, token, accept, info["bytes"])
                except urllib.error.HTTPError as e:
                    hint = " (if the repository is private, set GITHUB_TOKEN)" if e.code == 404 and not token else ""
                    sys.exit(f"download failed for {info['archive']}: HTTP {e.code}{hint}")
            h = hashlib.sha256()
            with open(arc, "rb") as f:
                for b in iter(lambda: f.read(1 << 20), b""):
                    h.update(b)
            if h.hexdigest() != info["sha256"]:
                os.remove(arc); sys.exit(f"checksum mismatch for {info['archive']}; removed the download, try again")
            with tarfile.open(arc) as t:
                t.extractall(ROOT, filter="data")
            if not a.keep_archives:
                os.remove(arc)
            print(f"ok {k} [{p}]")


if __name__ == "__main__":
    main()
