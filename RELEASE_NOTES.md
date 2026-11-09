# CFD Benchmark Suite v1.0

Open cardiovascular CFD cases for timing and checking GPU-accelerated solvers on identical meshes
and boundary conditions: five patient models from the Vascular Model Repository (AAA, two Fontan,
two aortas), a derived aorta and aortoiliac case, idealized AAA and Fontan, a tet and a hex pipe,
and the 1 M-tet pipe of svMultiPhysics issue #750. 39 case levels from about 1 M to 5.8 M linear
elements (tet4; hex8 for `pipe_hex`).

## Assets

One zip per case (`<case>.zip`, 17 to 373 MB), each holding every level of that case:

- `cases/<case>/<level>/mesh/`: svMultiPhysics mesh-complete (volume VTU and face VTPs)
- `cases/<case>/<level>/mfem/<case>_<level>.msh`: Gmsh 2.2 binary, read directly by MFEM
- the case files and the licence notices that apply to that case
- `SHA256SUMS`: checksums of the 12 zips

Unzip into a checkout of the repository, or let the downloader extract only what you need:

```sh
python tools/fetch_meshes.py --list
python tools/fetch_meshes.py 0031_H_ABAO_AAA/fine --solver mfem
```

Every level passes `tools/verify_case.py` (mesh and face consistency, quality gate, volume, BCs,
0D pressure, MFEM load and smoke solve). MFEM ASCII `.mesh` files are not shipped; the `.msh` is
the identical mesh, and `tools/msh_to_mfem.py` writes a `.mesh` on request.

## Licence and citation

Code: BSD 3-Clause (`LICENSE`). Data keep their sources' terms (`DATA_LICENSE.md`): VMR-derived cases under
the VMR README-COPYRIGHT (research and development use; notice included in every archive),
svMultiPhysics-derived cases under their BSD 3-Clause, synthetic cases under the suite's BSD 3-Clause. Please cite the suite
(`CITATION.cff`) and, for VMR cases, the VMR acknowledgement and Wilson et al. 2013,
doi:10.1115/1.4025983.
