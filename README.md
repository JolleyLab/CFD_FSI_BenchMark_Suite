# CFD Benchmark Suite

A standard, open set of cardiovascular CFD cases for testing and comparing CFD frameworks on the
same data. Each case fixes everything except the solver: geometry, meshes at several sizes,
boundary conditions, fluid properties and units. Results and timings from different codes,
hardware and institutions can then be compared directly.

The cases are patient models from the [Vascular Model Repository](https://www.vascularmodel.com)
(VMR), aortas, an abdominal aortic aneurysm and Fontan circulations, with the boundary conditions
published with them, plus synthetic cases (idealized AAA and Fontan, straight pipes) for
verification. Meshes range from about 1 M to 5.8 M linear elements. Every case is ready to run in
svMultiPhysics, and the meshes come in two neutral formats, VTK XML and Gmsh 2.2, which MFEM
and Gmsh read directly and most other finite element and finite volume codes read through Gmsh
or meshio.
Everything is regenerable from the scripts in `tools/`.

Typical uses:

- **Performance:** time per step, strong and weak scaling, CPU against GPU, on identical meshes
  from about 1 M to 5.8 M elements.
- **Verification and comparison:** run the same case in two frameworks and compare pressures,
  flows and wall shear stress; three mesh levels per case support a grid convergence study.
- **Regression testing:** the verifier and the small test case run in CI.

Jolley Lab, Children's Hospital of Philadelphia, 2026. Meshes built with the SlicerHeart and
SimVascular Slicer modules run headless; see [Credits](#credits) and `REFERENCES.md`.

## Quick start

```sh
git clone https://github.com/JolleyLab/CFD_FSI_BenchMark_Suite.git && cd CFD_FSI_BenchMark_Suite
python tools/fetch_meshes.py --list                              # cases, levels, element counts, zip sizes
python tools/fetch_meshes.py 0031_H_ABAO_AAA/medium --solver svmp
cd cases/0031_H_ABAO_AAA/medium && ./run.sh 8 20                # svMultiPhysics: 8 MPI ranks, 20 time steps
```

Meshes are not stored in git. They are assets of the
[GitHub Release](https://github.com/JolleyLab/CFD_FSI_BenchMark_Suite/releases) (v1.0): one zip per
case (17 to 373 MB) holding every level. `fetch_meshes.py` downloads the zip for each case you
name and extracts only the levels and format you ask for: `--solver svmp` (VTK XML
mesh-complete) or `--solver mfem` (Gmsh `.msh`). It needs only the Python standard library,
checks each zip's SHA-256 (`meshes.json`, `SHA256SUMS`), keeps the zip in `.downloads/` for
reuse, and needs no GitHub login. `pip install -r requirements.txt` is needed only for the
verifier and the build tools.

## Using a case with your framework

Each case level, `cases/<case>/<level>/`, gives you:

- **Mesh:** `mesh/mesh-complete.mesh.vtu` (volume, VTK XML) with one `mesh/mesh-surfaces/<face>.vtp`
  per boundary face, or `mfem/<case>_<level>.msh` (Gmsh 2.2 binary) with the same mesh. Element
  attribute 1 is the fluid; each boundary face has an integer id, named in
  `mfem/*_attributes.json` and in the Gmsh physical names (`cap_inflow`, `cap_<outlet>`, `wall`).
- **Boundary conditions:** in `case.json`, by face name. A periodic inflow waveform (also as
  `<cap>.flow`, time and flow rate; negative flow enters the domain, as in svMultiPhysics), a
  three-element Windkessel (RCR) or resistance on each outlet, and rigid no-slip walls.
- **Fluid and units:** CGS throughout (cm, g, s); density 1.06 g/cm³, viscosity 0.04 P.
- **Reference setup:** `solver.xml` is the svMultiPhysics input. It is a complete, readable record
  of the time step, number of steps and BC types to mirror in another code.

**MFEM:** `mfem::Mesh mesh("0031_H_ABAO_AAA_medium.msh");` loads the mesh with its attributes.
For MFEM's own ASCII format, `python tools/msh_to_mfem.py cases/<case>/<level>` writes
`mfem/<case>_<level>.mesh` (`--gzip` for `.mesh.gz`); it is the same mesh, array for array
(checked by `tools/compare_mfem_gmsh.py` on every level before release).

**svMultiPhysics:** `run.sh` uses `$SVMULTIPHYSICS`, else `svmultiphysics` on the PATH, else the
official Docker image `simvascular/solver:latest`, and writes the wall time to
`run_<ranks>x<steps>.log`.

**Check a case:** `python tools/verify_case.py cases/<case>/<level>`.

## Reporting results

So that results line up across platforms and institutions, report with each run:

- suite version and case level (for example `v1.0, 0031_H_ABAO_AAA/fine`) and the element count;
- framework and version, discretization (element order, stabilization), linear solver and
  preconditioner, tolerances, time step and number of steps;
- hardware (CPU or GPU model and count, nodes, MPI ranks) and wall time per time step, excluding
  setup, with setup time reported separately;
- for comparisons of results: outlet flows and pressures over the last cycle, and wall shear
  stress on the named faces.

## Cases

All elements are linear: **tet4** (4-node tetrahedra) or **hex8** (8-node hexahedra).

| case | elements | anatomy | source | published | coarse\* | medium | fine |
|---|---|---|---|---:|---:|---:|---:|
| `0031_H_ABAO_AAA` | tet4 | abdominal aortic aneurysm, M 73, CT; 1 inlet, 9 outlets (celiac, SMA, renals, internal/external iliacs) | VMR | 2.98 M | 1.00 M | 2.40 M | 4.97 M |
| `0054_H_PULMFON_TAT` | tet4 | Fontan (TCPC), M 3 y, MR; 1 inflow, 20 pulmonary outlets | VMR | 1.01 M | 1.85 M | 3.83 M | 4.78 M |
| `0065_H_PULMFON_SVD` | tet4 | Fontan (TCPC), M 9 y, MR; 1 inflow, 31 pulmonary outlets | VMR | 1.06 M | 1.74 M | 2.51 M | 5.55 M |
| `0001_H_AO_SVD` | tet4 | aorta of a Fontan patient, M 3 y, MR | VMR | 2.06 M | 1.05 M | 2.50 M | 5.16 M |
| `0011_H_AO_H` | tet4 | healthy adult thoracic aorta, F 23, MR | VMR | 2.99 M | 1.26 M | 2.83 M | 3.99 M |
| `aorta_0159_H_AO_H` | tet4 | adult thoracic aorta, CT | VMR surface; BCs derived here | | 0.96 M | 3.52 M | 4.40 M |
| `aortoiliac_svmp_iliac` | tet4 | distal aorta and iliacs (Stanford's 76 k "iliac" case, remeshed) | svMultiPhysics test geometry | | 1.16 M | 2.64 M | 5.16 M |
| `aaa_idealized` | tet4 | idealized infrarenal AAA with iliac bifurcation | synthetic | | 0.95 M | 2.94 M | 4.24 M |
| `tcpc_idealized` | tet4 | idealized extracardiac Fontan | synthetic | | 1.10 M | 3.66 M | 5.75 M |
| `pipe` | tet4 | straight pipe, d 2 cm | synthetic | | 1.31 M | 3.85 M | 5.54 M |
| `pipe_hex` | hex8 | straight pipe, d 2 cm, structured O-grid | synthetic | | 1.45 M | 2.83 M | 4.67 M |
| `stanford_reference/pipe_1M` | tet4 | the 1 M-tet pipe of svMultiPhysics issue #750, regenerated bit-identically | `make_pipe.py` (S. Dave, A. Lasso) | 1.00 M | | | |

**published**: the VMR project's own mesh, the one its boundary conditions were tuned on,
repackaged unchanged. It includes wall boundary layers (thin elements at the wall); the
generated levels are isotropic and come from one recipe, so they form a clean size series.

\* **coarse** levels suit performance and smoke tests; their smallest outlets (renals, distal
pulmonary branches) have about 4.5 elements across. Use medium, fine or published levels for
accuracy comparisons. Generated counts land within about ±35 % of target because Netgen's
size field is quantized; levels are always ordered by count.

Stanford's issue #750 AAA is 0031 trimmed to the inlet and two iliacs, about 1.5 M tets
([svMultiPhysics paper](https://arxiv.org/abs/2607.19631)); this suite keeps the full model.
`CATALOG.md` lists further VMR candidates; `SUMMARY.md` lists every level with its quality metrics.

## Each case level

```
cases/<case>/<level>/
  solver.xml, <cap>.flow        svMultiPhysics rigid-wall CFD input; BCs bound by face name
  run.sh                        one-command run (see Quick start)
  case.json                     provenance, licence, purpose, element counts, edge length, every BC value
  quality.md / quality.json     element quality report
  mfem/*_attributes.json        boundary attribute -> face name
  VMR-README-COPYRIGHT.txt      VMR cases only: the licence notice that travels with the data
  -- extracted by tools/fetch_meshes.py --solver svmp --
  mesh/mesh-complete.mesh.vtu   volume mesh (GlobalNodeID, GlobalElementID, ModelRegionID)
  mesh/mesh-complete.exterior.vtp, mesh/walls_combined.vtp, mesh/mesh-surfaces/<face>.vtp
  -- extracted by tools/fetch_meshes.py --solver mfem --
  mfem/<case>_<level>.msh       Gmsh 2.2 binary, physical groups = MFEM attributes
```

Units are CGS throughout: cm, g, s; density 1.06 g/cm³, viscosity 0.04 P; pressure in dyn/cm²
(1 mmHg = 1333.22 dyn/cm²), resistance in dyn·s/cm⁵, compliance in cm⁵/dyn.

## Boundary conditions

- **VMR cases** use the project's own tuned conditions unchanged: the periodic inflow waveform
  and period from its SimVascular job file (Dirichlet, parabolic profile, imposed flux) and a
  three-element Windkessel (RCR) on every outlet.
- **Derived and synthetic cases** use an analytic aortic waveform (or a low-pulsatility caval
  waveform for the Fontan) scaled to a stated mean flow, with RCR outlets tuned from a
  mean-pressure target and a flow split by outlet area (`R_i = ΔP / (f_i Q)`, `Rp` 6–9 % of
  `R_i`, compliance split by flow), or pulmonary resistances from a PVR in Wood units.
- **Walls** are rigid and no-slip.

The flow files follow the svMultiPhysics convention (negative = into the domain).

## Tests

`tools/verify_case.py` reads only the written files, independent of the build code, and prints
PASS / WARN / FAIL per check. Every level in this release passes.

| check | what it catches |
|---|---|
| mesh-complete consistency | one element type; `GlobalNodeID` 1..N; face files map onto volume nodes; closed, manifold exterior; face areas sum to exterior area |
| quality gate | inverted elements; dihedral angles outside 2–178°; more than 0.5 % of elements with scaled Jacobian < 0.2 (WARN); boundary-layer (published) meshes checked for inversion only |
| volume | element volume equals the volume enclosed by the exterior (holes, overlaps) |
| cap resolution | fewer than 10 elements across an inflow cap or 5 across an outlet (WARN) |
| boundary conditions | every cap has a condition; inflow enters the domain; flow-file header and period; RCR values positive |
| 0D pressure | mean pressure from inflow and outlet resistances within 50–130 mmHg systemic, 5–25 mmHg pulmonary (WARN) |
| flow regime, time step | peak inlet Re (WARN above 4000), Womersley number, inlet CFL for the set time step |
| caps | planar; facing out of the domain |
| MFEM | `.msh` structure, physical names and attributes match `*_attributes.json`; every boundary face owned by one element; loads in PyMFEM with no orientation fixes; integrated volume matches the VTK volume; Laplace smoke solve (inflow to outlets, levels up to 2 M elements) converges with inflow = outflow |

**CI.** `ci.yml` runs `tests/test_suite.py` on every push: release index and checksums, case
files present, tools compile, and a committed tiny pipe case (`tests/cases/pipe_tiny`, about
15 k tets) through the full verifier and an MFEM round trip. `smoke.yml` (manual only, from the
Actions tab) downloads every case zip from the release (about 3.4 GB), extracts the coarse levels and verifies them, including the PyMFEM
load test and smoke solve. CI does not run svMultiPhysics; `run.sh` is for local runs. Three levels per case support a Grid
Convergence Index study of solver outputs (Celik et al. 2008) once results exist.

## Rebuilding and extending

```sh
pip install vmtk netgen-mesher==6.2.2606 pytetwild vtk numpy
git clone https://github.com/vmtk/SlicerExtension-VMTK ext/SlicerExtension-VMTK
git clone https://github.com/SimVascular/SlicerSimVascular ext/SlicerSimVascular
pip install -e ext/SlicerSimVascular/SimVascularMeshPrep -e ext/SlicerSimVascular/SimVascularMultiPhysics -e ext/SlicerSimVascular/SimVascularROM
export CFDMESHGEN_LIB=ext/SlicerExtension-VMTK/CfdMeshGenerator

python tools/vmr_fetch.py 0031_H_ABAO_AAA            # VMR project -> specs/0031_H_ABAO_AAA.json (faces and BCs from the project)
python tools/build_case.py specs/0031_H_ABAO_AAA.json --res coarse,medium,fine
python tools/vmr_published.py specs/0031_H_ABAO_AAA.json downloads/0031_H_ABAO_AAA/0031_H_ABAO_AAA
python tools/make_synthetic.py && python tools/hex_pipe.py   # synthetic tier
```

`build_case.py` picks the surface edge length for an element target from the (edge length,
count) pairs recorded in the spec, waits for enough free memory, and keeps the attempt closest
to target. Mesher: Netgen (LGPL) by default; fTetWild (MPL) or TetGen (AGPL; needs a VMTK built
with it, as the Slicer extension is) on request.

## Not in this release

ALE/FSI meshes (deferred) and hexahedral meshes beyond the pipe. MFEM ASCII `.mesh` files are
not shipped (see MFEM above).

## Licence, acknowledgement and credits {#credits}

- **Code** (tools, tests, workflows, docs): BSD 3-Clause, see `LICENSE`.
- **Data** keep their sources' terms (`DATA_LICENSE.md`): VMR-derived cases under the VMR
  README-COPYRIGHT (research and development use; the notice travels with every copy and is in
  each VMR case folder and mesh archive); `aortoiliac_svmp_iliac` and `stanford_reference` under
  the svMultiPhysics BSD 3-Clause licence (`svMultiPhysics-LICENSE.txt` in each); synthetic cases
  under the suite's BSD 3-Clause licence.
- **Acknowledgement:** cite the suite and version (`CITATION.cff`); for VMR cases include the VMR
  acknowledgement (NLM R01LM013120 / NHLBI HHSN268201100035C) and/or Wilson et al. 2013,
  doi:10.1115/1.4025983; cite each case's dataset DOI and source study in `REFERENCES.md`. The
  exact wording is in `DATA_LICENSE.md`.
- **Software:** meshing by [SlicerExtension-VMTK](https://github.com/vmtk/SlicerExtension-VMTK)
  (CFD Mesh Generator; VMTK, Antiga et al. 2008) with Netgen (Schöberl 1997); `mesh-complete`,
  `solver.xml` and flow files by [SlicerSimVascular](https://github.com/SimVascular/SlicerSimVascular)
  (SimVascular Mesh Prep, SimVascular MultiPhysics); target solvers
  [svMultiPhysics](https://github.com/SimVascular/svMultiPhysics) and [MFEM](https://mfem.org).
  Full citations, each checked against OpenAlex, are in `REFERENCES.md`.
