# Data licence

The code in this repository is under the BSD 3-Clause licence (`LICENSE`). The case data (geometry,
meshes, boundary conditions) keep the terms of their sources:

| cases | source | terms | notice file shipped with each case level |
|---|---|---|---|
| `0001_H_AO_SVD`, `0011_H_AO_H`, `0031_H_ABAO_AAA`, `0054_H_PULMFON_TAT`, `0065_H_PULMFON_SVD`, `aorta_0159_H_AO_H` | Vascular Model Repository (vascularmodel.com) | VMR README-COPYRIGHT: research and development use; keep the notice with every copy of any portion; acknowledge (below) | `VMR-README-COPYRIGHT.txt` |
| `aortoiliac_svmp_iliac`, `stanford_reference` | svMultiPhysics repository (github.com/SimVascular/svMultiPhysics) | BSD 3-Clause; keep the notice | `svMultiPhysics-LICENSE.txt` |
| `aaa_idealized`, `tcpc_idealized`, `pipe`, `pipe_hex`, `tests/cases/pipe_tiny` | generated here (`tools/make_synthetic.py`, `tools/hex_pipe.py`) | BSD 3-Clause, as the code | `LICENSE` |

Meshes generated here from a source geometry are derived from it and carry its terms.

## Acknowledgement

When you use this suite in a publication or presentation, please:

1. cite the suite: Jolley Lab, *CFD Benchmark Suite*, version used (e.g. v1.0),
   https://github.com/JolleyLab/CFD_FSI_BenchMark_Suite (see `CITATION.cff`);
2. for VMR cases, include the acknowledgement the VMR asks for: "The data used herein was
   provided in whole or in part with Federal funds from the National Library of Medicine under
   Grant No. R01LM013120, and the National Heart, Lung, and Blood Institute, National Institutes
   of Health, Department of Health and Human Services, under Contract No. HHSN268201100035C",
   and/or cite Wilson NM, Ortiz AK, Johnson AB. The Vascular Model Repository: a public resource
   of medical imaging data and blood flow simulation results. *J Med Devices* 2013;7(4):040923,
   doi:10.1115/1.4025983, and/or reference www.vascularmodel.com;
3. cite each case's dataset DOI and source study, listed in `REFERENCES.md`.
