# Catalogue: open clinical geometries (VMR)

Source: `dataset/dataset-svprojects.csv` and `dataset-svresults.csv` in
[SimVascular/vascularmodel](https://github.com/SimVascular/vascularmodel) (the VMR website
repository; 316 projects). Every VMR project zip (`https://www.vascularmodel.com/svprojects/<NAME>.zip`)
is a SimVascular project: `Images/` (.vti), `Paths/`, `Segmentations/`, `Models/<model>.vtp +
.mdl` (closed surface with `ModelFaceID`, face names and cap/wall types), `Meshes/` (.vtu/.vtp
as published), `flow-files/inflow_3d.flow` (svSolver sign: negative into the model), and
`Simulations/<job>.sjb` (+ `rcrt.dat`, `*.svpre`) holding the tuned RCR / resistance values.
Projects marked "results" also have time-resolved 3D results (`svresults/<NAME>/<NAME>_3D_RIGID_VTU.zip`,
several GB) that can serve as reference solutions.

Licence (all VMR data, `README-COPYRIGHT` in each zip and in `surfaces/VMR-README-COPYRIGHT.txt`):
copyright Stanford / UC Regents / OSMSC; free use for research and development; the notice
must accompany any copy; acknowledge NLM R01LM013120 and NHLBI HHSN268201100035C and/or cite
Wilson, Ortiz, Johnson, J. Med. Devices 2013 (doi:10.1115/1.4025983); do not attempt
re-identification. This permits redistribution of derived meshes in an open suite with the
notice attached; it is not an OSI/CC licence, so say "VMR licence" rather than "CC-BY".

## Recommended 5-case clinical set

| # | role | VMR model | subject | why |
|---|---|---|---|---|
| 1 | healthy adult thoracic aorta | **0011_H_AO_H** (0159_H_AO_H also included) | F 23, MR, results available | arch + 3 branches, pulsatile RCR; the classic aorta CFD case (LaDisa 2011 control group) |
| 2 | pediatric aorta, single-ventricle | **0001_H_AO_SVD** | M 3 y, MR, results available | CHD aorta from the Fontan cohort (Marsden 2010), small vessels + high heart rate: a different time-step/mesh regime |
| 3 | abdominal aortic aneurysm | **0031_H_ABAO_AAA** | M 73, CT, results available | the AAA used in the svMultiPhysics paper and in Sujal's GPU table (1.58 M elements CFD, 571 k FSI); CFD and ALE-FSI (thick wall) |
| 4 | Fontan (TCPC), young | **0054_H_PULMFON_TAT** | M 3 y, MR, results available | Marsden 2010 patient C; two caval inflows, pulmonary resistance outlets, low Reynolds, 100-step runs in Sujal's table (475 k elements) |
| 5 | Fontan (TCPC), older | **0065_H_PULMFON_SVD** (alternative: 0163_H_PULMFON_SVD, TEVG, Schwarz 2021) | M 9 y, MR, results available | second Fontan anatomy for variability; the TEVG models have 6-month / 3-year pairs of the same patients |

Also included: `aortoiliac_svmp_iliac` (svMultiPhysics `tests/cases/fluid/iliac_artery`, the
shared iliac performance case) remeshed at the three sizes.

Candidates with their catalogue data:

| VMR model | legacy | sex / age | anatomy / disease / procedure | modality | mesh / sim / results | zip MB | DOI / paper |
|---|---|---|---|---|---|---:|---|
| 0011_H_AO_H | 0094_0001 | F / 23 | Aorta / Healthy | MR | yes / yes / 3D rigid VTP+VTU | 224 | 10.1111/j.1747-0803.2011.00553.x |
| 0007_H_AO_H | 0090_0001 | M / 13 | Aorta / Healthy | MR | yes / yes / 3D rigid | 240 | 10.1111/j.1747-0803.2011.00553.x |
| 0012_H_AO_H | 0095_0001 | M / 26 | Aorta / Healthy | MR | yes / yes / 3D rigid | 251 | 10.1111/j.1747-0803.2011.00553.x |
| 0159_H_AO_H | O690801_2007_aorta | – | Aorta / Healthy (aorta of coronary-SVG CT) | CT | model only (included) | 2 | 10.1007/s10439-024-03611-z |
| 0001_H_AO_SVD | 0063_1001 | M / 3 | Aorta / Single ventricle (Fontan pt, no arch reconstruction) | MR | yes / yes / 3D rigid | 150 | 10.1111/j.1747-0803.2010.00383.x |
| 0002_H_AO_SVD | 0064_1001 | F / 6 | Neo-aorta / Single ventricle / arch reconstruction | MR | yes / yes / 3D rigid | 144 | 10.1111/j.1747-0803.2010.00383.x |
| 0016_H_AO_COA | 0104_0001 | F / 11 | Aorta / Coarctation / end-to-end repair | MR | yes / yes / 3D rigid | 230 | 10.1111/j.1747-0803.2011.00553.x |
| 0031_H_ABAO_AAA | 0138_1001 | M / 73 | Abdominal aorta / AAA | CT | yes / yes / 3D rigid | 274 | 10.25740/nr862gy9188 |
| 0032_H_ABAO_AAA | 0139_1001 | F / 79 | Abdominal aorta / AAA | CT | yes / yes / 3D rigid | 257 | 10.25740/wy327qj3901 |
| 0033_H_ABAO_AAA | 0140_2001 | M / 76 | Abdominal aorta / AAA | CT | yes / yes / 3D rigid | 341 | 10.25740/ff987bh8208 |
| 0054_H_PULMFON_TAT | 0063_0001 | M / 3 | Pulmonary Fontan / tricuspid atresia | MR | yes / yes / 3D rigid | 90 | 10.1111/j.1747-0803.2010.00383.x |
| 0055_H_PULMFON_HLHS | 0064_0001 | F / 6 | Pulmonary Fontan / HLHS | MR | yes / yes / 3D rigid | 77 | 10.1111/j.1747-0803.2010.00383.x |
| 0057_H_PULMFON_TAT | 0075_0001 | F / 17 | Pulmonary Fontan / tricuspid atresia | MR | yes / yes / 3D rigid | 110 | 10.1111/j.1747-0803.2010.00383.x |
| 0064_H_PULMFON_SVD | 0125_0001 | M / 5 | Pulmonary Fontan | MR | yes / yes / 3D rigid | 90 | 10.25740/xp014wf4919 |
| 0065_H_PULMFON_SVD | 0126_0001 | M / 9 | Pulmonary Fontan | MR | yes / yes / 3D rigid | 96 | 10.25740/dx598vp4827 |
| 0161–0168_H_PULMFON_SVD | PT1–4, 6M/3Y | 4 pts, 2.5–6 y | Pulmonary Fontan / TEVG conduit | MR | yes / yes / – | 400–950 | 10.1038/s41536-021-00148-w |
| 0060–0063_H_PULMGLN_SVD | 0096–0099 | 3–5 y | Pulmonary Glenn | MR | yes / yes / 3D rigid | ~150 | 10.1115/1.4005377 |
| 0047_H_ABAO_AIOD | 0161_0001 | F / 67 | Abdominal aorta / aortoiliac occlusive disease, bypass | MR | yes / yes / 3D rigid | 265 | 10.3109/10929080500230445 |

The VMR also holds 15 more AAA projects contributed in 2023 (0179–0198_H_ABAO_AAA, images from
Ken Tran, models but mostly no meshes/simulations) and 87 coronary projects, not needed here.

## Other open sources checked

- **SimVascular/svMultiPhysics `tests/cases`**: iliac_artery (used here), pipes, cylinder, LV,
  shells. Real geometry but small; good CI tests, not production-size.
- **SimVascular/svFSI-Tests**: pipe FSI (ALE and CMM), vein-graft CMM (`07-fsi/cmm/02-vein-graft`,
  patient-specific stenosed vein graft, ~44 MB with meshes), 3D-0D coupled patient-specific
  pulmonary case (`04-fluid/04-3D0D-coupling-BC/patient-specific-closed-loop-case`, 19 outlets).
  Both are candidates for a CMM/closed-loop extension of the suite.
- **vmtk/vmtk-test-data**: `input/vmr-0159-aorta-surface.vtp` (VMR aorta, used here) and the
  VMTK tutorial aorta (`aorta-surface.vtp`, mm, abdominal segment); tutorial data of uncertain
  provenance, not used.
- **SimTK sv_tests** (Aortofemoral Normal 1/2, Coronary Normal, Healthy Pulmonary, Demo
  Project = VMR 0110 aorta) and **Stanford Digital Repository** mirrors of every VMR model
  (`purl.stanford.edu/<druid>`, listed in the CSV `SDR DOI` column): same data and licence as
  the VMR zips.

## Status

The recommended five (0031, 0054, 0065, 0001, 0011) are built: published mesh plus coarse,
medium and fine levels each (see README). Any other model above is added with
`tools/vmr_fetch.py <ID>` then `tools/build_case.py specs/<ID>.json`.
