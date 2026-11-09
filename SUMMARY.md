# Suite summary

Every case level in this release with its element quality (from quality.json). h is the target surface edge length (published VMR meshes: as published, with wall boundary layers, so their aspect ratios are large by design). Coarse levels are for performance and smoke tests.

| case | res | type | elements | points | h (cm) | AR max | min dihedral | max dihedral | min scaled J | SJ<0.2 | path |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0001_H_AO_SVD | coarse | tet4 | 1,045,575 | 188,764 | 0.0637 | 5.04 | 13.4 | 161 | 0.199 | 0.000% | cases/0001_H_AO_SVD/coarse |
| 0001_H_AO_SVD | fine | tet4 | 5,159,149 | 905,423 | 0.04137 | 4.57 | 14.6 | 159 | 0.192 | 0.000% | cases/0001_H_AO_SVD/fine |
| 0001_H_AO_SVD | medium | tet4 | 2,496,985 | 444,688 | 0.05198 | 5.23 | 14 | 161 | 0.145 | 0.000% | cases/0001_H_AO_SVD/medium |
| 0001_H_AO_SVD | published | tet4 | 2,055,236 | 357,362 | - | 411 | 0.506 | 179 | 0.00134 | 33.070% | cases/0001_H_AO_SVD/published |
| 0011_H_AO_H | coarse | tet4 | 1,259,119 | 223,439 | 0.1154 | 5.86 | 10.4 | 164 | 0.0965 | 0.006% | cases/0011_H_AO_H/coarse |
| 0011_H_AO_H | fine | tet4 | 3,993,392 | 705,691 | 0.07993 | 5.84 | 12.1 | 163 | 0.151 | 0.000% | cases/0011_H_AO_H/fine |
| 0011_H_AO_H | medium | tet4 | 2,833,752 | 500,284 | 0.0919 | 5.36 | 13.6 | 162 | 0.165 | 0.000% | cases/0011_H_AO_H/medium |
| 0011_H_AO_H | published | tet4 | 2,989,185 | 554,474 | - | 204 | 0.33 | 179 | 0.00266 | 44.397% | cases/0011_H_AO_H/published |
| 0031_H_ABAO_AAA | coarse | tet4 | 1,000,019 | 184,536 | 0.1274 | 4.98 | 12.1 | 161 | 0.123 | 0.011% | cases/0031_H_ABAO_AAA/coarse |
| 0031_H_ABAO_AAA | fine | tet4 | 4,966,851 | 889,107 | 0.09115 | 5.2 | 12 | 161 | 0.132 | 0.000% | cases/0031_H_ABAO_AAA/fine |
| 0031_H_ABAO_AAA | medium | tet4 | 2,396,588 | 432,758 | 0.1057 | 9.34 | 9.18 | 168 | 0.0696 | 0.004% | cases/0031_H_ABAO_AAA/medium |
| 0031_H_ABAO_AAA | published | tet4 | 2,983,028 | 550,959 | - | 727 | 0.111 | 180 | 0.000287 | 67.759% | cases/0031_H_ABAO_AAA/published |
| 0054_H_PULMFON_TAT | coarse | tet4 | 1,852,725 | 334,919 | 0.05985 | 5.15 | 11.3 | 161 | 0.149 | 0.000% | cases/0054_H_PULMFON_TAT/coarse |
| 0054_H_PULMFON_TAT | fine | tet4 | 4,776,242 | 841,793 | 0.04333 | 5.15 | 12.2 | 161 | 0.159 | 0.000% | cases/0054_H_PULMFON_TAT/fine |
| 0054_H_PULMFON_TAT | medium | tet4 | 3,834,675 | 686,846 | 0.05215 | 5.07 | 11.3 | 161 | 0.129 | 0.000% | cases/0054_H_PULMFON_TAT/medium |
| 0054_H_PULMFON_TAT | published | tet4 | 1,010,672 | 200,758 | - | 830 | 0.185 | 180 | 0.000189 | 51.170% | cases/0054_H_PULMFON_TAT/published |
| 0065_H_PULMFON_SVD | coarse | tet4 | 1,738,199 | 312,376 | 0.08604 | 12 | 6.49 | 166 | 0.0498 | 0.043% | cases/0065_H_PULMFON_SVD/coarse |
| 0065_H_PULMFON_SVD | fine | tet4 | 5,549,704 | 966,246 | 0.05464 | 13 | 3.71 | 169 | 0.0371 | 0.000% | cases/0065_H_PULMFON_SVD/fine |
| 0065_H_PULMFON_SVD | medium | tet4 | 2,509,520 | 451,759 | 0.0714 | 5.79 | 9.16 | 164 | 0.0634 | 0.008% | cases/0065_H_PULMFON_SVD/medium |
| 0065_H_PULMFON_SVD | published | tet4 | 1,060,146 | 206,104 | - | 143 | 0.597 | 179 | 0.00343 | 18.219% | cases/0065_H_PULMFON_SVD/published |
| aaa_idealized | coarse | tet4 | 947,570 | 163,600 | 0.1619 | 3.89 | 17 | 155 | 0.197 | 0.000% | cases/aaa_idealized/coarse |
| aaa_idealized | fine | tet4 | 4,240,401 | 732,000 | 0.09819 | 4.33 | 15.9 | 158 | 0.197 | 0.000% | cases/aaa_idealized/fine |
| aaa_idealized | medium | tet4 | 2,944,823 | 512,719 | 0.1113 | 5.12 | 13.8 | 161 | 0.185 | 0.000% | cases/aaa_idealized/medium |
| aorta_0159_H_AO_H | coarse | tet4 | 959,194 | 170,091 | 0.1654 | 4.83 | 14.9 | 159 | 0.176 | 0.000% | cases/aorta_0159_H_AO_H/coarse |
| aorta_0159_H_AO_H | fine | tet4 | 4,401,147 | 778,827 | 0.1105 | 4.8 | 12.4 | 160 | 0.195 | 0.000% | cases/aorta_0159_H_AO_H/fine |
| aorta_0159_H_AO_H | medium | tet4 | 3,519,334 | 555,089 | 0.1263 | 9.39 | 5.61 | 170 | 0.077 | 0.183% | cases/aorta_0159_H_AO_H/medium |
| aortoiliac_svmp_iliac | coarse | tet4 | 1,162,446 | 206,853 | 0.1412 | 4.42 | 15.6 | 158 | 0.196 | 0.000% | cases/aortoiliac_svmp_iliac/coarse |
| aortoiliac_svmp_iliac | fine | tet4 | 5,161,455 | 908,401 | 0.08345 | 4.92 | 13.2 | 160 | 0.112 | 0.001% | cases/aortoiliac_svmp_iliac/fine |
| aortoiliac_svmp_iliac | medium | tet4 | 2,644,731 | 467,559 | 0.1147 | 5.22 | 12.5 | 161 | 0.084 | 0.003% | cases/aortoiliac_svmp_iliac/medium |
| pipe | coarse | tet4 | 1,311,300 | 224,645 | 0.09142 | 3.57 | 16.9 | 153 | 0.202 | 0.000% | cases/pipe/coarse |
| pipe | fine | tet4 | 5,538,105 | 953,902 | 0.05282 | 4.21 | 15.3 | 157 | 0.228 | 0.000% | cases/pipe/fine |
| pipe | medium | tet4 | 3,854,674 | 667,076 | 0.0657 | 4.42 | 16.7 | 158 | 0.212 | 0.000% | cases/pipe/medium |
| pipe_hex | coarse | hex8 | 1,445,570 | 1,476,792 | 0.0285 | 1.54 | - | - | 0.707 | 0.000% | cases/pipe_hex/coarse |
| pipe_hex | fine | hex8 | 4,672,512 | 4,741,209 | 0.0191 | 1.51 | - | - | 0.707 | 0.000% | cases/pipe_hex/fine |
| pipe_hex | medium | hex8 | 2,828,496 | 2,877,736 | 0.0226 | 1.51 | - | - | 0.707 | 0.000% | cases/pipe_hex/medium |
| stanford_reference_pipe_1M | pipe_1M | tet4 | 1,002,240 | 181,875 | - | 23.4 | 2.3 | 175 | 0.0262 | 2.720% | cases/stanford_reference/pipe_1M |
| tcpc_idealized | coarse | tet4 | 1,100,781 | 192,232 | 0.08865 | 4.8 | 13.2 | 159 | 0.169 | 0.001% | cases/tcpc_idealized/coarse |
| tcpc_idealized | fine | tet4 | 5,749,055 | 924,396 | 0.06613 | 7.49 | 8.44 | 167 | 0.104 | 0.029% | cases/tcpc_idealized/fine |
| tcpc_idealized | medium | tet4 | 3,662,580 | 638,054 | 0.06878 | 4.48 | 13 | 158 | 0.162 | 0.000% | cases/tcpc_idealized/medium |
