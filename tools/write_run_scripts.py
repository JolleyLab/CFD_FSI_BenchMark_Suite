#!/usr/bin/env python3
"""Write run.sh into every case folder so a case runs with one command, no setup:

    cd cases/0031_H_ABAO_AAA/medium && ./run.sh 8 20     # 8 MPI ranks, 20 time steps

Order of preference for the solver: $SVMULTIPHYSICS, `svmultiphysics` on PATH, then the
published Docker image simvascular/solver:latest (svMultiPhysics built with Trilinos and PETSc;
see Docker/README.md in SimVascular/svMultiPhysics). The step count rewrites
Number_of_time_steps in a copy of solver.xml (solver.run.xml), so the shipped file stays as is.
Results land in <ranks>-procs/ as svMultiPhysics writes them; wall time goes to run_<ranks>x<steps>.log.
"""
import glob, os, stat

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

SCRIPT = r'''#!/usr/bin/env bash
# Run this cvbench case with svMultiPhysics.  Usage: ./run.sh [ranks=1] [steps=as in solver.xml]
# Solver: $SVMULTIPHYSICS, else `svmultiphysics` on PATH, else Docker image simvascular/solver:latest.
set -euo pipefail
cd "$(dirname "$0")"
ranks=${1:-1}
steps=${2:-}
xml=solver.xml
if [[ -n "$steps" ]]; then
  sed -E "s#<Number_of_time_steps>[^<]*</Number_of_time_steps>#<Number_of_time_steps>${steps}</Number_of_time_steps>#" solver.xml > solver.run.xml
  xml=solver.run.xml
fi
log="run_${ranks}x${steps:-full}.log"
exe=${SVMULTIPHYSICS:-$(command -v svmultiphysics || true)}
start=$(date +%s.%N)
if [[ -n "$exe" ]]; then
  if [[ "$ranks" -gt 1 ]]; then mpiexec -np "$ranks" "$exe" "$xml" | tee "$log"; else "$exe" "$xml" | tee "$log"; fi
else
  command -v docker >/dev/null || { echo "Need svmultiphysics on PATH, \$SVMULTIPHYSICS, or Docker." >&2; exit 1; }
  img=${SVMP_IMAGE:-simvascular/solver:latest}
  bin=${SVMP_IMAGE_BIN:-/build-trilinos/svMultiPhysics-build/bin/svmultiphysics}
  docker run --rm --platform linux/amd64 -v "$PWD":/case -w /case "$img" \
    bash -lc "mpiexec --allow-run-as-root -np $ranks $bin $xml" | tee "$log"
fi
end=$(date +%s.%N)
echo "wall_time_s $(awk "BEGIN{printf \"%.3f\", $end - $start}")  ranks $ranks  steps ${steps:-full}" | tee -a "$log"
'''


def main(roots=None):
    """Write run.sh beside every solver.xml under the given folders (default: cases/)."""
    n = 0
    roots = roots or [os.path.join(ROOT, "cases")]
    for xml in [x for r in roots for x in glob.glob(os.path.join(r, "**", "solver.xml"), recursive=True)]:
        path = os.path.join(os.path.dirname(xml), "run.sh")
        open(path, "w").write(SCRIPT)
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        n += 1
    print(f"wrote run.sh in {n} case folders")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
