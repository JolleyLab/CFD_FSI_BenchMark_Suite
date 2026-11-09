#!/usr/bin/env bash
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
