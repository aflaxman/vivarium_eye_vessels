#!/bin/bash
# Run a job file of calibration candidates across seeds, a few at a time.
#
# Usage: scripts/sweep/sweep_jobs.sh JOBFILE [SEEDS] [PARALLEL]
#   JOBFILE   lines of NAME<TAB>JSON, where JSON maps dotted spec keys to
#             override values, e.g.
#             v18_si18	{"particles.terminal_velocity": 0.18, "path_splitter.split_interval": 18}
#   SEEDS     comma-separated random seeds (default: the five calibration seeds)
#   PARALLEL  concurrent runs (default 4; one run is one core)
#
# Each NAME x SEED runs as its own process (multiseed_eval.py) so a crash costs
# one cell. Logs go to $SWEEP_DIR/NAME_sSEED.log and results to
# $SWEEP_DIR/NAME_sSEED_multiseed.json; a log ending in MULTI-DONE is complete
# and is skipped on re-runs, so an interrupted sweep can be relaunched as is.
# $SWEEP_DIR defaults to sweep_out/ under the repository root.
set -u
REPO=$(cd "$(dirname "$0")/../.." && pwd)
export REPO
export SWEEP_DIR=${SWEEP_DIR:-$REPO/sweep_out}
mkdir -p "$SWEEP_DIR"
cd "$REPO"
JOBS=$1
SEEDS=${2:-7,42,909,2024,123456}
PARALLEL=${3:-4}
tmp=$(mktemp)
while IFS=$'\t' read -r name json; do
  [ -z "$name" ] && continue
  for seed in ${SEEDS//,/ }; do
    printf '%s\t%s\t%s\n' "$name" "$seed" "$json" >> "$tmp"
  done
done < "$JOBS"
run_one() {
  IFS=$'\t' read -r name seed json <<< "$1"
  log=$SWEEP_DIR/${name}_s${seed}.log
  if grep -q MULTI-DONE "$log" 2>/dev/null; then return 0; fi
  python "$REPO/scripts/sweep/multiseed_eval.py" "$SWEEP_DIR" "${name}_s${seed}" "$json" "$seed" > "$log" 2>&1
  echo "JOB_DONE rc=$?" >> "$log"
}
export -f run_one
xargs -P "$PARALLEL" -d '\n' -I{} bash -c 'run_one "$1"' _ {} < "$tmp"
rm -f "$tmp"
echo "SWEEP_DONE"
