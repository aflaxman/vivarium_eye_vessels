#!/bin/bash
# Submit a job file of calibration candidates to Slurm, one array task per NAME x SEED.
#
# Usage: scripts/sweep/sweep_slurm.sh JOBFILE [SEEDS] [ACCOUNT] [PARTITION]
#   Same JOBFILE format as sweep_jobs.sh (NAME<TAB>JSON). Each task runs
#   multiseed_eval.py for one cell inside the project's conda env; logs and
#   results land in $SWEEP_DIR exactly as sweep_jobs.sh writes them, so
#   seed_matrix.py / sweep_table.py read both. Completed cells (MULTI-DONE in
#   the log) are skipped, so a partial sweep can be resubmitted as is.
set -eu
REPO=$(cd "$(dirname "$0")/../.." && pwd)
export SWEEP_DIR=${SWEEP_DIR:-$REPO/sweep_out}
mkdir -p "$SWEEP_DIR"
JOBS=$1
SEEDS=${2:-7,42,909,2024,123456}
ACCOUNT=${3:-proj_simscience}
PARTITION=${4:-all.q}
ENV=${CONDA_ENV:-vivarium_eye_vessels}
cells=$SWEEP_DIR/cells_$(date +%s).tsv
: > "$cells"
while IFS=$'\t' read -r name json; do
  [ -z "$name" ] && continue
  for seed in ${SEEDS//,/ }; do
    grep -q MULTI-DONE "$SWEEP_DIR/${name}_s${seed}.log" 2>/dev/null && continue
    printf '%s\t%s\t%s\n' "$name" "$seed" "$json" >> "$cells"
  done
done < "$JOBS"
n=$(wc -l < "$cells")
if [ "$n" -eq 0 ]; then echo "nothing to run"; exit 0; fi
sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -c 1 --mem=6G -t 4:00:00 \
  -J vev_sweep --array="1-$n" -o "$SWEEP_DIR/slurm_%A_%a.out" \
  --wrap "cd $REPO && line=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" $cells) && IFS=\$'\t' read -r name seed json <<< \"\$line\" && conda run --no-capture-output -n $ENV python scripts/sweep/multiseed_eval.py $SWEEP_DIR \${name}_s\${seed} \"\$json\" \$seed > $SWEEP_DIR/\${name}_s\${seed}.log 2>&1; echo JOB_DONE rc=\$? >> $SWEEP_DIR/\${name}_s\${seed}.log"
echo "submitted $n cells from $cells"
