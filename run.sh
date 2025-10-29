#!/bin/bash
#SBATCH -J neotoma_pollen
#SBATCH -p compregular
#SBATCH --qos=normal
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 30:00:00
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err
#SBATCH --array=1-20

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

START_ID=1
END_ID=100000
CHUNK=5000
THROTTLE=0.3
OUTROOT="neotoma_out_batches"
PY="python3"
SCRIPT="00_neotoma_download_extract_geo.py"
INCLUDE_SPORES=1
MAX_RETRIES=3

mkdir -p "$OUTROOT"

TASK=${SLURM_ARRAY_TASK_ID}
RANGE_START=$(( START_ID + (TASK-1)*CHUNK ))
RANGE_END=$(( RANGE_START + CHUNK - 1 ))
(( RANGE_START > END_ID )) && exit 0
(( RANGE_END > END_ID )) && RANGE_END=$END_ID

OUTDIR="${OUTROOT}/batch_${RANGE_START}_${RANGE_END}"
mkdir -p "$OUTDIR"

ARGS=( --id-start "$RANGE_START" --id-end "$RANGE_END" --outdir "$OUTDIR" --throttle "$THROTTLE" )
(( INCLUDE_SPORES == 0 )) && ARGS+=( --no-spores )

echo "Task ${TASK}: IDs ${RANGE_START}-${RANGE_END}"

n=0
until $PY "$SCRIPT" "${ARGS[@]}"; do
  rc=$?; n=$((n+1))
  (( n >= MAX_RETRIES )) && exit $rc
  sleep $(( 10 * n ))
  echo "Task ${TASK}: retry ${n}..."
done

touch "${OUTDIR}/DONE_${RANGE_START}_${RANGE_END}"
