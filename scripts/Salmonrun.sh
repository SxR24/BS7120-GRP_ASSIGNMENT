#!/usr/bin/env bash

set -euo pipefail

# Inputs/Outputs
FASTQ_DIR="/media/sf_GroupA_University2026_Project/FASTQ"
SALMON_INDEX="/media/sf_GroupA_University2026_Project/Refrence/salmon_index"
OUT_BASE="/media/sf_GroupA_University2026_Project/Salmon"

WORKDIR="/media/IOtemprun/"
mkdir -p "$WORKDIR"
mkdir -p "$OUT_BASE"

#  Parallelism settings
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-60}
THREADS_PER_SAMPLE=6
MAX_JOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

#  Build a sample list
SAMPLELIST="$WORKDIR/samples.txt"
: > "$SAMPLELIST"

for r1 in "$FASTQ_DIR"/*_R1*.fastq; do
  [ -e "$r1" ] || continue
  r2="${r1/_R1/_R2}"
  sample=$(basename "$r1")
  sample=${sample%%_R1*}
  echo "$sample|$r1|$r2" >> "$SAMPLELIST"
done

echo "Samples found: $(wc -l < "$SAMPLELIST")"

#  Worker function
run_one () {
  local line="$1"
  IFS="|" read -r sample r1 r2 <<< "$line"

  local outdir="$WORKDIR/$sample"
  mkdir -p "$outdir"

  echo "Running Salmon for $sample on $(hostname)"

  salmon quant \
    -i "$SALMON_INDEX" \
    -l A \
    -1 "$r1" \
    -2 "$r2" \
    -p "$THREADS_PER_SAMPLE" \
    --validateMappings \
    --gcBias \
    --seqBias \
    -o "$outdir" \
> "$outdir/salmon.stdout.log" \
    2> "$outdir/salmon.stderr.log"

  mkdir -p "$OUT_BASE/$sample"
  rsync -a "$outdir/" "$OUT_BASE/$sample/"
}

export -f run_one
export SALMON_INDEX OUT_BASE WORKDIR THREADS_PER_SAMPLE

# Parallel execution across samples
xargs -a "$SAMPLELIST" -I{} -P "$MAX_JOBS" bash -c 'run_one "$@"' _ {}

echo "All done. Results in: $OUT_BASE"