#!/usr/bin/env bash
set -euo pipefail

# TopHat2 batch align (paired-end SRR)

#  Inputs / Outputs
FASTQ_DIR="/media/sf_GroupA_University2026_Project/SRA/fastq"
INDEX_BASE="/media/sf_GroupA_University2026_Project/Refrence_genes/GRCh38_index"
OUT_BASE="/home/randolando/align"
export WORKDIR="/home/randolando/tmp_tophat"

# Default to 32 threads (can be overridden for target machine)
TOTAL_CPUS="${TOTAL_CPUS:-30}"

#  Threads per TopHat2 process
THREADS_PER_SAMPLE=6

#  How many samples to run concurrently
MAX_JOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

echo "Using TOTAL_CPUS=$TOTAL_CPUS"
echo "THREADS_PER_SAMPLE=$THREADS_PER_SAMPLE"
echo "MAX_JOBS=$MAX_JOBS"
echo "WORKDIR=$WORKDIR"

#  Build a sample list
SAMPLELIST="$WORKDIR/samples.txt"
: > "$SAMPLELIST"

for r1 in "$FASTQ_DIR"/SRR*_1.fastq*; do
  base=$(basename "$r1")
  srr=${base%%_1.fastq*}
  r2="$FASTQ_DIR/${srr}_2.fastq${base#*_1.fastq}"

  if [ ! -f "$r2" ]; then
    echo "Missing pair for $srr: expected $r2" >&2
    continue
  fi

  echo "$srr|$r1|$r2" >> "$SAMPLELIST"
done

echo "Samples found: $(wc -l < "$SAMPLELIST")"

#  Worker function
run_one () {
  local line="$1"
  local srr r1 r2
  srr="${line%%|*}"
  line="${line#*|}"
  r1="${line%%|*}"
  r2="${line#*|}"

  local outdir="$WORKDIR/$srr"
  mkdir -p "$outdir"

  echo "=== Running $srr on $(hostname) ==="
  echo "R1=$r1"
  echo "R2=$r2"
  
   #  Set working env to python 2.7
    source $(conda info --base)/etc/profile.d/conda.sh && conda activate tophat2_py2

  #  Run TopHat2
  tophat2 -p "$THREADS_PER_SAMPLE" -o "$outdir" "$INDEX_BASE" "$r1" "$r2" \
    > "$outdir/tophat.stdout.log" 2> "$outdir/tophat.stderr.log"

  # Copy results back to persistent storage
  mkdir -p "$OUT_BASE/$srr"
  rsync -a "$outdir/" "$OUT_BASE/$srr/"
}

export -f run_one
export FASTQ_DIR INDEX_BASE OUT_BASE WORKDIR THREADS_PER_SAMPLE

# Parallel execution across samples
# Uses xargs parallelism (Potentially change)
cat "$SAMPLELIST" | xargs -I{} -P "$MAX_JOBS" bash -lc 'run_one "$@"' _ {}

echo "All done. Results in: $OUT_BASE"
