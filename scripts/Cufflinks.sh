#!/usr/bin/env bash
set -euo pipefail

# Cufflinks batch processing

# Input/Output
ALIGN_DIR="/home/randolando/tmp_tophat/"           
GTF_FILE="/media/sf_GroupA_University2026_Project/Refrence_genes/gencode.v22.primary_assembly.annotation.gtf"
WORKDIR="${WORKDIR:-/home/randolando/tmp_cufflinks}"

# CPU/Threads allocation
TOTAL_CPUS="${TOTAL_CPUS:-6}"
MAX_JOBS=$(( TOTAL_CPUS / 3 ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

mkdir -p "$WORKDIR" "/home/randolando/cufflinks_out"

SAMPLELIST="$WORKDIR/samples.txt"
: > "$SAMPLELIST"

# Find all TopHat2 output directories with accepted_hits.bam
for sample_dir in "$ALIGN_DIR"/SRR*; do
  [ -d "$sample_dir" ] || continue
  
  srr=$(basename "$sample_dir")
  bam="$sample_dir/accepted_hits.bam"
  
  if [ ! -f "$bam" ]; then
    echo "WARNING: No accepted_hits.bam found for $srr, skipping" >&2
    continue
  fi
  
  echo "$srr|$bam" >> "$SAMPLELIST"
done

num_samples=$(wc -l < "$SAMPLELIST")
echo "Samples found: $num_samples"

if [ "$num_samples" -eq 0 ]; then
  echo "ERROR: No samples with accepted_hits.bam found" >&2
  exit 1
fi

# worker function
run_cufflinks () {
  local line="$1"
  IFS='|' read -r srr bam <<< "$line"
  
  local outdir="$WORKDIR/$srr"
  mkdir -p "$outdir"
  
  echo "Processing $srr on $(hostname) at $(date)"
  echo "BAM: $bam"
  
   #  Set working env to python 3.6
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate py36

  # Run Cufflinks
  cufflinks \
    -p 3 \
    -o "$outdir" \
    -G "$GTF_FILE" \
    $MULTI_READ_CORRECT \
    $FRAG_BIAS_CORRECT \
    "$bam" \
    > "$outdir/cufflinks.stdout.log" 2> "$outdir/cufflinks.stderr.log"
  
  local exit_code=$?
  
  if [ $exit_code -ne 0 ]; then
    echo "ERROR: Cufflinks failed for $srr (exit code: $exit_code)" >&2
    return $exit_code
  fi
  
  # Copy results to persistent storage
  mkdir -p "/home/randolando/cufflinks_out/$srr"
  rsync -a "$outdir/" "/home/randolando/cufflinks_out/$srr/"
  
  echo "=== Completed $srr at $(date) ==="
  return 0
}

export -f run_cufflinks
export WORKDIR GTF_FILE MULTI_READ_CORRECT FRAG_BIAS_CORRECT

# parallel execution
# Uses xargs parallelism

cat "$SAMPLELIST" | xargs -I{} -P "$MAX_JOBS" bash -c 'run_cufflinks "{}"'

# Completion

echo "All samples processed"
echo "Results directory: $OUT_BASE"