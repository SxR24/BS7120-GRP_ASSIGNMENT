#!/usr/bin/env bash
set -euo pipefail


# #############################################################################
# TopHat2 Batch Alignment (Paired‑End SRR)
# #############################################################################
# Purpose: Align paired‑end FASTQ files from download_and_convert_sra.ps1 to the
#          GRCh38 genome using TopHat2, with parallel job execution.
# Input:   Paired FASTQ files (*_1.fastq, *_2.fastq) in FASTQ_DIR,
#          GRCh38 Bowtie2 index at INDEX_BASE.
# Output:  Alignment results (BAM, logs) per sample in OUT_BASE.
# #############################################################################
# TopHat2 batch align (paired-end SRR)

#  Config and setup
FASTQ_DIR="/media/sf_GroupA_University2026_Project/SRA/fastq"
INDEX_BASE="/media/sf_GroupA_University2026_Project/Refrence_genes/GRCh38_index"
OUT_BASE="/home/randolando/align"
export WORKDIR="/home/randolando/tmp_tophat"

TOTAL_CPUS="${TOTAL_CPUS:-60}"
THREADS_PER_SAMPLE=15
MAX_JOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

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

  echo "Running $srr on $(hostname)"
  echo "R1=$r1"
  echo "R2=$r2"
  
   #  Set working env to python 2.7
    source $(conda info --base)/etc/profile.d/conda.sh && conda activate tophat2_py2

  #  Run TopHat2 and save results
  tophat2 -p "$THREADS_PER_SAMPLE" -o "$outdir" "$INDEX_BASE" "$r1" "$r2" \
    > "$outdir/tophat.stdout.log" 2> "$outdir/tophat.stderr.log"

  mkdir -p "$OUT_BASE/$srr"
  rsync -a "$outdir/" "$OUT_BASE/$srr/"
}

export -f run_one
export FASTQ_DIR INDEX_BASE OUT_BASE WORKDIR THREADS_PER_SAMPLE

# Parallel execution using xargs
cat "$SAMPLELIST" | xargs -I{} -P "$MAX_JOBS" bash -lc 'run_one "$@"' _ {}

echo "All done. Results in: $OUT_BASE"
