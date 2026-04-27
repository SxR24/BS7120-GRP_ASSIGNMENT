#!/usr/bin/env bash

#  Inputs/Outputs
FASTQ_DIR="/media/sf_GroupA_University2026_Project/SRA/fastq"
STAR_INDEX="/media/sf_GroupA_University2026_Project/Refrence_genes_Alt"
OUT_BASE="/media/sf_GroupA_University2026_Project/STARsalmon"

#  SSD working dir (fast)
WORKDIR="/media/IOtemprun/"
mkdir -p "$WORKDIR"
mkdir -p "$OUT_BASE"

#  Parallelism settings
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-60}
THREADS_PER_SAMPLE=12
MAX_JOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

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

  echo "Running STAR for $srr on $(hostname)"
  echo "R1=$r1"
  echo "R2=$r2"

  STAR \
    --runThreadN "$THREADS_PER_SAMPLE" \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$r1" "$r2" \
    --outFileNamePrefix "$outdir/" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts TranscriptomeSAM \
    --genomeLoad LoadAndKeep \
    > "$outdir/star.stdout.log" 2> "$outdir/star.stderr.log"

  mkdir -p "$OUT_BASE/$srr"
  rsync -a "$outdir/" "$OUT_BASE/$srr/"
}

export -f run_one
export FASTQ_DIR STAR_INDEX OUT_BASE WORKDIR THREADS_PER_SAMPLE

# Pre-load genome into shared memory
STAR --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit

# Parallel execution across samples
cat "$SAMPLELIST" | xargs -I{} -P "$MAX_JOBS" bash -c 'run_one "$@"' _ {}

# Cleanup shared memory genome
STAR --genomeDir "$STAR_INDEX" --genomeLoad Remove

echo "All done. Results in: $OUT_BASE"