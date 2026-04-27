#!/usr/bin/env bash


# #############################################################################
# STAR Alignment Pipeline with Shared Genome Memory
# #############################################################################
# Purpose: Align paired‑end FASTQ reads to the GRCh38 genome using STAR in
#          parallel, leveraging shared memory genome loading for efficiency.
#          Quantification via GeneCounts and TranscriptomeSAM.
# Input:   Paired FASTQ files (SRR*_1.fastq, SRR*_2.fastq) 
#          in FASTQ_DIR (from download_and_convert_sra.ps1),
#          STAR index at STAR_INDEX (from STARindex.sh).
# Output:  Sorted BAMs, counts, and logs per sample in OUT_BASE.
# #############################################################################

#  Config and setup
FASTQ_DIR="/media/sf_GroupA_University2026_Project/SRA/fastq"
STAR_INDEX="/media/sf_GroupA_University2026_Project/Refrence_genes_Alt"
OUT_BASE="/media/sf_GroupA_University2026_Project/STARsalmon"

WORKDIR="/media/IOtemprun/" #  SSD working dir (speed)
mkdir -p "$WORKDIR"
mkdir -p "$OUT_BASE"

TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-60}
THREADS_PER_SAMPLE=12
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

# Pre-load genome into shared memory (comment to disable)
STAR --genomeDir "$STAR_INDEX" --genomeLoad LoadAndExit

# Parallel execution using xargs
cat "$SAMPLELIST" | xargs -I{} -P "$MAX_JOBS" bash -c 'run_one "$@"' _ {}

# Cleanup shared memory genome
STAR --genomeDir "$STAR_INDEX" --genomeLoad Remove

echo "All done. Results in: $OUT_BASE"