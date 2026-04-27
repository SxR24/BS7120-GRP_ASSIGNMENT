#  Inputs/Outputs
FASTQ_DIR="/media/sf_GroupA_University2026_Project/SRA/fastq2"
OUT_DIR="/media/sf_GroupA_University2026_Project/rRNA_results/rRNA_results"

mkdir -p "$OUT_DIR"

#  Uses bowtie2 to check rRNA levels
for fq1 in "$FASTQ_DIR"/*_1.fastq; do
    sample=${fq1%_1.fastq}
    fq2="${sample}_2.fastq"
    base=$(basename "$sample")

    bowtie2 -x /home/i/ia256/rRNAindex/rRNA_index \
        -1 "$fq1" -2 "$fq2" \
        -S "$OUT_DIR/${base}_rRNA.sam" \
        -p 15 --very-sensitive

    mapped=$(samtools view -c -F 4 "$OUT_DIR/${base}_rRNA.sam")
    total=$(samtools view -c "$OUT_DIR/${base}_rRNA.sam")
    percent=$(echo "scale=2; $mapped/$total*100" | bc)

    echo "$base: $mapped / $total reads mapped to rRNA ($percent%)" >> "$OUT_DIR/rRNA_summary.txt"
done