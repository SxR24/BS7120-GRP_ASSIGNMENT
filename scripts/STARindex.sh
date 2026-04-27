#!/usr/bin/env bash

# Creates STAR index, 99 overhang because 100bp
STAR \
  --runThreadN 30 \
  --runMode genomeGenerate \
  --genomeDir D:\GroupA_University2026_Project\Refrence_genes \
  --genomeFastaFiles D:\GroupA_University2026_Project\Refrence_genes\GRCh38_index.fa \
  --sjdbGTFfile D:\GroupA_University2026_Project\Refrence_genes\Homo_sapiens.GRCh38.115.gtf \
  --sjdbOverhang  99