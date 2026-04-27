#!/usr/bin/env bash


# #############################################################################
# STAR Genome Index Generation
# #############################################################################
# Purpose: Build a STAR index from the GRCh38 reference genome and GTF
#          annotation, using an overhang of 99 (for read length 100 bp).
# Input:   GRCh38_index.fa (ensemble), Homo_sapiens.GRCh38.115.gtf (ensemble)
# Output:  STAR index files in D:\GroupA_University2026_Project\Refrence_genes
# #############################################################################

# Generate the STAR index
STAR \
  --runThreadN 30 \
  --runMode genomeGenerate \
  --genomeDir D:\GroupA_University2026_Project\Refrence_genes \
  --genomeFastaFiles D:\GroupA_University2026_Project\Refrence_genes\GRCh38_index.fa \
  --sjdbGTFfile D:\GroupA_University2026_Project\Refrence_genes\Homo_sapiens.GRCh38.115.gtf \
  --sjdbOverhang  99