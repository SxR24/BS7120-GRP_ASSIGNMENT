#!/usr/bin/env bash

# ########################################################################################################
# rRNA Index Creation for Bowtie2
# ########################################################################################################
# Purpose: Generate a Bowtie2 index from a FASTA file of human rRNA sequences.
# Input:   4V6X_human_rRNAs.fa sourced from the ARF R package)
# Output:  A set of index files with the prefix 'rRNA_index' for use in rRNAcheck.sh
# ########################################################################################################

# Build the Bowtie2 index
bowtie2-build /home/i/ia256/Downloads/4V6X_human_rRNAs.fa rRNA_index