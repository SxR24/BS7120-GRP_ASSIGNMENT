#!/usr/bin/env bash


# #####################################################################################
# Salmon Index Creation
# #####################################################################################
# Purpose: Build a Salmon index from a gentrome and a
#          decoy list (gentromebuilder.sh outputs) to improve quantification accuracy.
# Input:   gentrome.fa.gz, decoys.txt
# Output:  Salmon index directory at salmon_index
# #####################################################################################

# Build the Salmon index using the gentrome and decoys
salmon index \
  -t /lustre/alice3/scratch/alice/i/ia256/GroupA_University2026_Project/Refrence/gentrome.fa.gz \
  -d decoys.txt \
  -i /lustre/alice3/scratch/alice/i/ia256/GroupA_University2026_Project/Refrence/salmon_index \
  -p 15 \
  -k 31