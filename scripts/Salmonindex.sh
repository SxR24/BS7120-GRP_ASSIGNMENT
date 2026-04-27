#!/usr/bin/env bash

# Salmon index creator
salmon index \
  -t /lustre/alice3/scratch/alice/i/ia256/GroupA_University2026_Project/Refrence/gentrome.fa.gz \
  -d decoys.txt \
  -i /lustre/alice3/scratch/alice/i/ia256/GroupA_University2026_Project/Refrence/salmon_index \
  -p 15 \
  -k 31