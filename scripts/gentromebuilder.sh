#!/usr/bin/env bash

# concatenating the transcriptome and the genome (Gentrome)

grep "^>" <(gunzip -c /home/i/ia256/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat /home/i/ia256/Downloads/gencode.v49.transcripts.fa.gz /home/i/ia256/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > gentrome.fa.gz