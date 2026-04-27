#!/usr/bin/env bash


# #############################################################################################
# Gentrome Generation (Transcriptome + Genome Concatenation)
# #############################################################################################
# Purpose: Create a decoy‑aware "gentrome" file for transcriptome indexing by
#          extracting genome sequence IDs (decoys) and concatenating the
#          transcriptome and genome FASTA files.
# Input:   gencode.v49.transcripts.fa.gz (transcriptome from gencode)
#          Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (latest genome from ensemble)
# Output:  decoys.txt (list of genome sequence IDs)
#          gentrome.fa.gz (concatenated transcriptome + genome)
# #############################################################################################

# Extract genome sequence IDs for the decoy list
grep "^>" <(gunzip -c /home/i/ia256/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# Concatenate transcriptome and genome into one file
cat /home/i/ia256/Downloads/gencode.v49.transcripts.fa.gz /home/i/ia256/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > gentrome.fa.gz