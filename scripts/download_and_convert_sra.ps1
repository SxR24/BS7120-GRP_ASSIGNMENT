# #############################################################################
# SRA Download and FASTQ Conversion
# #############################################################################
# Purpose: Download SRA files via prefetch from an accession list and convert
#          them to paired‑end FASTQ files using fasterq‑dump.
# Input:   SRR_Acc_List.txt (SRA accession numbers from PRJNA304502)
# Output:  Paired FASTQ files (*_1.fastq, *_2.fastq) in the fastq directory
# #############################################################################

# Download SRA files and prepare output directories
$env:Path = "C:\sratoolkit.3.3.0-win64\bin;" + $env:Path
C:\sratoolkit.3.3.0-win64\bin\prefetch.exe --output-directory D:\SRA --option-file "C:\sratoolkit.3.3.0-win64\bin\SRR_Acc_List.txt"
mkdir "D:\GroupA_University2026_Project\SRA\fastq" -ErrorAction SilentlyContinue
mkdir "D:\GroupA_University2026_Project\SRA\tmp" -ErrorAction SilentlyContinue

# Convert all downloaded .sra files to paired FASTQ
Get-ChildItem "D:\GroupA_University2026_Project\SRA" -Recurse -Filter *.sra | ForEach-Object {
    fasterq-dump --split-files --outdir "D:\GroupA_University2026_Project\SRA\fastq" --temp "D:\GroupA_University2026Project\SRA\tmp" $.FullName
}
