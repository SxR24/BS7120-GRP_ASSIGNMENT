# SRA Download and FASTQ Conversion
#
# Add SRA Toolkit to PATH
$env:Path = "C:\sratoolkit.3.3.0-win64\bin;" + $env:Path

# Download SRA files listed in SRR_Acc_List.txt to D:\SRA
C:\sratoolkit.3.3.0-win64\bin\prefetch.exe --output-directory D:\SRA --option-file "C:\sratoolkit.3.3.0-win64\bin\SRR_Acc_List.txt"

# Create output directories
mkdir "D:\GroupA_University2026_Project\SRA\fastq" -ErrorAction SilentlyContinue
mkdir "D:\GroupA_University2026_Project\SRA\tmp" -ErrorAction SilentlyContinue

# Convert all downloaded .sra files to paired FASTQ files
Get-ChildItem "D:\GroupA_University2026_Project\SRA" -Recurse -Filter *.sra | ForEach-Object {
    fasterq-dump --split-files --outdir "D:\GroupA_University2026_Project\SRA\fastq" --temp "D:\GroupA_University2026Project\SRA\tmp" $.FullName
}
