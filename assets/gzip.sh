#!/bin/bash

# Define the directory where the fastq files are located
directory="/home/cbr16/Documents/WeeYeZhi/input"

# List of fastq files
files=(
    "SRR9038731.fastq"
    "SRR9038732.fastq"
    "SRR11266556_1.fastq"
    "SRR11266556_2.fastq"
    "SRR11266555_1.fastq"
    "SRR11266555_2.fastq"
    "SRR11266554_1.fastq"
    "SRR11266554_2.fastq"
    "SRR9038729_1.fastq"
    "SRR9038729_2.fastq"
    "SRR9038730_1.fastq"
    "SRR9038730_2.fastq"
    "SRR9038731_1.fastq"
    "SRR9038731_2.fastq"
    "SRR9038732_1.fastq"
    "SRR9038732_2.fastq"
    "SRR9038733_1.fastq"
    "SRR9038733_2.fastq"
    "SRR9038734_1.fastq"
    "SRR9038734_2.fastq"
    "SRR9690969_1.fastq"
    "SRR9690969_2.fastq"
    "SRR9690970_1.fastq"
    "SRR9690970_2.fastq"
    "SRR9690971_1.fastq"
    "SRR9690971_2.fastq"
    "SRR9690972_1.fastq"
    "SRR9690972_2.fastq"
    "SRR9690973_1.fastq"
    "SRR9690973_2.fastq"
    "SRR9690974_1.fastq"
    "SRR9690974_2.fastq"
)

# Change to the specified directory
cd "$directory" || { echo "Directory not found!"; exit 1; }

# Loop through each file and gzip it
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "Compressing $file..."
        gzip "$file"
        echo "$file has been compressed to $file.gz"
    else
        echo "File $file not found!"
    fi
done

echo "All files have been processed."