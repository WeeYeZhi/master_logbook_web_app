#!/bin/bash

# Directory containing the FASTQ files
FASTQ_DIR="/media/raid/Wee/WeeYeZhi/raw_RNA_seq_from_NCBI_SRA"
# Directory to store the output files
OUTPUT_DIR="/media/raid/Wee/WeeYeZhi/output/fastqc_raw_RNA_seq_from_NCBI_SRA_results"

# Loop through all .fastq.gz files in the current directory
for FILE in "$FASTQ_DIR"/*.fastq.gz; do
    # Check if there are any .fastq.gz files
    if [[ -e "$FILE" ]]; then
        # Extract the base name without extension
        BASENAME=$(basename "$FILE" .fastq.gz)
        
        # Process the file with falco
        fastqc -f fastq -t 16 "$FILE" -o "$OUTPUT_DIR"
    else
        echo "No .fastq.gz files found."
        break
    fi
done

echo "Quality check completed for all the .fastq.gz files of raw CPB transcriptomic reads (derived from NCBI SRA)."