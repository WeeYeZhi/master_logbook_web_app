#!/bin/bash

# Directory containing the FASTQ files
FASTQ_DIR="/home/cbr16/Documents/WeeYeZhi/output/fastpresults"
# Directory to store the output files
OUTPUT_DIR="/home/cbr16/Documents/WeeYeZhi/output/processed/falcoresults"

# Loop through all .fastq.gz files in the current directory
for FILE in "$FASTQ_DIR"/*.fastq.gz; do
    # Check if there are any .fastq.gz files
    if [[ -e "$FILE" ]]; then
        # Extract the base name without extension
        BASENAME=$(basename "$FILE" .fastq.gz)
        
        # Process the file with falco
        falco "$FILE" -o "$OUTPUT_DIR" \
            -D "$OUTPUT_DIR/${BASENAME}_fastqc_data.txt" \
            -R "$OUTPUT_DIR/${BASENAME}_fastqc_report.html" \
            -S "$OUTPUT_DIR/${BASENAME}_summary.txt"
    else
        echo "No .fastq.gz files found."
        break
    fi
done

echo "Quality check completed for all .fastq.gz files."