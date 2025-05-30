#!/bin/bash

# Directory containing the input FASTQ files (paired-end)
FASTQ_DIR="/media/raid/Wee/WeeYeZhi/raw_RNA_seq_from_NCBI_SRA"

# Directory to store the output trimmed files
OUTPUT_DIR="/media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results"

# Path to the adapter fasta file used for clipping
ADAPTERS="/media/raid/Wee/WeeYeZhi/raw_RNA_seq_from_NCBI_SRA/adapter.fa"

# Number of threads to use
THREADS=16

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all forward reads (_R1.fastq.gz) in the input directory
for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    # Derive the corresponding reverse read filename
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    # Check if the paired reverse read exists
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Paired file for $R1 not found, skipping."
        continue
    fi

    # Extract the base sample name
    BASENAME=$(basename "$R1" "_1.fastq.gz")

    # Define output filenames for paired and unpaired reads
    OUT_P1="$OUTPUT_DIR/${BASENAME}_1_paired.fastq.gz"
    OUT_U1="$OUTPUT_DIR/${BASENAME}_1_unpaired.fastq.gz"
    OUT_P2="$OUTPUT_DIR/${BASENAME}_2_paired.fastq.gz"
    OUT_U2="$OUTPUT_DIR/${BASENAME}_2_unpaired.fastq.gz"

    # Define trim log file path
    TRIMLOG="$OUTPUT_DIR/${BASENAME}_trimmomatic.log"

    echo "Processing sample $BASENAME ..."

    # Run Trimmomatic PE with specified parameters and trimlog
    trimmomatic PE -threads "$THREADS" \
      -trimlog "$TRIMLOG" \
      "$R1" "$R2" \
      "$OUT_P1" "$OUT_U1" \
      "$OUT_P2" "$OUT_U2" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      SLIDINGWINDOW:4:5 \
      LEADING:5 \
      TRAILING:5 \
      MINLEN:25

    echo "Finished processing $BASENAME"
done

echo "Adapter removal completed for all paired-end .fastq.gz files via trimmomatic."
