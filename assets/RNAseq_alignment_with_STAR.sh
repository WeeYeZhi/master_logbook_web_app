#!/bin/bash

# Define paths to the genome index and GTF file
GENOME_DIR="/path/to/genomeDir"
GTF_FILE="/path/to/annotations.gtf"
OUTPUT_DIR="/path/to/output"

# List of SRA IDs
SRA_IDS=("SRR11266554" "SRR11266555" "SRR11266556" "SRR9690969" "SRR9690970" "SRR9690971" "SRR9690972" "SRR9690973" "SRR9690974")

# Run STAR for each SRA ID
for SRA in "${SRA_IDS[@]}"; do
    echo "Processing $SRA..."

    # Create a separate directory for each sample's output
    SAMPLE_DIR="${OUTPUT_DIR}/${SRA}"
    mkdir -p $SAMPLE_DIR

    # Run STAR alignment
    STAR --runThreadN 40 \
         --genomeDir $GENOME_DIR \
         --sjdbGTFfile $GTF_FILE \
         --readFilesIn ${SRA}_1.fastq ${SRA}_2.fastq \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${SAMPLE_DIR}/${SRA}_ \
         --quantMode GeneCounts \
         --twopassMode Basic

    echo "$SRA alignment complete."
done

echo "STAR alignment complete for all samples."
