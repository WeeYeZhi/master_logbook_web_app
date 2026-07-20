#!/bin/bash

# Array of SRA sample names
SAMPLES=(
    "SRR11266554"
    "SRR11266555"
    "SRR11266556"
    "SRR9690969"
    "SRR9690970"
    "SRR9690971"
    "SRR9690972"
    "SRR9690973"
    "SRR9690974"
)

# Number of threads to use for samtools view and sort
THREADS=16

# Loop over each sample
for sample in "${SAMPLES[@]}"
do
    # Define input and output paths
    SAM_FILE="/media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_mapping_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/${sample}/${sample}.sam"
    BAM_FILE="/media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_mapping_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/${sample}/${sample}.bam"
    SORTED_BAM_FILE="/media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_mapping_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/${sample}/${sample}.sorted.bam"
    FLAGSTAT_FILE="/media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_mapping_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/${sample}/${sample}.flagstat"

    echo "Processing sample: $sample"

    # Check if input SAM file exists
    if [[ ! -f "$SAM_FILE" ]]; then
        echo "Warning: SAM file for sample $sample not found at $SAM_FILE. Skipping."
        continue
    fi

    # Convert SAM to BAM (compress)
    samtools view -@ $THREADS -Sb "$SAM_FILE" > "$BAM_FILE"
    if [[ $? -ne 0 ]]; then
        echo "Error: samtools view failed for sample $sample."
        continue
    fi

    # Sort the BAM file
    samtools sort -@ $THREADS -o "$SORTED_BAM_FILE" "$BAM_FILE"
    if [[ $? -ne 0 ]]; then
        echo "Error: samtools sort failed for sample $sample."
        continue
    fi

    # Index the sorted BAM file
    samtools index "$SORTED_BAM_FILE"
    if [[ $? -ne 0 ]]; then
        echo "Error: samtools index failed for sample $sample."
        continue
    fi

    # Generate mapping statistics using the sorted BAM file
    samtools flagstat "$SORTED_BAM_FILE" > "$FLAGSTAT_FILE"
    if [[ $? -ne 0 ]]; then
        echo "Error: samtools flagstat failed for sample $sample."
        continue
    fi

    echo "Finished processing sample: $sample"
done

echo "All samples processed."