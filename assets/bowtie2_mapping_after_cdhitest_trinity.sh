!/bin/bash

# Define the sample names array
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

# Define the bowtie2 index prefix (basename of the bowtie index files)
INDEX_PREFIX="/media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_build_results/trinity_assembly_after_cdhitest_trinity/trinity_assembly_after_cdhitest_trinity"

# Number of threads to use per bowtie2 run
THREADS=16

# Loop through each sample and run bowtie2
for SAMPLE in "${SAMPLES[@]}"; do
    # Input files: trimmed paired-end reads
    READ1="/media/cbr14/Two/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/${SAMPLE}_1_paired.fastq.gz"
    READ2="/media/cbr14/Two/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/${SAMPLE}_2_paired.fastq.gz"

    # Output directory and SAM output file
    OUTDIR="/media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_mapping_results/trinity_assembly_after_cdhitest_trinity/${SAMPLE}"
    OUTSAM="${OUTDIR}/${SAMPLE}.sam"

    # Create output directory if it does not exist
    mkdir -p "${OUTDIR}"

    echo "Running bowtie2 for sample ${SAMPLE}..."
    bowtie2 -p ${THREADS} -q --no-unal -x "${INDEX_PREFIX}" -1 "${READ1}" -2 "${READ2}" -S "${OUTSAM}"

    if [ $? -eq 0 ]; then
        echo "bowtie2 mapping for ${SAMPLE} completed successfully."
    else
        echo "Error running bowtie2 for ${SAMPLE}" >&2
        exit 1
    fi
done

echo "All samples have been processed."