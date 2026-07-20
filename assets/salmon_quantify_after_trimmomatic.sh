#!/bin/bash

echo "Starting Salmon quantification pipeline..."

echo "Running Salmon to quantify the transcript abundance of SRR11266554..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR11266554_after_trimmomatic
echo "Completed quantification for SRR11266554."

echo "Running Salmon to quantify the transcript abundance of SRR11266555..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR11266555_after_trimmomatic
echo "Completed quantification for SRR11266555."

echo "Running Salmon to quantify the transcript abundance of SRR11266556..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR11266556_after_trimmomatic
echo "Completed quantification for SRR11266556."

echo "Running Salmon to quantify the transcript abundance of SRR9690969..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR9690969_after_trimmomatic
echo "Completed quantification for SRR9690969."

echo "Running Salmon to quantify the transcript abundance of SRR9690970..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR9690970_after_trimmomatic
echo "Completed quantification for SRR9690970."

echo "Running Salmon to quantify the transcript abundance of SRR9690971..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR9690971_after_trimmomatic
echo "Completed quantification for SRR9690971."

echo "Running Salmon to quantify the transcript abundance of SRR9690972..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR9690972_after_trimmomatic
echo "Completed quantification for SRR9690972."

echo "Running Salmon to quantify the transcript abundance of SRR9690973..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR9690973_after_trimmomatic
echo "Completed quantification for SRR9690973."

echo "Running Salmon to quantify the transcript abundance of SRR9690974..."
salmon quant -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic \
-l ISR \
-1 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz \
-2 /media/raid/Wee/WeeYeZhi/output/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz \
-p 16 \
--validateMappings \
-o /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_quantify/salmon_quantify_after_trimmomatic/quant_SRR9690974_after_trimmomatic
echo "Completed quantification for SRR9690974."

echo "All Salmon transcripts abundance estimation runs for all the reads have been completed and finished."