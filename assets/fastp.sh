#!/bin/bash

# Directory containing the FASTQ files
FASTQ_DIR="/home/cbr16/Documents/WeeYeZhi/input"
# Directory to store the output files
OUTPUT_DIR="/home/cbr16/Documents/WeeYeZhi/input/../output/fastpresults"

# Loop through each file in the FASTQ directory
for file in "$FASTQ_DIR"/*_1.fastq.gz; do
	# Get the base name of the file (without the _R1.fastq part)
	base_name=$(basename "$file" "_1.fastq.gz")

	# Define input files for read1 and read2
        input_R1="${FASTQ_DIR}/${base_name}_1.fastq.gz"
        input_R2="${FASTQ_DIR}/${base_name}_2.fastq.gz"

        # Define output files for read1 and read2
        output_R1="${OUTPUT_DIR}/${base_name}_1.fastq.gz"
        output_R2="${OUTPUT_DIR}/${base_name}_2.fastq.gz"

        #Define report files
	html_report="${OUTPUT_DIR}/${base_name}_report.html"
	json_report="${OUTPUT_DIR}/${base_name}_report.json"

	# Run fastp on the pair of files with overrepresentation analysis
	fastp -i "$input_R1" -I "$input_R2" -o "$output_R1" -O "$output_R2" -n 2 -f 15 -q 20 -l 70 --correction --detect_adapter_for_pe --html "$html_report" --json "$json_report"

	# Check if fastp is executed successfully
	if [ $? -eq 0 ]; then
        	echo "Processed: $base_name"
	else
        	echo "Failed to process: $base_name" >&2
	fi
done

echo "All files have been processed."
