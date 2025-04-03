# This script loops through multiple samples of processed & clean reads and performs reads alignment with the reference genome using STAR

# List of sample names (SRR IDs for each paired-end dataset)
samples=("SRR11266556" "SRR11266555" "SRR11266554" "SRR9038729" "SRR9038731" "SRR9038733" "SRR9038730" "SRR9038732" "SRR9038734" "SRR9690969" "SRR9690970" "SRR9690971" "SRR9690972" "SRR9690973" "SRR9690974")

# Define the temporary directory
tmp_dir="/home/cbr16/Documents/WeeYeZhi/tmp_dir_alignment"

for sample in "${samples[@]}"; do
	echo "Processing sample: $sample"  # Debug statement to indicate which sample is being processed

	# Check if input files exist
	if [[ ! -f /home/cbr16/Documents/WeeYeZhi/output/fastpresults/${sample}_1.fastq.gz || ! -f /home/cbr16/Documents/WeeYeZhi/output/fastpresults/${sample}_2.fastq.gz ]]; then
		echo "Error: Input files for sample $sample not found."
		continue  # Skip to the next sample if input files do not exist
	fi

	# Running STAR for alignment
	STAR --runMode alignReads \
		   --runThreadN 8 \
		   --readFilesIn /home/cbr16/Documents/WeeYeZhi/output/fastpresults/${sample}_1.fastq.gz /home/cbr16/Documents/WeeYeZhi/output/fastpresults/${sample}_2.fastq.gz \
		   --readFilesCommand zcat \
		   --genomeDir /home/cbr16/Documents/WeeYeZhi/output/GenomeIndexresults/GenomeIndex/ \
		   --outFileNamePrefix /home/cbr16/Documents/WeeYeZhi/tmp_dir_alignment/${sample}_ \
		   --outSAMtype BAM SortedByCoordinate \
		   --outSAMunmapped Within \
		   --outSAMattributes Standard \
		   --outTmpDir ${tmp_dir}${sample}_tmp/  # Specify unique tmp dir for each sample

	# Check if STAR executed successfully
	if [[ $? -ne 0 ]]; then
		echo "Error: STAR alignment failed for sample $sample."
	else
		echo "Successfully processed sample: $sample."
	fi
done
