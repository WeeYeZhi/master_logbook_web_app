from PIL import Image
from pathlib import Path
import requests
import streamlit as st
import pandas as pd
from pygments.lexers.sql import language_re
from streamlit_lottie import st_lottie
from streamlit_option_menu import option_menu
from streamlit_timeline import timeline
import matplotlib.pyplot as plt

# Find more emojis here: https://www.webfx.com/tools/emoji-cheat-sheet/
# Find more animations here: https://lottiefiles.com/search?category=animations&utm_source=search&utm_medium=platform

st.set_page_config(page_title="My LogBook", page_icon="📚", layout="wide")

#----PATH SETTINGS----
current_dir = Path(__file__).parent if "__file__" in locals() else Path.cwd()
falco1_file = current_dir / "assets" / "falco1.sh"
falco2_file = current_dir / "assets" / "falco2.sh"
gzip_file = current_dir / "assets" / "gzip.sh"
fastp_file = current_dir / "assets" / "fastp.sh"
genomeestimate_file = current_dir / "assets" / "genome_estimate.py"
longstitch_file = current_dir / "assets" / "run_longstitch.sh"
braker3_file = current_dir / "assets" / "braker3_perl_module_installation.sh"
star_file = current_dir / "assets" / "RNAseq_alignment_with_STAR.sh"
deseq2rmd_file = current_dir / "assets" / "deseq2.Rmd"
CPB_pic = current_dir / "assets" / "CPB.png"

# ---- HEADER SECTION ----
with st.container():
    left_column, right_column = st.columns((1, 1))
    with left_column:
        st.title("Identification of inhibitors against cocoa pod borer (*Conopomorpha cramerella*) developmental proteins using bioinformatics approach")
        st.header("LogBook 👨🔬‍")
        st.write("###")
    with right_column:
        CPB_pic = Image.open(CPB_pic)
        st.image(CPB_pic, width=600)

# ----SIDE BAR MENU ----
with st.sidebar:
    selected = option_menu(
        menu_title="Methodology",
        options=["Phase 1: Sequence-Based Analysis", "Phase 2: Structure-Based Analysis", "Phase 3: Molecular Docking & Dynamics Simulation", "Additional Notes"],
    )

#----CONTENT SECTION----

# Phase 1: Sequence-Based Analysis

if selected == "Phase 1: Sequence-Based Analysis":
    with st.container():
        st.write("---")
        st.header("Sequence-Based Analysis 🧬")
        st.write("###")
        st.write("**1. code as the substituted cbr15 user (avoid working in the root as it will damage the OS)**")
        st.write("✔️ change from the root user to the cbr15 user")
        st.code("su cbr15", language="bash")
        st.code("cd ~", language="bash")
        st.code("export CONDA_PREFIX='/media/Raid/Home/anaconda3'", language="bash")
        st.code("export PATH='/media/Raid/Home/anaconda3/bin:$PATH'", language="bash")
        st.code("conda activate base", language="bash") # after activating the base environment, you can begin to create the other envs and install the required bioinformatics tools

        st.write('###')

        st.write("**2. Download the raw sequencing data of CPB from NCBI SRA database**")
        st.write("✔️install the sra-toolkit within the NGSWee environment")
        st.code("sudo apt -y install sra-toolkit", language='bash')
        st.write("✔️prefetch all the 15 .sra files by using the prefetch tool available in the sra-toolkit")
        st.code(
            "prefetch SRR11266556 SRR11266555 SRR11266554 SRR9038729 SRR9038731 SRR9038733 SRR9038730 SRR9038732 SRR9038734 SRR9690969 SRR9690970 SRR9690971 SRR9690972 SRR9690973 SRR9690974",
            language='bash')
        st.write("✔️convert all the 15 files extension one by one from .sra to .fastq format")
        st.code("fasterq-dump SRR11266556 SRR11266555 SRR11266554 SRR9038729 SRR9038731 SRR9038733 SRR9038730 SRR9038732 SRR9038734 SRR9690969 SRR9690970 SRR9690971 SRR9690972 SRR9690973 SRR9690974",language='bash')
        st.write("✔️zip the two large .fastq files (paired-end sequencing data) for each sample into one .gz file by using the gzip bash script")
        # ----LOAD GZIP BASH SCRIPT----
        # Check if the file exists before reading
        if gzip_file.exists():
            with open(gzip_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Gzip Bash Script",
                data=script_byte,
                file_name=gzip_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{gzip_file.name} does not exist.")

        st.write("###")

        st.write("**3. Check the base quality of the 30 raw fastq files by using Falco**")
        st.write("✔️install Falco within Linux terminal")
        st.code("conda install -c bioconda falco", language='bash')
        st.write("✔️check the base quality of all the .fastq.gz files one by one by using the falco bash script")
        # ----LOAD FALCO BASH SCRIPT----
        # Check if the file exists before reading
        if falco1_file.exists():
            with open(falco1_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Falco Bash Script",
                data=script_byte,
                file_name=falco1_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{falco1_file.name} does not exist.")

        st.write("###")

        st.write("**4. Trim and clean the 30 raw reads using fastp**")
        st.write("✔️install fastp within Linux terminal")
        st.code("sudo apt -y install fastp", language='bash')
        st.write(
            "✔️remove low quality reads with Phred score < 30, remove short reads with length < 70bp, remove adapters & remove ambiguous bases (N) up to a maximum of 2 by using the fastp bash script")
        # ----LOAD FASTP BASH SCRIPT----
        # Check if the file exists before reading
        if fastp_file.exists():
            with open(fastp_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Fastp Bash Script",
                data=script_byte,
                file_name=fastp_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{fastp_file.name} does not exist.")

        st.write("###")

        st.write("**5. Check the base quality of the 30 trimmed .fastq.gz files once again using falco**")
        st.write("✔️check the base quality of the trimmed files using the bash script")
        # ----LOAD FALCO2 BASH SCRIPT----
        # Check if the file exists before reading
        if falco2_file.exists():
            with open(falco2_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Falco Bash Script",
                data=script_byte,
                file_name=falco2_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{falco2_file.name} does not exist.")

        st.write("###")

        st.write("**6. Check the base quality of the two short illumina .fastq files provided by LKM by using Falco**")
        st.code("falco Conopomorpha_raw_1.fastq.gz", language="bash")
        st.code("falco Conopomorpha_raw_2.fastq.gz", language="bash")

        st.write("###")

        st.write("**7. Trim and clean the two short illumina .fastq files using fastp**")
        st.code("nohup fastp -i /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1.fastq.gz -I /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2.fastq.gz -o /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq.gz -O /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq.gz -n 2 -f 15 -q 20 -l 70 --correction --detect_adapter_for_pe --html /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/fastp_report.html --json /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/fastp_report.json > fastp_output.log 2>&1 &", language="bash") # apply relaxed fastp parameter by setting -l 70 (originally wanted to retain more reads, but seems like a lot of sequencing errors and kmer artifacts remained)
        st.code("fastp -i /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1.fastq.gz -I /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2.fastq.gz -o /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/trimmed_Conopomorpha_1.fastq.gz -O /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/trimmed_Conopomorpha_2.fastq.gz -n 2 -f 15 -q 20 -l 150 --correction --detect_adapter_for_pe --html /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/fastp_report.html --json /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/fastp_report.json", language="bash") # apply stricter fastp parameter to trim more low-quality short reads by setting -l 150 to remove sequencing errors and kmer artifacts
        st.code("nohup fastp -i /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1.fastq.gz -I /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2.fastq.gz -o /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/trimmed_Conopomorpha_1.fastq.gz -O /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/trimmed_Conopomorpha_2.fastq.gz -n 2 -f 15 -q 20 -l 200 --correction --detect_adapter_for_pe --html /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/fastp_report.html --json /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results_stricter/fastp_report.json > fastp_output.log 2>&1 &", language="bash") # apply even stricter fastp parameter to trim more low-quality short reads by setting -l 200 to remove more sequencing errors and kmer artifacts (make only high-quality reads are retained)

        st.write("###")

        st.write("**8. Check the base quality of the two short Illumina .fastq files using Falco once again**")
        st.code("falco trimmed_Conopomorpha_1.fastq.gz", language="bash")
        st.code("falco trimmed_Conopomorpha_2.fastq.gz", language="bash")

        st.write("###")

        st.write("**9. Alternatively, you can also trim the two short Illumina .fastq files using Cutadapt**")
        st.code("nohup cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC -q 20,20 -m 70 --max-n 2 -j 48 --poly-a --no-indels --trim-n --report full -o trimmed_Conopomorpha_1.fastq.gz -p trimmed_Conopomorpha_2.fastq.gz /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1.fastq.gz /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2.fastq.gz > cutadapt.log 2>&1 &", language="bash")
        st.code("nohup cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGC -q 20,20 -m 150 --max-n 2 -j 48 --poly-a --no-indels --trim-n --report full -o trimmed_Conopomorpha_1.fastq.gz -p trimmed_Conopomorpha_2.fastq.gz /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1.fastq.gz /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2.fastq.gz > cutadapt.log 2>&1 &", language="bash")

        st.write("###")

        st.write("**10. Estimate the haploid genome size of the insect using jellyfish**")
        st.write("✔️ perform k-mer counting")
        st.code("nohup jellyfish count -m 19 -s 50G -t 48 -C /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq -o k_19_mer_counts.jf > jellyfish_k19_output.log 2>&1 &", language="bash")
        st.code("nohup jellyfish count -m 21 -s 50G -t 48 -C /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq -o k_21_mer_counts.jf > jellyfish_k21_output.log 2>&1 &", language="bash")
        st.code("nohup jellyfish count -m 22 -s 50G -t 48 -C /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq -o k_22_mer_counts.jf > jellyfish_k22_output.log 2>&1 &", language="bash")
        st.code("nohup jellyfish count -m 27 -s 50G -t 48 -C /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq -o k_27_mer_counts.jf > jellyfish_k27_output.log 2>&1 &", language="bash")
        st.code("nohup jellyfish count -m 31 -s 50G -t 48 -C /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq -o k_31_mer_counts.jf > jellyfish_k31_output.log 2>&1 &", language="bash")
        st.write("✔️ generate a k-mer frequency histogram")
        st.code("jellyfish histo k_19_mer_counts.jf > k_19_mer_counts.histo", language="bash")
        st.code("jellyfish histo k_21_mer_counts.jf > k_21_mer_counts.histo", language="bash")
        st.code("jellyfish histo k_22_mer_counts.jf > k_22_mer_counts.histo", language="bash")
        st.code("jellyfish histo k_27_mer_counts.jf > k_27_mer_counts.histo", language="bash")
        st.code("jellyfish histo k_31_mer_counts.jf > k_31_mer_counts.histo", language="bash")
        st.write("✔️ estimate haploid genome size of the insect")
        st.code("python3 genome_estimate.py", language="bash") # you can check the version of python using 'python3 --version'
        # ----LOAD  BASH SCRIPT----
        # Check if the file exists before reading
        if genomeestimate_file.exists():
            with open(genomeestimate_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Genome Estimate Script",
                data=script_byte,
                file_name=genomeestimate_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{falco2_file.name} does not exist.")

        st.write("###")

        st.write("**11. Check & compare the quality of the raw & processed LKM PacBio read using NanoPlot**")
        st.code("nohup NanoPlot -t 48 --fastq /media/Raid/Wee/WeeYeZhi/resources_from_LKM/pacbio_long_read/PacBio.fq.gz --info_in_report --plots dot kde --legacy hex -o nanoplot_processed_pacbio_read > nanoplot.log 2>&1 &", language="bash") # check the quality of the raw PacBio read
        st.code("nohup NanoPlot -t 48 --fastq /media/Raid/Wee/WeeYeZhi/processed_pacbio/proovread/proovread/siamaera_output_2.fq --info_in_report --plots dot kde --legacy hex -o nanoplot_processed_pacbio_read > nanoplot.log 2>&1 &", language="bash") # check the quality of the processed PacBio read
        st.markdown("[Visit NanoPlot Github Documentation](https://github.com/wdecoster/NanoPlot?tab=readme-ov-file)")
        st.markdown("[Read & Cite NanoPlot Publication](https://doi.org/10.1093/bioinformatics/btad311)")

        st.write("###")

        st.write("12. If the quality of the PacBio read is bad, then proovread 2.14.1 will be used to correct the PacBio read using the short Illumina read")
        st.code("conda create -n proovread", language="bash") # use the cbr15 INBIOSIS computer (no need to su cbr15)
        st.code("conda activate proovread", language="bash")
        st.code("git clone https://github.com/BioInf-Wuerzburg/proovread", language="bash")
        st.code("cd proovread", language="bash")
        st.code("make", language="bash")
        st.code("sudo apt-get install liblog-log4perl-perl", language="bash") # make sure that you install all the required dependencies of proovread for it to run successfully without encountering any errors later
        st.code("sudo apt-get install libfile-which-perl", language="bash")
        st.code("./bin/proovread --help", language="bash") # to verify the installation of proovread
        st.code("nohup ./bin/proovread -t 48 -s /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1/trimmed_Conopomorpha_raw_1.fastq -s /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2/trimmed_Conopomorpha_raw_2.fastq -l /media/Raid/Wee/WeeYeZhi/resources_from_LKM/pacbio_long_read/PacBio.fq -o corrected_PacBio.fq > proovread.log 2>&1 &", language="bash")
        st.code("perl -I /media/Raid/Wee/WeeYeZhi/processed_pacbio/proovread/lib ./bin/SeqFilter --defaults", language="bash") # display all the available command lines (help) of SeqFilter
        st.code("perl ./bin/siamaera --help", language="bash") # display all the available command lines (help) of siamaera
        st.code("nohup perl ./bin/SeqFilter --in proovread/proovread.untrimmed.fq --min-length 200 --substr proovread/proovread.chim.tsv --out proovread/seqfilter_output.fq > proovread/seqfilter_output.log 2>&1 &", language="bash") # run the SeqFilter step to trim the reads
        st.code("nohup /usr/bin/env perl ./bin/siamaera < proovread/seqfilter_output.fq > proovread/siamaera_output.fq 2> proovread/siamaera_output.log &") # run the siamaera step (final step of proovread) to further trim the chimeric reads
        st.markdown("[Visit proovread GitHub Page](https://github.com/BioInf-Wuerzburg/proovread/blob/master/README.org)")
        st.markdown("[Visit proovread GitHub Poster](https://github.com/BioInf-Wuerzburg/proovread/blob/master/media/proovread-poster.pdf)")
        st.markdown("[Visit proovread Publication](https://doi.org/10.1093/bioinformatics/btu392)")
        st.markdown("[Visit bwa Publication](https://doi.org/10.1093/bioinformatics/bts280)")
        st.markdown("[Visit BLASR Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-238)")
        st.markdown("[Visit SHRiMP Publication](https://doi.org/10.1371/journal.pcbi.1000386)")


        st.write("###")

        st.write("**13. Perform genome assembly to construct a complete reference genome of CPB using SPAdes (assemble short reads and long reads together)**")
        st.write("✔️create a virtual environment called SPAdes")
        st.code("conda create -n SPAdes", language="bash")
        st.write("✔️activate the SPAdes environment")
        st.code("conda activate SPAdes", language="bash")
        st.write("✔️install SPAdes within Linux terminal")
        st.code("wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz", language="bash")
        st.code("tar -xzf SPAdes-4.0.0-Linux.tar.gz", language="bash")
        st.code("cd SPAdes-4.0.0-Linux/bin/", language="bash")

        # Alternatively, you can install spades using the sudo apt command

        st.code("sudo apt -y install spades", language="bash")
        st.write("✔️test whether you've installed SPAdes successfully")
        st.code("spades.py --version", language="bash")
        st.write("✔️display the list of command options of SPAdes")
        st.code("spades.py -h", language="bash")
        st.write("✔️rename the .fna file to .fasta file")
        st.code("mv GCA_012932125.1_ASM1293212v1_genomic.fna GCA_012932125.1_ASM1293212v1_genomic.fasta", language="bash")
        st.write("✔️run SPAdes to perform genome assembly of CPB step by step")
        st.write("Run the error correction step to generate the corrected reads")
        st.code("nohup spades.py --only-error-correction --pe1-1 /home/cbr16/Documents/WeeYeZhi/resources_from_NCBI/all_R1.fastq.gz --pe1-2 /home/cbr16/Documents/WeeYeZhi/resources_from_NCBI/all_R2.fastq.gz --pe2-1 /home/cbr16/Documents/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1/trimmed_Conopomorpha_raw_1.fastq.gz --pe2-2 /home/cbr16/Documents/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2/trimmed_Conopomorpha_raw_2.fastq.gz --pacbio /home/cbr16/Documents/WeeYeZhi/resources_from_LKM/pacbio_long_read/PacBio.fq.gz -o /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/SPAdes_CPB_genome/error_correction --threads 8 > /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/SPAdes_CPB_genome/error_correction/error_correction.log 2>&1 &", language="bash")
        st.write("Next, use the corrected reads to proceed with the assembly")
        st.code("nohup spades.py --only-assembler --pe1-1 /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/error_correction/corrected/all_R1.fastq.00.0_0.cor.fastq.gz --pe1-2 /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/error_correction/corrected/all_R2.fastq.00.0_0.cor.fastq.gz --pe2-1 /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/error_correction/corrected/trimmed_Conopomorpha_raw_1.fastq.00.0_0.cor.fastq.gz --pe2-2 /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/error_correction/corrected/trimmed_Conopomorpha_raw_2.fastq.00.0_0.cor.fastq.gz --pacbio /home/cbr16/Documents/WeeYeZhi/resources_from_LKM/pacbio_long_read/PacBio.fq.gz --trusted-contigs /home/cbr16/Documents/WeeYeZhi/resources_from_LKM/hybrid_genome_assembly_of_CPB/CPB_insect_draft_assembly.v4.fa -o /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/assembly_output --threads 8 > /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/assembly_output.log 2>&1 &", language="bash")
        st.write("Run the mismatch correction step")
        st.code("nohup spades.py --only-mismatch-correction -o /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/assembly_output --threads 8 > /home/cbr16/Documents/WeeYeZhi/output/SPAdesresults/mismatch_correction.log 2>&1 &", language="bash")
        st.write("Alternatively, if you have a strong HPC, you can directly run the code below to perform genome assembly all at once")
        st.write("set your working directory as /home/cbr16/Documents/WeeYeZhi/ and run SPAdes in this working directory")
        st.code("nohup spades.py --pe1-1 /media/Raid/Wee/WeeYeZhi/concatenated_processed_RNA_seq/all_R1.fastq.gz --pe1-2 /media/Raid/Wee/WeeYeZhi/concatenated_processed_RNA_seq/all_R2.fastq.gz --pe2-1 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1/trimmed_Conopomorpha_raw_1.fastq.gz --pe2-2 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2/trimmed_Conopomorpha_raw_2.fastq.gz --pacbio /media/Raid/Wee/WeeYeZhi/resources_from_LKM/pacbio_long_read/PacBio.fq.gz --trusted-contigs /media/Raid/Wee/WeeYeZhi/resources_from_LKM/hybrid_genome_assembly_of_CPB/CPB_insect_draft_assembly.v4.fa -o /media/Raid/Wee/WeeYeZhi/output/SPAdesresults/SPAdes_CPB_genome/assembly_output --threads 48 --careful > /media/Raid/Wee/WeeYeZhi/output/SPAdesresults/SPAdes_CPB_genome/assembly_output/output.log 2>&1 &", language="bash")
        st.markdown("[Visit SPAdes GitHub Page](https://github.com/ablab/spades/blob/v4.0.0/README.md)")
        st.markdown("[Visit SPAdes Assembly Toolkit](https://ablab.github.io/spades/)")
        st.markdown("[Read SPAdes De Novo Assembler Publication](https://doi.org/10.1002/cpbi.102)")
        st.markdown("[Read HybridSPAdes Publication](https://doi.org/10.1093/bioinformatics/btv688)")

        st.write("---")

        st.write("Run 3 iterations of the SPAdes pipeline with varying k-mer sizes, resulting in three different assembly versions 'v1.0', 'v1.1' and 'v1.2' (make 3 different assembly folders for 'v1.0', 'v1.1' and 'v1.2'")
        st.code("nohup spades.py --pe1-1 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1.fastq.gz --pe1-2 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2.fastq.gz --pacbio /media/Raid/Wee/WeeYeZhi/resources_from_LKM/pacbio_long_read/PacBio.fq.gz -o /media/Raid/Wee/WeeYeZhi/output/SPAdesresults/SPAdes_hybrid_genome_assembly_k213355_1 --threads 32 --memory 400 --careful -k 21,33,55 &", language="bash")
        st.code("nohup spades.py --pe1-1 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq.gz --pe1-2 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq.gz --pacbio /media/Raid/Wee/WeeYeZhi/processed_pacbio/proovread/proovread/siamaera_output_2.fq -o /media/Raid/Wee/WeeYeZhi/output/SPAdesresults/SPAdes_hybrid_genome_assembly_k21335577 --threads 32 --memory 400 --careful -k 21,33,55,77 &", language="bash")
        st.code("nohup spades.py --pe1-1 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_1/trimmed_Conopomorpha_1.fastq.gz --pe1-2 /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_fastp_2/trimmed_Conopomorpha_2.fastq.gz --pacbio /media/Raid/Wee/WeeYeZhi/processed_pacbio/proovread/proovread/siamaera_output_2.fq -o /media/Raid/Wee/WeeYeZhi/output/SPAdesresults/SPAdes_hybrid_genome_assembly_k127 --threads 32 --memory 400 --careful -k 127 &", language="bash")
        st.write("❗Just in case if spades crashes due to insufficient memory, you can resume the spades run by restarting the run from the last checkpoint/stage")
        st.code("nohup spades.py -t 32 --memory 450 --restart-from last -o /media/Raid/Wee/WeeYeZhi/output/SPAdesresults/SPAdes_raw_hybrid_genome_assembly_k213355 > spades_continue.log 2>&1 &", language="bash")
        st.write("---")

        st.write("**Additional note**")
        st.write("❗SPAdes is a versatile toolkit designed for assembly and analysis of sequencing data. SPAdes is primarily developed for Illumina sequencing data, but can be used for IonTorrent as well. Most of SPAdes pipelines support hybrid mode, i.e. allow using long reads (PacBio and Oxford Nanopore) as a supplementary data.")
        st.write("❗Only files with extension .fq, .fastq, .bam, .fa, .fasta, .fq.gz, .fastq.gz, .bam.gz, .fa.gz, .fasta.gz are supported by SPAdes")
        st.write("---")

        st.write("###")

        st.write("**14. Index the genome assembly by using bwa-mem2**")
        st.write("✔️create a virtual environment called bwa-mem2")
        st.code("conda create -n bwa-mem2", language="bash")
        st.write("✔️activate the bwa-mem2 environment")
        st.code("conda activate bwa-mem2", language="bash")
        st.write("✔️install bwa-mem2 within the environment")
        st.code("conda install -c bioconda bwa-mem2", language="bash")
        st.write("✔️index the genome assembly provided by LKM to allow efficient alignment of the genome assembly with short reads later on")
        st.code("nohup bwa-mem2 index /media/Raid/Wee/WeeYeZhi/resources_from_LKM/hybrid_genome_assembly_of_CPB/CPB_insect_draft_assembly.v4.fa > bwa_index.log 2>&1 &")
        st.markdown("[Visit BWA GitHub Page](https://github.com/lh3/bwa/blob/master/README.md)")
        st.markdown("[Visit BWA-MEM2 GitHub Page](https://github.com/bwa-mem2/bwa-mem2)")

        st.write("###")

        # Continue working in the same 'bwa-mem2' environment

        st.write("**15. Correct the genome assembly provided by LKM using short reads**")
        st.write("✔️align short genomic reads with the genome assembly using bwa-mem2 to produce a SAM file")
        st.code("nohup bwa-mem2 mem -t 48 /media/Raid/Wee/WeeYeZhi/resources_from_LKM/hybrid_genome_assembly_of_CPB/CPB_insect_draft_assembly.v4.fa /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_1/trimmed_Conopomorpha_raw_1.fastq /media/Raid/Wee/WeeYeZhi/output/Illumina_reads_LKM/fastp_results/Conopomorpha_raw_2/trimmed_Conopomorpha_raw_2.fastq > CPB_raw_hybrid_assembly_output.sam 2> CPB_raw_hybrid_assembly_output_sam.log &", language="bash")
        st.write("✔️install samtools within the 'samtools' conda environment")
        st.code("conda install bioconda::samtools", language="bash")
        st.write("✔️convert the SAM file into BAM file")
        st.code("nohup samtools view -@ 48 -Sb -o /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/samtobam_conversion/CPB_raw_hybrid_assembly_output.bam /media/Raid/Wee/WeeYeZhi/output/bwa-mem2results/CPB_raw_LKM_hybrid_assembly/short_read_alignment/CPB_raw_hybrid_assembly_output.sam > /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/samtobam_conversion/CPB_raw_hybrid_assembly_output_samtools_view.log 2>&1 &", language="bash")
        st.write("✔️sort the BAM file based on its genomic coordinates")
        st.code("nohup samtools sort -@ 48 -O bam -o /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/bam_sorting/CPB_raw_hybrid_assembly_output_sorted.bam /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/samtobam_conversion/CPB_raw_hybrid_assembly_output.bam > /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/bam_sorting/CPB_raw_hybrid_assembly_output_samtools_sort.log 2>&1 &")
        st.write("✔️index the BAM file")
        st.code("nohup samtools index /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/bam_sorting/CPB_raw_hybrid_assembly_output_sorted.bam > /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/bam_indexing/CPB_raw_hybrid_assembly_output_samtools_index.log 2>&1 &", language="bash")
        st.write("✔️generate the statistics summary of the sorted BAM file")
        st.code("nohup samtools flagstat output_sorted.bam > output_samtools_flagstat.log 2>&1 &", language="bash")
        st.code("nohup samtools coverage -o coverage.txt output_sorted.bam > output_samtools_coverage.log 2>&1 &", language="bash")
        st.code("nohup samtools depth output_sorted.bam > samtools_depth.log 2>&1 &", language="bash")
        st.write("✔️install pilon within the 'pilon' conda environment")
        st.code("java -version", language="bash") #check whether your conda environment has Java installed already
        st.code("sudo apt update", language="bash") #update the conda environment
        st.code("sudo apt install openjdk-11-jdk", language="bash") #install the Java Development Kit (JDK) 11 if your conda environment doesn't have Java installed
        st.code("wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar -O pilon.jar", language="bash") # Download and install Pilon
        st.code("java -jar /media/Raid/Wee/WeeYeZhi/output/Pilonresults/pilon.jar --version", language="bash") # After installing Pilon, check its version
        st.write("✔️polish the genome assembly by correcting SNPs, small insertions, deletions (indels) and large structural variations using short reads with Pilon")
        st.code("nohup java -Xmx400G -jar /media/Raid/Wee/WeeYeZhi/output/Pilonresults/pilon.jar --genome /media/Raid/Wee/WeeYeZhi/resources_from_LKM/hybrid_genome_assembly_of_CPB/CPB_insect_draft_assembly.v4.fa --fix all --changes --bam /media/Raid/Wee/WeeYeZhi/output/samtoolsresults/CPB_raw_LKM_hybrid_assembly/bam_sorting/CPB_raw_hybrid_assembly_output_sorted.bam --output polished > pilon_stdout.log 2> pilon_stderr.log &") # try to allocate more memory for Pilon to run to avoid encountering OutofMemoryError
        st.markdown("[Visit samtools GitHub Page](https://github.com/rnnh/bioinfo-notebook/blob/master/docs/samtools.md)")
        st.markdown("[Visit Pilon GitHub Page](https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage)")
        st.markdown("[Visit Pilon Step-by-step Installation Page](https://kalonjilabs.com/posts/How-to-Install-Pilon/)")


        st.write("###")

        st.write("**16. Align long read with the genome assembly to output the alignment statistics later (Optional) (No need to run this before running Longstitch as Longstitch already incorporated the minimap2 pipeline)**")
        st.write("✔️create a new environment named as 'minimap2' and activate it")
        st.code("conda create -n minimap2", language="bash")
        st.code("conda activate minimap2", language="bash")
        st.write("✔️install minimap2 within the 'minimap2' environment")
        st.code("git clone https://github.com/lh3/minimap2", language="bash")
        st.code("cd minimap2 && make", language="bash")
        st.write("✔️check the version of minimap2")
        st.code("./minimap2 --version", language="bash") #the official version of minimap2 is 2.28
        st.write("✔️run minimap2 to align the long PacBio CLR read with the insect genome assembly")
        st.code("nohup ./minimap2 -ax map-pb your_genome_assembly_polished_by_Pilon.fa pacbio.fq.gz -t 48 > aln.sam 2> minimap2.log &", language="bash")
        st.markdown("[Visit minimap2 GitHub Page](https://github.com/lh3/minimap2?tab=readme-ov-file)")
        st.markdown("[Visit minimap2 User Manual Page](https://lh3.github.io/minimap2/minimap2.html)")

        st.write("###")

        st.write("**17. Correct the genome assembly provided by LKM further using long read**")
        st.write("✔️create a new environment named as 'longstitch' with a compatible version of Python to install longstitch")
        st.code("conda create -n longstitch", language="bash")
        st.code("conda install -c bioconda longstitch", language="bash")
        st.write("✔️check the version of Longstitch, tigmint and ntLink that you are using respectively within the 'longstitch' conda environment")
        st.code("longstitch", language="bash")
        st.write("✔️check whether all the dependencies of Longstitch have been installed successfully within the 'longstitch' conda environment") # you can also use this to check the version of each dependency
        st.code("conda list | grep make", language="bash")
        st.code("conda list | grep tigmint", language="bash")
        st.code("conda list | grep ntlink", language="bash")
        st.code("conda list | grep abyss", language="bash")
        st.code("conda list | grep arcs", language="bash")
        st.code("conda list | grep links", language="bash")
        st.code("conda list | grep samtools", language="bash")
        st.write("✔️run longstitch (tigmint-long & ntLink) pipeline through bash script to scaffold the CPB genome assembly")
        st.write("❗make sure that you place your executable bash script, genome assembly file & long read file in the same working directory before you run the bash script")
        st.code("chmod +x run_longstitch.sh", language="bash")
        st.code("nohup bash run_longstitch.sh > longstitch_output.log 2>&1 &", language="bash") # run bash script in the background
        st.write("✔️check whether the bash script is still running normally in the background")
        st.code("ps aux | grep run_longstitch.sh", language="bash")
        st.write("✔️if there's need for you to terminate the bash script")
        st.code("pkill -f -9 run_longstitch.sh", language="bash")
        # ----LOAD  BASH SCRIPT----
        # Check if the file exists before reading
        if longstitch_file.exists():
            with open(longstitch_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Longstitch Script",
                data=script_byte,
                file_name=longstitch_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{longstitch_file.name} does not exist.")

        st.markdown("[Visit LongStitch GitHub Page](https://github.com/bcgsc/LongStitch/blob/master/README.md)")
        st.markdown("[Visit tigmint GitHub Page](https://github.com/bcgsc/tigmint/blob/master/README.md)")
        st.markdown("[Visit ntLink GitHub Page](https://github.com/bcgsc/ntLink/blob/master/README.md)")
        st.markdown("[Visit arks GitHub Page](https://github.com/bcgsc/arcs/blob/master/README.md)")

        st.write("###")

        st.write("**18. Evaluate the completeness of the genome assembly after scaffolding**")
        st.write("✔️create a virtual environment called 'busco' & install BUSCO within the environment")
        st.code("conda create -n busco busco -c bioconda -c conda-forge -c defaults", language="bash")
        st.write("try to run wget or git clone to get the latest version of busco from official github page")
        st.write("✔️activate the BUSCO environment")
        st.code("conda activate busco", language="bash")
        st.write("✔️determine the lineage file suitable to be used for CPB genome by listing the lineage datasets available in BUSCO first, followed by referring to the NCBI BioProject of CPB (taxonomy) to check its taxonomy")
        st.code("busco --list-datasets", language="bash")
        st.code("busco --list-datasets | grep -i lepidoptera_odb12", language="bash")
        st.write("✔️run BUSCO to evaluate the completeness of the improved CPB genome")
        st.code("export NUMEXPR_MAX_THREADS=48", language="bash")
        st.code("nohup busco -m genome -i /media/Raid/Wee/WeeYeZhi/resources_from_LKM/hybrid_genome_assembly_of_CPB/CPB_insect_draft_assembly.v4.fa -c 48 -l lepidoptera_odb12 -o CPB_raw_hybrid_assembly_busco > CPB_raw_hybrid_assembly_busco_output.log 2>&1 &", language="bash")
        st.markdown("[Visit BUSCO User Guide Page](https://busco.ezlab.org/busco_userguide.html)")

        st.write("###")

        st.write("❗After running Pilon, remember to remove the word 'pilon' from the list of scaffold sequences within the polished.fasta file and save the polished.fasta file as polished_modified.fasta. After that, only then, you proceed with the next analysis like BUSCO analysis for example")
        st.code("sed 's/|pilon//' polished.fasta > polished_modified.fasta")

        st.write("###")

        st.write("**19. Evaluate the quality of the genome assembly after scaffolding**")
        st.write("✔️create a virtual environment called 'quast'")
        st.code("conda create -n quast", language="bash")
        st.write("✔️activate the 'quast' environment")
        st.code("conda activate quast", language="bash")
        st.write("✔️install quast within the Linux terminal")
        st.code("conda install -c bioconda quast")
        st.write("✔️display the command-line options available within the quast")
        st.code("quast.py -h", language="bash")
        st.write("✔️run quast to evaluate the quality of the draft assembly and the scaffolded assembly")
        st.code("nohup quast.py CPB_assembly.fa --report-all-metrics --large --eukaryote --threads 48 -o quast_output > quast_output.log 2>&1 &", language="bash")
        st.markdown("[Visit Quast GitHub Page](https://github.com/ablab/quast/blob/master/README.md)")
        st.markdown("[Visit Quast User Manual Page](https://quast.sourceforge.net/docs/manual.html)")

        st.write("###")

        st.write("**20. Perform softmasking for the repeat regions of the genome assembly using RepeatMasker before running BRAKER3**")
        st.write("✔️install RepeatMasker via bioconda")
        st.code("conda install -c bioconda repeatmasker", language="bash") # this should install the other required dependencies like perl & dfam database
        st.write("✔️display the help menu of Repeatmasker to make sure you have installed it correctly")
        st.code("RepeatMasker -h", language="bash")
        st.write("✔️double check whether perl has been installed and install perl if it's not installed")
        st.code("""perl --version
conda install -c bioconda perl""", language="bash")
        st.write("✔️list down all the available repeat libraries used by RepeatMasker & check whether Dfam database is readily available")
        st.code("RepeatMasker -list", language="bash")
        st.write("✔️run RepeatMasker to softmask the genome assembly to convert the repeat regions from uppercase letters to lowercase letters to easily distinguish real genes from repeat regions so that BRAKER3 won't mistake the repeat regions of the genome as coding regions during gene & protein prediction")
        st.code("RepeatMasker -pa 8 -species insect -xsmall -dir masked_output genome.fasta", language="bash") # the xsmall flag ensures softmasking (lowercase repeats)
        st.markdown("[Visit RepeatMasker GitHub Page](https://github.com/Dfam-consortium/RepeatMasker)")
        st.markdown("[Visit RepeatMasker DockerHub Page](https://hub.docker.com/r/dnalinux/repeatmasker)")
        st.markdown("[Visit RepeatMasker User Manual Page](https://www.repeatmasker.org/webrepeatmaskerhelp.html)")
        st.markdown("[Visit RepeatMasker WebServer](https://www.repeatmasker.org/cgi-bin/WEBRepeatMasker)")

        st.write("###")

        st.write("**21. Predict the list of coding genes and proteins of the CPB genome using BRAKER3**")
        st.write("✔️switch to a non-root user")
        st.code("su cbr15", language="bash")
        st.write("✔️after switching, get your user ID & group ID")
        st.code("""
        id -u
        id -g""", language="bash")
        st.write("✔️check the version of docker that you are using")
        st.code("docker --version", language="bash")
        st.write("✔️start running docker desktop and check its status to ensure it's actively running in the background.")
        st.code("systemctl --user start docker-desktop", language="bash")
        st.code("systemctl --user status docker-desktop", language="bash")
        st.write("✔️pull/download the braker3 docker image from DockerHub")
        st.code("docker pull teambraker/braker3", language="bash")
        st.write("✔️check & verify whether you've pulled the braker3 docker image correctly & successfully")
        st.code("docker images | grep braker3", language="bash")
        st.write("✔️run the braker3 docker container as shell")
        st.code("docker run --user 1000:1000 --rm -it -v /media/Raid/Wee/WeeYeZhi/output/braker3:/data teambraker/braker3:latest bash", language="bash")
        st.write("✔️double check & make sure all the perl dependencies of braker3 are already installed inside the braker3 docker container (via bash script or via anaconda environment)")
        # ----LOAD  BASH SCRIPT----
        # Check if the file exists before reading
        if braker3_file.exists():
            with open(braker3_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Braker3 Perl Module Installation Script",
                data=script_byte,
                file_name=braker3_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{braker3_file.name} does not exist.")

        st.code("""
        conda install -c anaconda perl (already pre-instaled inside the braker3 docker container)
        conda install -c anaconda biopython (equivalent to conda/mamba packages: biopython)
        conda install -c bioconda perl-app-cpanminus
        conda install -c bioconda perl-file-spec
        conda install -c bioconda perl-hash-merge (equivalent to libhash-merge-perl)
        conda install -c bioconda perl-module-load-conditional
        conda install -c bioconda perl-posix
        conda install -c bioconda perl-file-homedir
        conda install -c bioconda perl-parallel-forkmanager (equivalent to perl-parallel-forkmanager)
        conda install -c bioconda perl-scalar-util-numeric (equivalent to libscalar-util-numeric-perl)
        conda install -c bioconda perl-yaml (equivalent to libyaml-perl)
        conda install -c bioconda perl-class-data-inheritable (equivalent to libclass-data-inheritable-perl)
        conda install -c bioconda perl-exception-class (equivalent to libexception-class-perl)
        conda install -c bioconda perl-test-pod (equivalent to libtest-pod-perl)
        conda install -c bioconda perl-file-which (equivalent to libfile-which-perl)
        conda install -c bioconda perl-mce (equivalent to libmce-perl)
        conda install -c bioconda perl-threaded
        conda install -c bioconda perl-list-util (equivalent to libscalar-list-utils-perl)
        conda install -c bioconda perl-math-utils
        conda install -c bioconda cdbtools
        conda install -c eumetsat perl-yaml-xs
        conda install -c bioconda perl-data-dumper""", language="bash")
        st.write("✔️ensure you have all the perl and python scripts of braker3 ready and make sure all of them are executable (you should get the expected output as shown below)")
        st.code("ls -l *.pl *.py", language="bash")
        st.code("""Expected output
        -rwxr-xr-x 1 katharina katharina  18191 Mai  7 10:25 align2hints.pl
        -rwxr-xr-x 1 katharina katharina   6090 Feb 19 09:35 braker_cleanup.pl
        -rwxr-xr-x 1 katharina katharina 408782 Aug 17 18:24 braker.pl
        -rwxr-xr-x 1 katharina katharina   5024 Mai  7 10:25 downsample_traingenes.pl
        -rwxr-xr-x 1 katharina katharina   5024 Mai  7 10:23 ensure_n_training_genes.py
        -rwxr-xr-x 1 katharina katharina   4542 Apr  3  2019 filter_augustus_gff.pl
        -rwxr-xr-x 1 katharina katharina  30453 Mai  7 10:25 filterGenemark.pl
        -rwxr-xr-x 1 katharina katharina   5754 Mai  7 10:25 filterIntronsFindStrand.pl
        -rwxr-xr-x 1 katharina katharina   7765 Mai  7 10:25 findGenesInIntrons.pl
        -rwxr-xr-x 1 katharina katharina   1664 Feb 12  2019 gatech_pmp2hints.pl
        -rwxr-xr-x 1 katharina katharina   2250 Jan  9 13:55 log_reg_prothints.pl
        -rwxr-xr-x 1 katharina katharina   4679 Jan  9 13:55 merge_transcript_sets.pl
        -rwxr-xr-x 1 katharina katharina  41674 Mai  7 10:25 startAlign.pl""", language="bash")
        st.write("✔️if you don't get the expected output, please make sure all the perl & python scripts are present and executable")
        st.code("chmod a+x *.pl *.py", language="bash")
        st.write("✔️Navigate to the working directory in which all the BRAKER perl scripts reside and add this working directory to your $PATH environment variable (to run your braker3 perl scripts successfully from anywhere)")
        st.code("export PATH=$(pwd):$PATH", language="bash")
        st.write("✔️check if the AUGUSTUS_CONFIG_PATH is already set inside the docker container")
        st.code("docker run --rm teambraker/braker3 bash -c 'echo $AUGUSTUS_CONFIG_PATH' OR echo $AUGUSTUS_CONFIG_PATH")
        st.write("✔️test your braker3 installation and validate whether you have correctly set up the braker3 docker container by testing test3.sh (before running braker3)")
        st.code("bash /opt/BRAKER/example/docker-tests/test3.sh", language="bash") # refer to braker3 github documentation
        st.write("✔️ensure all the mandatory softwares and tools have been installed inside the braker3 docker container")
        mandatory_softwares = ["GeneMark-ETP (already present inside container)",
                               "AUGUSTUS (already present inside container)",
                               "Python3 (already installed by default on Ubuntu)",
                               "Bamtools (already present inside container)",
                               "NCBI BLAST+ or DIAMOND (already present inside container)",
                               "StringTie2 (already present inside container)",
                               "BEDTools (already present inside container)",
                               "GffRead (already present inside container)"]
        for mandatory_software in mandatory_softwares:
            st.write(f"- {mandatory_software}")
        st.write("✔️check whether you need to install all the following optional tools of BRAKER3")
        optional_softwares = ["Samtools (install if you arent sure whether your files are formatted correctly)",
                              "BioPython (already present inside container)",
                              "cdbfasta (already present inside container)",
                              "Spaln (deprecated already, no need to install separately & use it)",
                              "GUSHR",
                              "Tools from UCSC (already present inside container)",
                              "MakeHub (already present inside containerz0",
                              "SRA Toolkit (already present inside container)",
                              "HISAT2 (already present inside container)",
                              "compleasm (need to install manually inside docker container to get the best gene model with highest busco completeness score)",
                              "pandas (install pandas python package manually)",
                              "libc-bin (install libc-bin manually)"]
        for optional_software in optional_softwares:
            st.write(f"- {optional_software}")
        st.write("✔️execute braker3 in the terminal")
        st.code("time perl braker.pl --workingdir=BRAKER3 --genome=/opt/BRAKER/example/gstenome.fa --bam=/opt/BRAKER/example/RNAseq.bam --prot_seq=/opt/BRAKER/example/proteins.fa --AUGUSTUS_BIN_PATH=/usr/bin/ --AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts/ --threads=16  --species=conopomorpha_cramerella", language="bash") # command line from chatgpt, it may not be correct, double check if you want to use it
        st.code("braker.pl --species=conopomorpha_cramerella --genome=genome.fasta --rnaseq_sets_ids=SRA_ID1, SRA_ID2 --makehub --email=weeyz02@gmail.com --threads=16 --rnaseq_sets_dirs=/path/to/local/fastq/files/ --gff3 --ab_initio --busco_lineage=arthropoda_odb10 --workingdir=/your/directory/to/store/braker3/results", language="bash") # run braker3 with provided RNA-seq data (providing SRA IDs & location of fastq files)
        st.code("braker.pl --species=conopomorpha_cramerella --genome=genome.fasta --bam=file.bam --UTR=on --makehub --email=weeyz02@gmail.com --threads=16 --gff3 --ab_initio --busco_lineage=arthropoda_odb10 --workingdir=/your/directory/to/store/braker3/results") # run braker3 with provided aligned RNA-seq data done by HISAT2 tool and with UTR parameter on to predict untranslated regions (UTRs). Running UTR prediction sometimes improves coding sequence prediction accuracy, but not always. If you try this feature, carefully compare results with and without UTR parameters. Please note that we generally assume that bam files were generated with HiSat2 because that is the aligner that would also be executed by BRAKER3 with fastq input. If you want for some reason to generate the bam files with STAR, use the option --outSAMstrandField intronMotif of STAR to produce files that are compatible wiht StringTie in BRAKER3.
        st.write("✔️view & inspect the braker3 gene prediction results using the UCSC Genome Browser")
        st.write("---")
        st.write("**Additional note when running BRAKER3 pipeline**")
        st.write("❗Running UTR prediction with the UTR parameter on (--UTR=on) sometimes improves coding sequence prediction accuracy, but not always. If you try this feature, carefully compare results with and without UTR parameters. Please note that we generally assume that bam files were generated with HiSat2 because that is the aligner that would also be executed by BRAKER3 with fastq input. If you want for some reason to generate the bam files with STAR, use the option --outSAMstrandField intronMotif of STAR to produce files that are compatible wiht StringTie in BRAKER3.")
        st.write("❗Decide whether or not you need to specify the use of the flag, --ab_initio, to run & output both the Augustus's external evidence-based/hints-based prediction result & intrinsic/ab-initio (based on raw input sequence) prediction result & then compare whether using external evidences for instance using RNA-seq data & protein data improves the gene prediction accuracy. You can also save & avoid having computational overhead to not running ab_initio")
        st.write("❗Decide whether or not you need to specify the use of the flag, --makehub --email=weeyz02@gmail.com,  a track data hub for visualizing results with the UCSC Genome Browser will be generated using MakeHub (https://github.com/Gaius-Augustus/MakeHub). This visualization step is highly recommended if you wish to pubish your annotated genome in high-quality journals")
        st.write("❗Decide whether or not you need to specify the use of the flag, --skipAllTraining, to perform only AUGUSTUS predictions, using pre-trained, already existing parameters to save a lot of time & computational resources without the need to train the augustus models. This is used when you dont have any RNA-seq data & protein data as external evidences for your species of interest. This can be faster since it avoids the training step, but the predictions might not be as accurate or species-specific as those obtained after custom training.")
        st.write("❗Decide whether or not you need to specify the use of the flag, --busco_lineage, to select the best gene model with the highest BUSCO completeness scores & lowest BUSCO missing scores with the help & implementation of compleasm. If needed, remember to install compleasm within the BRAKER3 docker container manually")
        st.write("❗The config/ directory from AUGUSTUS can be accessed with the variable AUGUSTUS_CONFIG_PATH. BRAKER3 requires this directory to be in a writable location, so if that is not the case, copy this directory to a writable location, e.g.: cp -r /root/mambaforge/envs/braker3/config/ /absolute_path_to_user_writable_directory/ export AUGUSTUS_CONFIG_PATH=/absolute_path_to_user_writable_directory/config Due to license and distribution restrictions, GeneMark-ETP and ProtHint should be additionally installed for BRAKER3 to fully work. These packages can be either installed as part of the BRAKER3 environment, or the PATH variable should be configured to point to them. The GeneMark key should be located in /root/.gm_key and GENEMARK_PATH should include the path to the GeneMark executables gmes_petap.pl or gmetp.pl.")
        st.write("❗If AUGUSTUS is properly installed and its paths are correctly set in your system's environment variables, you don't need to explicitly specify --AUGUSTUS_BIN_PATH and --AUGUSTUS_SCRIPTS_PATH in the BRAKER command. You can check if AUGUSTUS is properly configured by running:")
        st.write("❗Decide whether or not you need to align RNA-seq data with the genome assembly to output the .bam file using STAR2 or HISAT2")
        st.write("❗Remember to cite all the tools and dependencies of BRAKER3 that you used while executing the BRAKER3 pipeline, don't cite BRAKER3 itself only. Refer to the BRAKER3 GitHub documentation to know how to cite the use of BRAKER3 professionally")
        st.write("❗Check whether you have used NCBI BLAST or DIAMOND to remove redundant gene structures by checking the output braker.log file so that you know which tool you need to cite as BRAKER either uses NCBI BLAST or DIAMOND")
        st.code("which augustus", language="bash")
        st.code("echo $AUGUSTUS_CONFIG_PATH", language="bash")
        st.write("❗To check whether BRAKER3 has been properly installed in the terminal, run the following command:")
        st.code("braker3 --help", language="bash")
        st.code("braker.pl --help", language="bash")
        st.write("❗Refer to the common problems posted at BRAKER3 GitHub Page (https://github.com/Gaius-Augustus/BRAKER#common-problems) and try to get the solutions if you are encountering any errors while executing the BRAKER3 pipeline")
        st.write("❗List of output files produced by BRAKER3 pipeline")
        braker3_outputs = ["braker3.gtf (final set of gene predictions in gtf format that can be further used for downstream analysis, such as RNA-seq analysis, GO enrichment analysis, pathway analysis)",
                            "braker.codingseq (final gene set with coding sequences in FASTA format)",
                            "braker.aa (final gene set with protein sequences in FASTA format)",
                            "braker.gff3 (final set of gene predictions in gff3 format (only produced if the flag --gff3 was specified in the BRAKER3 execution command line) (this annotated genome file produced by BRAKER3 can be further used for downstream analysis, such as visualization of genome annotation or further annotation of genome)",
                            "AUGUSTUS.gtf (set of genes predicted by AUGUSTUS)",
                            "GeneMark.gtf (set of genes predicted by GeneMark)",
                            "hintsfile.gff (the extrinsic evidence data extracted from RNAseq.bam and/or protein data)",
                            "braker_original (set of genes predicted by BRAKER (TSEBRA merge) before compleasm was used to improve BUSCO completeness)",
                            "bbc (output folder of best_by_compleasm.py script from TSEBRA that is used to improve BUSCO completeness in the final output of BRAKER)"]
        for braker3_output in braker3_outputs:
            st.write(f"- {braker3_output}")
        st.write("---")
        st.markdown("[Visit BRAKER3 GitHub Page](https://github.com/Gaius-Augustus/BRAKER/blob/master/README.md)")
        st.markdown("[Visit BRAKER3 Container](https://hub.docker.com/r/teambraker/braker3)")

        st.write("---")
        st.write("###")

        st.write("✔️concatenate all the short reads 1 (r1) and short reads 2 (r2) together respectively (short paired-end reads derived from the NCBI SRA database and short illumina reads provided by LKM")
        st.code("cat SRR*_1.fastq.gz Conopomorpha_raw_1.fastq.gz > all_R1.fastq.gz", language="bash")
        st.code("cat SRR*_2.fastq.gz Conopomorpha_raw_2.fastq.gz > all_R2.fastq.gz", language="bash")
        st.write("✔️use BRAKER3 with the galaxy version")
        st.write("✔️create a new virtual environment called 'braker3'")
        st.code("conda create -n braker3", language="bash")
        st.write("✔️activate the braker3 environment")
        st.code("conda activate braker3", language="bash")
        st.write("✔️install braker3 in the 'braker3' virtual environment")
        st.code("conda install -c bioconda braker3", language="bash") # install BRAKER3
        st.write("✔️check the current directory of the AUGUSTUS config in your system")
        st.code("echo $AUGUSTUS_CONFIG_PATH", language="bash")
        st.write("✔️create a new directory for the AUGUSTUS config in your home directory & copy the AUGUSTUS config to a writable location") # make the AUGUSTUS Config Directory writable
        st.code("cp -r /root/mambaforge/envs/braker3/config/ ~/augustus_config/", language="bash")
        st.write("✔️set the new AUGUSTUS config path")
        st.code("export AUGUSTUS_CONFIG_PATH=~/augustus_config/config", language="bash")
        st.write("✔️verify the new AUGUSTUS config path")
        st.code("echo $AUGUSTUS_CONFIG_PATH", language="bash")
        st.write("✔️install and set up GeneMark-ETP")
        st.markdown("[register & download GeneMark-ETP & get the license key](http://topaz.gatech.edu/GeneMark/license_download.cgi)")
        st.write("✔️extract the downloaded file")
        st.code("tar -xzf gmes_linux_64_4.tar.gz")
        st.code("cd gmes_linux_64_4", language="bash")
        st.write("✔️set GeneMark-ETP path")
        st.code("", language="bash")

        st.write("**22. Check the BUSCO completeness score for the protein sequences predicted by BRAKER3 in the form of braker.aa file**")
        st.code("export NUMEXPR_MAX_THREADS=48", language="bash")
        st.code("nohup busco -m protein -i /path/to/your/braker.aa -c 48 -l lepidoptera_odb12 -o CPB_raw_hybrid_assembly_busco > CPB_raw_hybrid_assembly_busco_output.log 2>&1 &", language="bash")

        st.write("**23. Compute the QUAST metrics for the genome assembly with the additional braker.gff3 file produced by BRAKER3 to know how well the predicted genes are supported by your genome assembly")
        st.code("nohup quast.py /path/to/your/assembled_CPB_genome.fa --report-all-metrics --large --eukaryote --threads 48 -o quast_output --gff /path/to/braker.gff3 > quast_output.log 2>&1 &", language="bash")

# Phase 2: Structure-Based Analysis

if selected == "Phase 2: Structure-Based Analysis":
    with st.container():
        st.write("---")
        st.header("Structure-Based Analysis ⚛")
        st.write("###")
        st.write("**1. Align RNA-seq data with the genome assembly using STAR2**")
        st.write("✔️Create & activate the 'star' conda environment")
        st.code("""
        conda create -n star
        conda activate star""", language="bash")
        st.write("✔️navigate to the desired working directory and clone the star's remote repository")
        st.code("git clone https://github.com/alexdobin/STAR.git", language="bash")
        st.write("✔️install the necessary compilers (g++ compiler, which is a c++ compiler) and build tools (make) within the conda environment")
        st.code("conda install -c conda-forge gxx_linux-64 make", language="bash") #
        st.write("✔️build star from source after cloning the repository")
        st.code("""
        cd STAR/source
        make STAR""", language="bash")
        st.write("✔️double check whether you have installed STAR successfully by checking its version")
        st.code("""
        ls STAR/source/STAR
        STAR --version
        STAR --help""", language="bash")
        st.write("✔️generate the index file of the genome assembly using STAR2")
        st.code("""
        STAR # execute STAR aligner
        --runThreadN 40 # specify the number of threads used for genome indices generation run
        --runMode genomeGenerate # direct STAR to run genome indices generation job
        --sjdbGTFfile Homo_sapiens.GRCh38.97.gtf # specify the path to the file with annotated transcripts in the standard GTF format (produced by BRAKER3)
        --genomeDir GRCh38 # specify the directory to store the output indexed genome file
        --genomeFastaFiles GRCh38.dna.primary.fa # provide your genome of interest in the form of fasta (.fa) format
        --sjdbOverhang 99 # specify the length of the genomic sequence on each side of the splice junction that will be used for constructing the splice junction database during the genome indexing step. In this case, 99 is specified if the maximum length of your RNA-seq data is 100. If the maximum length of your RNA-seq data is 150, please specify 149. This is highly recommended for you to specify. The default value is 100.
        """, language="bash")
        st.write("✔️align the clean RNA-seq data derived from the NCBI SRA database one by one with the indexed genome assembly")
        st.code("""
        STAR # execute STAR aligner
        --runThreadN 40 # specify the number of threads
        --genomeDir GRCh38 # specify the directory where the indexed genome file produced by the previous STAR indexing process is located
        --sjdbGTFfile Homo_sapiens.GRCh38.97.gtf # provide the path to the GTF file for gene annotations (produced by BRAKER3)
        --readFilesIn /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz # input FastQ files for paired-end reads
        --readFilesCommand zcat # specify this flag only if your input RNA-seq reads is compressed. For gzipped files, specify 'readFilesCommand zcat' or 'readFilesCommand gunzip -c'. Whereas, for bzip2-compressed files, specify 'readFilesCommand bunzip2 -c'. Bear in mind that you do not need to specify this flag if your input RNA-seq reads are not compressed.
        --outSAMtype BAM SortedByCoordinate # set the output type to BAM and sort by coordinates
        --outSAMunmapped Within # output unmapped reads within the BAM file
        --outSAMattributes Standard # specify the use of standard attributes to include in the output BAM file
        --quantMode GeneCounts # enable gene count mode for quantification. This is optional (no need to specify if you are going to use featureCounts to count the total number of mapped reads per gene)
        --outFileNamePrefix AlignmentSample1 # define the prefix of the output file names to be "AlignmentSample1..."
        --twopassMode Basic # set the two-pass mode to Basic for alignment (enable STAR to align your input RNA-seq data twice to improve alignment accuracy) (For the most sensitive & novel discovery of splice junctions, it is recommended to run STAR in the 2-pass mode.) 
        """, language="bash")
        st.write("✔️Alternatively, perform batch processing by aligning all the clean RNA-seq data derived from the NCBI SRA database at one single time with the indexed genome assembly")
        # ----LOAD STAR BASH SCRIPT----
        # Check if the file exists before reading
        if star_file.exists():
            with open(star_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download STAR Bash Script",
                data=script_byte,
                file_name=star_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        st.markdown("[Visit STAR GitHub Page](https://github.com/alexdobin/STAR?tab=readme-ov-file)")
        st.markdown("[Visit STAR User Manual Page](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)")

        st.write("###")

        st.write("**2. Alternatively, align RNA-seq data with genome assembly using HISAT2**")
        st.write("✔️Create & activate the 'hisat2' conda environment")
        st.code("""
        conda create -n hisat2
        conda activate hisat2""", language="bash")
        st.write("✔️navigate to the desired working directory and clone the HISAT2's remote repository")
        st.code("git clone https://github.com/DaehwanKimLab/hisat2.git", language="bash")
        st.write("✔️install and build hisat2 tools within the conda environment")
        st.code("""
        cd hisat2
        make
        """, language="bash")
        st.write("✔️build and create the index for the genome assembly using HISAT2 build")
        st.code("""
        hisat2-build # used to build an index for the reference genome assembly so that it can be used later for aligning sequencing reads using HISAT2. This index allows HISAT2 to quickly and efficiently map sequencing reads to the reference genome during the alignment process.
        genome.fa # specify the path to the genome assembly
        genome # specify the prefix of the output indexed genome file
        """, language="bash")
        st.write("✔️After indexing the genome assembly, align one single paired-end read with the reference genome assembly using HISAT2")
        st.code("""
        hisat2 # execute HISAT2 aligner
        -x genome # specify the prefix of the indexed genome file produced by the previous HISAT2-build pipeline
        -1 reads_1.fq # specify the forward read of the paired-end read
        -2 reads_2.fq # specify the reverse read of the paired-end read
        -S output.sam # specify the output file format for the alignment to save the alignment results in the file, 'output.sam'
        """, language="bash")
        st.markdown("[Visit HISAT2 GitHub Page](https://github.com/DaehwanKimLab/hisat2?tab=readme-ov-file)")
        st.markdown("[Visit HISAT2 User Manual Page](https://daehwankimlab.github.io/hisat2/manual/)")
        st.markdown("[Read HISAT2 Publication](https://www.nature.com/articles/s41587-019-0201-4)")
        st.markdown("[Read HISAT Publication](https://www.nature.com/articles/nmeth.3317)")
        st.markdown("[Read this publication to compare between STAR2 and HISAT2](https://ieeexplore.ieee.org/document/10178793)")

        st.write("###")

        st.write("**3. After aligning RNA-seq data with the genome assembly, count the total number of mapped reads per gene using featureCounts**")
        st.code("""
        featureCounts # execute featureCounts, which is read-counting tool from the Subread package that counts how many aligned reads (from BAM files) fall within annotated genomic features (like genes or exons).
        -T 40 # specify the number of threads (CPU cores)
        -p # tell featureCounts that the input BAM files are from paired-end sequencing whereas it counts read pairs as one fragment instead of counting each read separately. (Only specify this if your RNA-seq data is considered paired-end read) (If your read is considered single-end, no need to specify this flag)
        -t exon # tell featureCounts to look for 'feature_type = exon' entry in the GTF file, and count the reads overlapping those exon regions only
        -g gene_id # tell featureCounts to group exons by their gene id, using the gene_id attribute in the GTF file (all exons with the same gene_id will be grouped together)(check your gtf annotation file & decide whether or not you want to group exons by gene_name or gene_id)
        -F GTF # tell featureCounts that the provided genome annotation file in .gtf format
        -a Homo_sapiens.GRCh38.97.gtf # specify the input genome annotation file in the GTF format (produced by previous BRAKER3 pipeline)
        -o countmatrix.txt S1.bam ... Sn.bam # specify the name of the output file as this output file will contain a table with read counts for each gene across all the RNA-seq samples
        """, language="bash")
        st.markdown("[Visit featureCounts User Manual Page](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)")
        st.markdown("[Visit featureCounts Demonstration Video](https://asciinema.org/a/306584?autoplay=1)")
        st.markdown("[Read featureCounts Publication](https://academic.oup.com/bioinformatics/article/30/7/923/232889?login=false)")

        st.write("###")

        st.write("**4. Perform differential expression gene (DEG) analysis using R**")
        st.write("✔️install BiocManager in RStudio") #BiocManager is the official R package that is used to install Bioconductor packages, manage Bioconductor versions, & handle dependencies correctly between CRAN & Bioconductor
        st.code("""if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")""", language="r")
        st.write("✔️install DESeq2, biomaRt, rtracklayer, and GenomicFeatures via BiocManager")
        st.code("""
        BiocManager::install("DESeq2") # perform DEG analysis
        BiocManager::install("biomaRt") # perform gene annotation & ID conversion
        BiocManager::install("rtracklayer") # work with genomic ranges and annotations (like GTF)
        BiocManager::install("GenomicFeatures") # work with genomic ranges and annotations (like GTF)
        """, language="r")
        st.write("✔️install tidyverse and ggrepel via CRAN")
        st.code("""
        install.packages("tidyverse") # manipulate data
        install.packages("ggrepel") # visualize data
        """)
        st.write("✔️load all the installed R libraries into RStudio")
        st.code("""
        library(DESeq2)
        library(biomaRt)
        library(rtracklayer)
        library(GenomicFeatures)
        library(tidyverse)
        library(dplyr)
        library(stringr)
        library(ggplot2)
        library(ggrepel)
        """, language="bash")
        st.write("✔️create a R markdown script in RStudio to run DEG analysis automatically")
        # ----LOAD DESeq2 R SCRIPT----
        # Check if the file exists before reading
        if deseq2rmd_file.exists():
            with open(deseq2rmd_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download DESeq2 R Script",
                data=script_byte,
                file_name=deseq2rmd_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{deseq2rmd_file.name} does not exist.")
        st.markdown("[Visit DESeq2 Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)")
        st.markdown("[Visit DESeq2 GitHub Page](https://github.com/thelovelab/DESeq2)")
        st.markdown("[Visit DESeq2 Tutorial Manual 1](https://lashlock.github.io/compbio/R_presentation.html)")
        st.markdown("[Visit DESeq2 Tutorial Manual 2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start)")
        st.markdown("[Read DESeq2 Publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)")
        st.markdown("[Visit biomaRt Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)")
        st.markdown("[Visit biomaRt GitHub Page](https://github.com/Huber-group-EMBL/biomaRt)")
        st.markdown("[Visit biomaRtr GitHub Page](https://github.com/ropensci/biomartr)")
        st.markdown("[Read biomaRt Publication](https://www.nature.com/articles/nprot.2009.97)")
        st.markdown("[Visit rtracklayer Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)")
        st.markdown("[Visit rtracklayer GitHub Page](https://github.com/lawremi/rtracklayer)")
        st.markdown("[Read rtracklayer Publication](https://academic.oup.com/bioinformatics/article/25/14/1841/225816?login=false)")
        st.markdown("[Visit GenomicFeatures Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)")
        st.markdown("[Visit GenomicFeatures GitHub Page](https://github.com/Bioconductor/GenomicFeatures)")
        st.markdown("[Read GenomicFeatures Publication](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118)")
        st.markdown("[Visit tidyverse Page](https://www.tidyverse.org/packages/)")
        st.markdown("[Visit tidyverse GitHub Page](https://github.com/tidyverse)")
        st.markdown("[Read tidyverse Publication](https://joss.theoj.org/papers/10.21105/joss.01686)")
        st.markdown("[Visit ggrepel Tutorial Manual](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)")
        st.markdown("[Visit ggrepel GitHub Page](https://github.com/slowkow/ggrepel)")
        st.markdown("[Integrated DEG and Pathway Analysis](https://bioinformatics.sdstate.edu/idep/)")
        st.markdown("[IDEP GitHub Page](https://github.com/gexijin/idepGolem)")


# Phase 3: Molecular Docking & Dynamics Simulation

if selected == "Phase 3: Molecular Docking & Dynamics Simulation":
    with st.container():
        st.write("---")
        st.header("Molecular Docking & Dynamics Simulation 🖥️🧪")
        st.write("###")
        st.write("You can use logMD to visualize the trajectory of your protein-ligand complex easily (logMD functions the same as VMD)")
        st.write("generative AI drug design method, DrugHive")
        st.markdown("[Visit the logmd GitHub Page](https://github.com/log-md/logmd)")
        st.markdown("[Try logmd here](https://colab.research.google.com/drive/12adhXXF1MQIzh_vEwKX9r_iF6jV-CNHE#scrollTo=N2_uubn_2qGM)")
        st.markdown("[Try logmd here](https://rcsb.ai/logmd/3d090180)")

# Additional Note

if selected == "Additional Notes":
    with st.container():
        st.write("---")
        st.header("Additional Note ❗")
        st.write("###")
        st.write("1. Kill the process")
        st.code("pkill -9 -f spades",language="bash")  # kill all 'spades' related processes without the need to specify the ID of the process
        st.code("kill -9 103839", language="bash")  # kill the main spades process with the process ID, '103839'

        st.write("###")

        st.write("2. Display all the list of created conda environments")
        st.code("conda info --envs", language="bash")

        st.write("###")

        st.write("3. To check the content of the output.log file whether the bioinformatics pipeline is running in the background, you can execute the code below and exit using Ctrl + C")
        st.write("❗Never press Ctrl + Z as it will stop the process instantly")
        st.code("tail -f output.log", language="bash")
        st.code("jobs -l", language="bash") # only work within the same working terminal, if you close or clear the terminal, you can no longer use this to display the running process in the background
        st.code("more output.log", language="bash")
        st.code("less output.log", language="bash")
        st.code("cat output.log", language="bash")
        st.code("head output.log", language="bash")
        st.code("ps aux | grep longstitch", language="bash")  # check whether the process is still running in the background (doesnt matter whether you already close or clear the terminal, you can still display the running process in the background
        st.code("ps -p 83012", language="bash")  # 83012 is the example ID of your process # to check whether the process is running actively in the background

        st.write("###")

        st.write("4. To remove a conda environment that has been previously created")
        st.code("conda env remove --name name_of_created_environment", language="bash")

        st.write("###")

        st.write("5. To check the version of the bioinformatics tools that you have downloaded within your conda environment")
        st.code("conda list | grep spades", language=" bash")

        st.write("###")

        st.write("6. Delete and clean all the files within the working directory")
        st.code("rm -rf /path/to/your/directory/*", language="bash")

        st.write("###")

        st.write("7. Remember to stop running the Docker desktop if you dont use it to save up CPU resources. Double check whether the Docker desktop already stops running in the background. After checking the status, exit by running 'q'.")
        st.code("systemctl --user stop docker-desktop", language="bash")
        st.code("systemctl --user status docker-desktop", language="bash")
        st.code("q", language="bash")

        st.write("###")

        st.write("8. If you want to show/teach others how to code within Ubuntu/Linux terminal")
        st.code("PS1='$ '", language="bash") # run this command to change how the Linux terminal looks like so that it saves & continues displaying both the previous command lines and results, allowing you to continue coding
        st.code("#", language="bash") # specify this '#" to add a comment before executing any command line/code

        st.write("###")

        st.write("9. End the Ubuntu/Linux terminal session")
        st.code("exit", language="bash")

        st.write("###")

        st.write("10. Clear the terminal")
        st.code("clear", language="bash")
