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
fastqc_raw_RNA_seq_file = current_dir / "assets" / "fastqc_raw_RNA_seq.sh"
fastqc_trimmed_RNA_seq_file = current_dir / "assets" / "fastqc_trimmed_RNA_seq.sh"
gzip_file = current_dir / "assets" / "gzip.sh"
deseq2rmd_file = current_dir / "assets" / "deseq2.Rmd"
gromacs_file = current_dir / "assets" / "Gromacs_codes.txt"
generate_busco_plot_python_script = current_dir / "assets" / "generate_plot.py"
generate_busco_plot_Rscript = current_dir / "assets" / "busco_figure.R"
trimmomatic_file = current_dir / "assets" / "trimmomatic.sh"
trinity_samples_sequenced_by_LKM = current_dir / "assets" / "trinity_samples_sequenced_by_LKM.txt"
trinity_samples_sequenced_by_INBIOSIS = current_dir / "assets" / "trinity_samples_sequenced_by_INBIOSIS.txt"
trinity_samples_sequenced_by_LKM_and_INBIOSIS = current_dir / "assets" / "trinity_samples_sequenced_by_LKM_and_INBIOSIS.txt"
transrate_file = current_dir / "assets" / "transrate.sh"
salmon_transcript_quantification_after_trimmomatic = current_dir / "assets" / "salmon_quantify_after_trimmomatic.sh"
salmon_transcript_quantification_after_trimmomatic_list = current_dir / "assets" / "salmon_transcript_quantification_files_after_trimmomatic.list"
salmon_transcript_quantification_after_cdhitest_trinity = current_dir / "assets" / "salmon_quantify_after_cdhitest_trinity.sh"
salmon_transcript_quantification_after_cdhitest_trinity_list = current_dir / "assets" / "salmon_transcript_quantification_files_after_cdhitest_trinity.list"
salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS = current_dir / "assets" / "salmon_quantify_after_cdhitest_trinity_NCBI_FCS.sh"
salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_list = current_dir / "assets" / "salmon_transcript_quantification_files_after_cdhitest_trinity_NCBI_FCS.list"
salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering = current_dir / "assets" / "salmon_quantify_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.sh"
salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_list = current_dir / "assets" / "salmon_transcript_quantification_files_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.list"
swissprot_fasta_file = current_dir / "assets" / "uniprot_sprot.fasta.gz"
pymol_movie_script = current_dir / "assets" / "movie01_script.pml"
CPB_pic = current_dir / "assets" / "CPB.png"

# ---- HEADER SECTION ----
with st.container():
    left_column, right_column = st.columns((1, 1))
    with left_column:
        st.title("Identification of inhibitors against cocoa pod borer (*Conopomorpha cramerella*) development-related proteins using bioinformatics approach")
        st.header("LogBook 👨🔬‍")
        st.write("###")
    with right_column:
        CPB_pic = Image.open(CPB_pic)
        st.image(CPB_pic, width=600)

# ----SIDE BAR MENU ----
with st.sidebar:
    selected = option_menu(
        menu_title="Methodology",
        options=["Phase 1: Sequence-Based Analysis", "Phase 2: Transcriptomic and Structural-Based Analysis", "Phase 3: Molecular Docking and Dynamics Simulation", "Additional Notes"],
    )

#----CONTENT SECTION----

# Phase 1: Sequence-Based Analysis

if selected == "Phase 1: Sequence-Based Analysis":
    with st.container():
        st.write("---")
        st.header("Phase 1: Sequence-Based Analysis 🧬")

        st.write("###")

        st.write("**1. Download the raw CPB transcriptomic data from NCBI SRA database**")
        st.write("✔️create a 'sratoolkit' conda environment, activate the environment and install sratoolkit via conda")
        st.code("""
        conda create -n sratoolkit
        conda activate sratoolkit
        conda install eyepsychonaut::sra-toolkit
        """, language="bash")
        st.write("✔️prefetch all the 15 .sra files by using the prefetch tool available in the sra-toolkit")
        st.code("""
        prefetch SRR11266556 SRR11266555 SRR11266554 SRR9038729 SRR9038731 SRR9038733 SRR9038730 SRR9038732 SRR9038734 SRR9690969 SRR9690970 SRR9690971 SRR9690972 SRR9690973 SRR9690974
        """,language="bash")
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

        st.write("**2. Check the quality of all the raw transcriptomic reads of CPB derived from the NCBI SRA database via FASTQC**")
        st.write("✔️create a 'fastqc' conda environment, activate the environment and install fastqc via bioconda")
        st.code("""
        conda create -n fastqc
        conda activate fastqc
        conda install bioconda::fastqc
        """, language="bash")
        st.write("✔️verify the installation of fastqc")
        st.code("""
        which fastqc
        fastqc -h
        fastqc -v
        """, language="bash")
        st.write("✔️run FASTQC on the transcriptomic reads of CPB derived from NCBI SRA database")
        st.code("""
        fastqc -f fastq -t 16 -o output_dir SRR11266556_1.fastq.gz # example of fastqc command to run fastqc on a single .fastq file
        """, language="bash")
        st.write("✔️Alternatively, automate the FASTQC run process in large batch mode via bash script")
        st.code("""
        dos2unix fastqc_raw_RNA_seq.sh
        chmod +x fastqc_raw_RNA_seq.sh
        nohup bash fastqc_raw_RNA_seq.sh > fastqc_raw_RNA_seq_from_NCBI_SRA_results_output.log 2>&1 &
        """, language="bash")
        # ----LOAD  BASH SCRIPT----
        # Check if the file exists before reading
        if fastqc_raw_RNA_seq_file.exists():
            with open(fastqc_raw_RNA_seq_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download FASTQC Bash Script",
                data=script_byte,
                file_name=fastqc_raw_RNA_seq_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{fastqc_raw_RNA_seq_file.name} does not exist.")
        st.markdown("[Visit FASTQC GitHub User Manual Page (navigate to the data section)](https://github.com/s-andrews/FastQC/blob/master/fastqc)")
        st.markdown("[Visit FASTQC Conda Installation Page](https://anaconda.org/bioconda/fastqc)")
        st.markdown("[Visit FASTQC Reports & Documentation Page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)")
        st.markdown("[Read and analyze Illumina read with good FASTQC result](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)")

        st.write("###")

        st.write("**3. Aggregate all the FASTQC results of raw CPB transcriptomic reads into one single HTML report for easier data comparison, analysis and visualization via MultiQC**")
        st.write("✔️create a 'multiqc' conda environment, activate the environment and install MultiQC via bioconda")
        st.code("""
        conda create -n multiqc
        conda activate multiqc
        conda install bioconda::multiqc
        """, language="bash")
        st.code("""
        which multiqc
        multiqc -h
        """, language="bash")
        st.write("✔️run MultiQC")
        st.code("""
        multiqc . -o /media/raid/Wee/WeeYeZhi/output/multiqc_raw_RNA_seq_from_NCBI_SRA_results # run this multiqc command within the directory containing all the previous FASTQC results of the reads
        """, language="bash")
        st.markdown("[Visit MultiQC Conda Installation Page](https://anaconda.org/bioconda/multiqc)")
        st.markdown("[Visit MultiQC GitHub Page](https://github.com/MultiQC/MultiQC?tab=readme-ov-file)")
        st.markdown("[Visit MultiQC User Manual Page](https://docs.seqera.io/multiqc/getting_started/running_multiqc)")

        st.write("###")

        st.write("**4. Trim off the adapter sequences from the transcriptomic Illumina paired end read via trimmomatic (if you want to trim it externally first before running trinity pipeline)**")
        st.write("✔️create a 'trimmomatic' conda environment, activate the environment and install trimmomatic via bioconda")
        st.code("""
        conda create -n trimmomatic
        conda activate trimmomatic
        conda install bioconda::trimmomatic
        """, language="bash")
        st.write("✔️verify the installation of trimmomatic")
        st.code("""
        which trimmomatic
        trimmomatic -h
        """, language="bash")
        st.write("✔️run trimmomatic to trim adapter sequences off from a transcriptomic Illumina paired end read")
        st.code("""
        trimmomatic PE -threads 16 input_R1.fastq input_R2.fastq output_R1_paired.fastq output_R1_unpaired.fastq output_R2_paired.fastq output_R2_unpaired.fastq ILLUMINACLIP:/path/to/adapter.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 # run trimmomatic to trim a single Illumina paired end read file only one at a time
        """, language="bash")
        st.write("✔️Alternatively, run trimmomatic to perform adapter trimming process in large batches via bash script")
        st.code("""
        dos2unix trimmomatic.sh
        chmod +x trimmomatic.sh
        nohup bash trimmomatic.sh > trimmomatic_output.log 2>&1 &
        """, language="bash")
        # ----LOAD  BASH SCRIPT----
        # Check if the file exists before reading
        if trimmomatic_file.exists():
            with open(trimmomatic_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Trimmomatic Bash Script",
                data=script_byte,
                file_name=trimmomatic_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{trimmomatic_file.name} does not exist.")
        st.markdown("[Visit Trimmomatic Conda Installation Page](https://anaconda.org/bioconda/trimmomatic)")
        st.markdown("[Visit Trimmomatic GitHub User Manual Page](https://github.com/usadellab/Trimmomatic)")

        st.write("###")

        st.write("**5. Check the quality of all the trimmed transcriptomic reads of CPB derived from the NCBI SRA database via FASTQC**")
        st.write("✔ activate the fastqc environment")
        st.code("""conda activate fastqc""", language="bash")
        st.write("✔️automate the FASTQC run process in large batch mode via bash script")
        st.code("""
               dos2unix fastqc_trimmed_RNA_seq.sh
               chmod +x fastqc_trimmed_RNA_seq.sh
               nohup bash fastqc_trimmed_RNA_seq.sh > fastqc_trimmed_RNA_seq_from_NCBI_SRA_results_output.log 2>&1 &
               """, language="bash")
        # ----LOAD  BASH SCRIPT----
        # Check if the file exists before reading
        if fastqc_trimmed_RNA_seq_file.exists():
            with open(fastqc_trimmed_RNA_seq_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download FASTQC Bash Script",
                data=script_byte,
                file_name=fastqc_trimmed_RNA_seq_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{fastqc_trimmed_RNA_seq_file.name} does not exist.")

        st.write("###")

        st.write("**6. Aggregate all the FASTQC results of trimmed CPB transcriptomic reads into one single HTML report for easier data comparison, analysis and visualization via MultiQC**")
        st.write("✔️activate the multiqc environment")
        st.code("""conda activate multiqc""", language="bash")
        st.write("✔️run multiqc")
        st.code("""
        multiqc . -o /media/raid/Wee/WeeYeZhi/output/multiqc_trimmed_RNA_seq_from_NCBI_SRA_results # run this multiqc command within the directory containing all the previous FASTQC results of the trimmed reads
        """, language="bash")

        st.write("###")

        st.write("**7. Run Trinity to construct CPB transcriptome assembly using reads sequenced by 'LKM','INBIOSIS' and 'LKM and INBIOSIS'**")
        st.write("✔️create a 'trinity' conda environment, activate the environment and install trinity via bioconda")
        st.code("""
        conda create -n trinity
        conda activate trinity
        conda install bioconda::trinity
        """, language="bash")
        st.write("✔️verify the installation of Trinity")
        st.code("""
        Trinity --show_full_usage_info
        which Trinity
        """, language="bash")
        st.write("✔️run Trinity to perform comprehensive transcriptome assembly to assemble all the 3 raw transcriptomic Illumina paired end reads (sequenced by LKM) after trimming the reads via Trimmomatic")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM/trinity_samples_sequenced_by_LKM.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM/trinity_assembly_done_by_LKM_output.log 2>&1 & 
        """, language="bash")
        st.write("✔️run Trinity to perform transcriptome assembly to assemble all the 6 raw transcriptomic Illumina paired end reads (sequenced by INBIOSIS) after trimming the reads via Trimmomatic")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS/trinity_samples_sequenced_by_INBIOSIS.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS/trinity_assembly_done_by_INBIOSIS_output.log 2>&1 & 
        """, language="bash")
        st.write("✔️run Trinity to perform comprehensive transcriptome assembly to assemble all the 9 raw transcriptomic Illumina paired end reads (sequenced by LKM and INBIOSIS)")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_samples_sequenced_by_LKM_and_INBIOSIS.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_CPB_transcriptome_assembly > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_CPB_transcriptome_assembly/trinity_output.log 2>&1 & 
        """, language="bash")
        # ----LOAD TRINITY SAMPLE FILE sequenced BY LKM----
        # Check if the file exists before reading
        if trinity_samples_sequenced_by_LKM.exists():
            with open(trinity_samples_sequenced_by_LKM, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Trinity Sample File Sequenced by LKM",
                data=script_byte,
                file_name=trinity_samples_sequenced_by_LKM.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{trinity_samples_sequenced_by_LKM.name} does not exist.")
        # ----LOAD TRINITY SAMPLE FILE sequenced BY INBIOSIS----
        # Check if the file exists before reading
        if trinity_samples_sequenced_by_INBIOSIS.exists():
            with open(trinity_samples_sequenced_by_INBIOSIS, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Trinity Sample File Sequenced by INBIOSIS",
                data=script_byte,
                file_name=trinity_samples_sequenced_by_INBIOSIS.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{trinity_samples_sequenced_by_INBIOSIS.name} does not exist.")
        # ----LOAD TRINITY SAMPLE FILE sequenced BY LKM & INBIOSIS----
        # Check if the file exists before reading
        if trinity_samples_sequenced_by_LKM_and_INBIOSIS.exists():
            with open(trinity_samples_sequenced_by_LKM_and_INBIOSIS, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Trinity Sample File Sequenced by LKM & INBIOSIS",
                data=script_byte,
                file_name=trinity_samples_sequenced_by_LKM_and_INBIOSIS.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{trinity_samples_sequenced_by_LKM_and_INBIOSIS.name} does not exist.")
        st.markdown("[Visit Trinity Conda Installation Page](https://anaconda.org/bioconda/trinity)")
        st.markdown("[Visit Trinity GitHub Installation Page](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity)")
        st.markdown("[Visit Trinity GitHub Wikipedia Page](https://github.com/trinityrnaseq/trinityrnaseq/wiki)")
        st.markdown("[Visit Trinity GitHub User Manual Page](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)")
        st.markdown("[Visit Trinity User Manual Page](https://stab.st-andrews.ac.uk/wiki/index.php/Trinity)")
        st.markdown("[Read Trinity User Manual Documentation](https://stab.st-andrews.ac.uk/wiki/index.php?title=Trinity&action=pdfbook&format=single)")
        st.markdown("[Visit Trinity User Manual Page 2](https://eagle.fish.washington.edu/whale/fish546/Trinity_r2013-08-14_analysis1-2014-02-08-20-44-13.233/bin/trinityrnaseq_r2013_08_14/docs/advanced_trinity_guide.html)")
        st.markdown("[Visit Trinity GitHub Progress Monitoring Page](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Progress-Monitoring)")
        st.markdown("[Read how to determine the orientation of the transcriptomic read libraries whether it is in the reverse-forward (RF) or forward-reverse (FR) orientation to specify the --SS_lib_type flag](https://knowledge.illumina.com/library-preparation/rna-library-prep/library-preparation-rna-library-prep-reference_material-list/000002238)")
        st.markdown("[Download list of adapter sequences contained within the Illumina paired end read to run trimmomatic for adapter removal purposes](https://github.com/usadellab/Trimmomatic/tree/main/adapters)")
        st.markdown("[Read how to create adapter files to perform adapter removal via trimmomatic](https://www.biostars.org/p/250425/)")

        st.write("###")

        st.write("**8. Cluster highly similar transcripts together to remove redundant and duplicated transcripts from the Trinity assembly via CD-HIT-EST**")
        st.write("✔️create a 'cd-hit' conda environment, activate the environment and install cd-hit via bioconda")
        st.code("""
        conda create -n cd-hit
        conda activate cd-hit
        conda install bioconda::cd-hit
        """, language="bash")
        st.write("✔️verify the installation of cd-hit")
        st.code("""
        which cd-hit
        which cd-hit-est
        """, language="bash")
        st.write("✔️run cd-hit-est to cluster highly similar nucleotide sequences of the Trinity assembly together (avoids inflating transcript counts and improves annotation accuracy.)")
        st.code("""
        nohup cd-hit-est -i /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_assembly_done_by_LKM_and_INBIOSIS.Trinity.fasta -o /media/raid/Wee/WeeYeZhi/output/cd_hit_est_results/cd_hit_est_after_trimmomatic -c 0.95 -n 5 -T 4 > cd_hit_est_after_trimmomatic_output.log 2>&1 & # specify fewer number of threads (use -T 4 in this case) if possible to avoid encountering OutofMemory error
        """, language="bash")
        st.markdown("[Visit CD-HIT Conda Installation Page](https://anaconda.org/bioconda/cd-hit)")
        st.markdown("[Visit CD-HIT GitHub User Manual Page](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.pdf)")
        st.markdown("[Visit CD-HIT User Manual Page](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide)")

        st.write("###")

        st.write("**9. Extract the longest transcript isoform per gene from the transcriptome assembly and use the filtered assembly to compute ExN50 value later**")
        st.write("✔️activate the trinity conda environment")
        st.code("""
        conda activate trinity
        """, language="bash")
        st.write("✔️extract the longest transcript isoform per gene from the transcriptome assembly before & after running CD-HIT-EST")
        st.code("""
        nohup perl /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/get_longest_isoform_seq_per_trinity_gene.pl /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_assembly_done_by_LKM_and_INBIOSIS.Trinity.fasta > /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_trimmomatic.fasta 2> /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/error_after_trimmomatic_output.log &
        nohup perl /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/get_longest_isoform_seq_per_trinity_gene.pl /media/raid/Wee/WeeYeZhi/output/cd_hit_est_results/cd_hit_est_after_trimmomatic/cd_hit_est_after_trimmomatic.fasta > /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta 2> /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/error_after_cdhitest_trinity_output.log &
        """, language="bash")

        st.write("###")

        st.write("**10. Remove adapter and vector contamination via NCBI FCS-Adaptor & remove foreign biological contamination via NCBI-FCS-GX from the already deduplicated CPB transcriptome assembly**")
        st.write("✔️create a conda environment, named ncbi-fcs-gx, activate the environment & install ncbi-fcs-gx")
        st.code("""
        conda create -n ncbi-fcs-gx
        conda activate ncbi-fcs-gx
        conda install bioconda::ncbi-fcs-gx
        conda install conda-forge::curl # install curl to download the scripts later
        """, language="bash")
        st.write("✔️remove adapter & vector contamination via NCBI FCS-Adaptor")
        st.write("-screen the transcriptome assembly to obtain a list of detected adapter and vector sequences")
        st.code("""
        curl -LO https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh # download the run_fcsadaptor.sh script
        chmod 755 run_fcsadaptor.sh # change the permission of the downloaded script so that you can read, write and execute the script
        nohup bash run_fcsadaptor.sh --fasta-input /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta --output-dir /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/run_fcsadaptor_results --euk > /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/run_fcsadaptor_results/run_fcsadaptor_output.log 2>&1 & # run the run_fcsadaptor.sh script to screen the transcriptome assembly
        """, language="bash")
        st.write("-clean the transcriptome assembly using the list of the detected adapter and vector sequences")
        st.code("""
        curl -LO https://github.com/ncbi/fcs/raw/main/dist/fcs.py # download the fcs.py python script 
        cat /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta | nohup python3 ./fcs.py clean genome --action-report /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/run_fcsadaptor_results/fcs_adaptor_report.txt --output /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/fcspy_results/clean.fasta --contam-fasta-out /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/fcspy_results/contam.fasta > fcspy_output.log 2>&1 & # run the fcs.py script to clean the transcriptome
        """, language="bash")
        st.write("✔️remove foreign biological contamination via NCBI FCS-GX")
        st.write("-download the NCBI FCS-GX database")
        st.code("""
        conda activate curl && curl -LO https://github.com/ncbi/fcs/raw/main/dist/fcs.py # activate the curl conda environment first before downloading the fcs.py runner script
        SOURCE_DB_MANIFEST="https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest" # set the environment variable for SOURCE_DB_MANIFEST 
        LOCAL_DB="/media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/database" # set the environment variable for LOCAL_DB
        python3 fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/gxdb" # download the NCBI FCS-GX database
        ls -lh "$LOCAL_DB/gxdb" # verify the download
        """, language="bash")
        st.write("-screen the transcriptome assembly to obtain a list of detected foreign biological contaminants")
        st.code("""
        GXDB_LOC="/media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/database" # set environment variable for GXDB_LOC
        nohup python3 ./fcs.py screen genome --fasta /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/fcspy_results/clean.fasta --out-dir /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/screen_genome_results --gx-db "$GXDB_LOC/gxdb" --tax-id 538958 > /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/screen_genome_results/screen_genome_output.log 2>&1 & # screen the transcriptome assembly to detect foreign biological contaminant sequences
        """, language="bash")
        st.write("-clean the transcriptome assembly to remove foreign biological contaminants")
        st.code("""
        cat /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_adaptor_results/fcspy_results/clean.fasta | nohup python3 ./fcs.py clean genome --action-report /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/screen_genome_results/clean.538958.fcs_gx_report.txt --output /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/clean_genome_results/clean.fasta --contam-fasta-out /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/clean_genome_results/contam.fasta > /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/clean_genome_results/clean_genome_output.log 2>&1 & # clean the transcriptome assembly to remove foreign biological contaminants
        """, language="bash")
        st.write("✔️remove the tilde suffix from the assembly fasta file after running NCBI FCS-adaptor and NCBI FCS-GX")
        st.code("""
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/clean_genome_results/clean.fasta > /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/ncbi_fcs_gx_results/clean_genome_results/modified_clean.fasta
        """, language="bash")

        st.write("###")

        st.write("**11. Evaluate the BUSCO completeness score of the CPB transcriptome assembly for 4 times (before deduplication, after deduplication, after NCBI FCS & after filtering bad contigs via Transrate)**")
        st.write("✔️create 'busco' conda environment, activate the environment and install busco via conda")
        st.code("""
        conda create -n busco
        conda activate busco
        conda install bioconda::busco
        """, language="bash")
        st.write("✔️verify the installation of busco")
        st.code("""
        which busco
        """, language="bash")
        st.write("✔️determine the lineage file suitable to be used for CPB genome by listing the lineage datasets available in BUSCO first, followed by referring to the NCBI BioProject of CPB (taxonomy) to check its taxonomy")
        st.code("busco --list-datasets", language="bash")
        st.code("busco --list-datasets | grep -i lepidoptera_odb12", language="bash")
        st.write("✔️run BUSCO to evaluate the completeness of the CPB transcriptome assembly")
        st.code("""
        nohup busco -m transcriptome -i /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM/trinity_assembly_done_by_LKM.Trinity.fasta -c 16 -l lepidoptera_odb12 -o /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_results_of_transcriptome_assembly_done_by_LKM > /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_results_of_transcriptome_assembly_done_by_LKM/CPB_transcriptome_assembly_busco_done_by_LKM_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly assembled using reads sequenced by LKM
        nohup busco -m transcriptome -i /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS/trinity_assembly_done_by_INBIOSIS.Trinity.fasta -c 16 -l lepidoptera_odb12 -o /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_results_of_transcriptome_assembly_done_by_INBIOSIS > /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_results_of_transcriptome_assembly_done_by_INBIOSIS/busco_trinity_assembly_done_by_INBIOSIS_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly assembled using reads sequenced by INBIOSIS
        nohup busco -m transcriptome -i /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_assembly_done_by_LKM_and_INBIOSIS.Trinity.fasta -c 16 -l lepidoptera_odb12 -o /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_results_of_transcriptome_assembly_done_by_LKM_and_INBIOSIS > /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_results_of_transcriptome_assembly_done_by_LKM_and_INBIOSIS/CPB_transcriptome_assembly_busco_done_by_LKM_and_INBIOSIS_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly assembled using reads sequenced by both LKM & INBIOSIS
        nohup busco -m transcriptome -i /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta -c 16 -l lepidoptera_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after deduplication
        nohup busco -m transcriptome -i /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta -c 16 -l lepidoptera_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after NCBI FCS
        nohup busco -m transcriptome -i /media/raid/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -c 16 -l lepidoptera_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after removing bad contigs via TransRate
        """, language="bash")
        st.write("✔️draw BUSCO plot for comparison using python script (built inside the busco conda package))")
        st.code("""
        which busco # to find the location of BUSCO executable 
        head -n 5 /home/cbr15/anaconda3/envs/busco/bin/generate_plot.py # display the first 5 lines to see whether there is a shebang line, "#!/usr/bin/env python3" at the first line of the file. If have, then no need to specify the flag, "python3" while running the command to draw the BUSCO plot
        nohup python3 /home/cbr15/anaconda3/envs/busco/bin/generate_plot.py -wd /media/raid/Wee/WeeYeZhi/output/Buscoresults/BUSCO_summaries_transcriptome_lepidoptera_odb12 -rt specific --no_r > busco_plot_output.log 2>&1 & # output only the R code to be run inside the RStudio later to draw the BUSCO plot # need to specify the flag, "python3" if the file doesn't have the shebang line
        """, language="bash")
        # ----LOAD THE GENERATE_BUSCO_PLOT PYTHON SCRIPT----
        # Check if the file exists before reading
        if generate_busco_plot_python_script.exists():
            with open(generate_busco_plot_python_script, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download generate_plot.py",
                data=script_byte,
                file_name=generate_busco_plot_python_script.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{generate_busco_plot_python_script.name} does not exist.")
        st.write(
            "✔️generate the BUSCO plot within RStudio by using the generated busco_figure.R script using your own laptop")
        # ----LOAD THE GENERATE_BUSCO_PLOT RSCRIPT----
        # Check if the file exists before reading
        if generate_busco_plot_Rscript.exists():
            with open(generate_busco_plot_Rscript, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download busco_figure.R",
                data=script_byte,
                file_name=generate_busco_plot_Rscript.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{generate_busco_plot_Rscript.name} does not exist.")

        st.write("###")

        st.write("**12. Compute the basic statistics of the CPB transcriptome assembly via Trinity**")
        st.write("✔️activate the 'trinity' conda environment")
        st.code("""
        conda activate trinity
        """, language="bash")
        st.write("✔️verify the location and existence of TrinityStats.pl script")
        st.code("""
        which TrinityStats.pl
        """, language="bash")
        st.write("✔️run the TrinityStats.pl script")
        st.code("""
        perl /home/cbr15/anaconda3/envs/trinity/bin/TrinityStats.pl /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_assembly_done_by_LKM_and_INBIOSIS.Trinity.fasta > stats_trinity_assembly_after_trimmomatic.log
        perl /home/cbr15/anaconda3/envs/trinity/bin/TrinityStats.pl /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta > stats_trinity_assembly_after_cdhitest_trinity.log
        perl /home/cbr15/anaconda3/envs/trinity/bin/TrinityStats.pl /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta > stats_trinity_assembly_after_cdhitest_trinity_NCBI_FCS.log
        perl /home/cbr15/anaconda3/envs/trinity/bin/TrinityStats.pl /media/raid/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta > stats_trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.log
        """, language="bash")

        st.write("###")

        st.write("**13. Perform reference-free quality assessment of de novo CPB transcriptome assembly for 4 times via TransRate (before deduplication, after deduplication, & after NCBI-FCS & after removing bad contigs via TransRate)**")
        st.write("✔️install transrate via docker container")
        st.code("""
        docker pull genevia/transrate:v1.0.3_orp # pull the transrate docker image from the DockerHub registry
        docker images # verify that the transrate docker image has been successfully pulled
        """, language="bash")
        st.write("✔️run transrate to evaluate the CPB transcriptome assembly (before deduplication) via docker")
        st.code("""
        docker run --user 1000:1000 -d --name orp_transrate -v /media/raid/Wee/WeeYeZhi/output:/data genevia/transrate:v1.0.3_orp /bin/bash -c "transrate --assembly=/data/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_assembly_done_by_LKM_and_INBIOSIS.Trinity.fasta --left=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz --right=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz --threads 16 --output=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_trimmomatic" 
        """, language="bash")
        st.write("✔️run transrate to evaluate the CPB transcriptome assembly (after deduplication) via docker")
        st.code("""
        docker run --user 1000:1000 -d --name orp_transrate_after_deduplication -v /media/raid/Wee/WeeYeZhi/output:/data genevia/transrate:v1.0.3_orp /bin/bash -c "transrate --assembly=/data/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta --left=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz --right=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz --threads 16 --output=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity" 
        """, language="bash")
        st.write("✔️run transrate to evaluate the CPB transcriptome assembly (after NCBI-FCS) via Docker")
        st.code("""
        docker run --user 1000:1000 -d --name orp_transrate_after_NCBI_FCS -v /media/raid/Wee/WeeYeZhi/output:/data genevia/transrate:v1.0.3_orp /bin/bash -c "transrate --assembly=/data/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta --left=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz --right=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz --threads 16 --output=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS" 
        """, language="bash")
        st.write("✔️run transrate to evaluate the CPB transcriptome assembly (after filtering & removing bad contigs via Transrate) via Docker")
        st.code("""
        docker run --user 1000:1000 -d --name orp_transrate_after_NCBI_FCS_transrate_filtering -v /media/raid/Wee/WeeYeZhi/output:/data genevia/transrate:v1.0.3_orp /bin/bash -c "transrate --assembly=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --left=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz --right=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz --threads 16 --output=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering" 
        """, language="bash")
        st.write("✔️show the log file content of the transrate docker that is running in the background")
        st.code("""
        docker logs orp_transrate # show the log file content of the run by using the designated name of the docker container before deduplication
        docker logs orp_transrate_after_deduplication # show the log file content of the run by using the designated name of the docker container after deduplication 
        docker logs orp_transrate_after_NCBI_FCS # show the log file content of the run by using the designated name of the docker container after deduplication & NCBI-FCS
        docker logs orp_transrate_after_NCBI_FCS_transrate_filtering # show the log file content of the run by using the designated name of the docker container after deduplication & NCBI-FCS & removing bad contigs via Transrate      
        docker ps -a # list all the current actively running and stopped docker containers. Those that exited with code 0 is the successful run
        """, language="bash")
        st.write("✔️remove the tilde suffix (~) inside the transcriptome assembly fasta headers before submission into the TSA database")
        st.code("""
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' good.clean.fasta > modified.good.clean.fasta
        """, language="bash")
        st.markdown("[Visit TransRate GitHub Page](https://github.com/blahah/transrate)")
        st.markdown("[Visit TransRate Conda Installation Page](https://anaconda.org/bioconda/transrate)")
        st.markdown("[Visit TransRate User Manual Page](https://hibberdlab.com/transrate/)")
        st.markdown("[Visit TransRate Step by Step User Manual Page](https://hibberdlab.com/transrate/getting_started.html)")

        st.write("###")

        st.write("**14. Run Salmon to perform transcript abundance estimation**")
        st.write("✔️create a 'salmon' conda environment, activate the environment and install salmon via bioconda")
        st.code("""
        conda create -n salmon
        conda activate salmon
        conda install bioconda::salmon
        """, language="bash")
        st.write("✔️verify the installation of Salmon")
        st.code("""
        which salmon
        salmon index -h
        salmon quant -h
        """, language="bash")
        st.write("✔️build a salmon index for the CPB transcriptome assembly in mapping-based mode")
        st.code("""
        nohup salmon index -t /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_trimmomatic.fasta -p 16 -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/transcripts_index_after_trimmomatic -k 31 > /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_trimmomatic/salmon_index_after_trimmomatic_output.log 2>&1 &
        nohup salmon index -t /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta -p 16 -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_cdhitest_trinity/transcripts_index_after_cdhitest_trinity -k 31 > /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_cdhitest_trinity/salmon_index_after_cdhitest_trinity_output.log 2>&1 &
        nohup salmon index -t /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta -p 16 -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_cdhitest_trinity_NCBI_FCS/transcripts_index_after_cdhitest_trinity_NCBI_FCS -k 31 > /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_cdhitest_trinity_NCBI_FCS/salmon_index_after_cdhitest_trinity_NCBI_FCS_output.log 2>&1 &
        nohup salmon index -t /media/raid/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -p 16 -i /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/transcripts_index_after_cdhitest_trinity_NCBI_FCS_transrate_filtering -k 31 > /media/raid/Wee/WeeYeZhi/output/salmon_results/salmon_index/salmon_index_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/salmon_index_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 &
        """, language="bash")
        st.write("✔️perform transcript abundance quantification against the transcriptome assembly for three times (before running deduplication, after running deduplication, and after running NCBI-FCS)")
        st.code("""
        dos2unix salmon_quantify_after_trimmomatic.sh
        chmod +x salmon_quantify_after_trimmomatic.sh
        nohup bash salmon_quantify_after_trimmomatic.sh > salmon_quantify_after_trimmomatic_output.log 2>&1 &
        """, language="bash")
        # ----LOAD SALMON TRANSCRIPT QUANTIFICATION BASH AFTER TRIMMOMATIC FILE----
        # Check if the file exists before reading
        if salmon_transcript_quantification_after_trimmomatic.exists():
            with open(salmon_transcript_quantification_after_trimmomatic, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification After Trimmomatic Bash File",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_trimmomatic.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_trimmomatic.name} does not exist.")
        st.code("""
        dos2unix salmon_quantify_after_cdhitest_trinity.sh
        chmod +x salmon_quantify_after_cdhitest_trinity.sh
        nohup bash salmon_quantify_after_cdhitest_trinity.sh > salmon_quantify_after_cdhitest_trinity_output.log 2>&1 &
        """, language="bash")
        # ----LOAD SALMON TRANSCRIPT QUANTIFICATION BASH AFTER CDHITEST & TRINITY FILE----
        # Check if the file exists before reading
        if salmon_transcript_quantification_after_cdhitest_trinity.exists():
            with open(salmon_transcript_quantification_after_cdhitest_trinity, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification After CDHITEST & Trinity Bash File",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_cdhitest_trinity.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_cdhitest_trinity.name} does not exist.")
        st.code("""
        dos2unix salmon_quantify_after_cdhitest_trinity_NCBI_FCS.sh
        chmod +x salmon_quantify_after_cdhitest_trinity_NCBI_FCS.sh
        nohup bash salmon_quantify_after_cdhitest_trinity_NCBI_FCS.sh > salmon_quantify_after_cdhitest_trinity_NCBI_FCS_output.log 2>&1 &
        """, language="bash")
        # ----LOAD SALMON TRANSCRIPT QUANTIFICATION BASH AFTER NCBI FCS FILE----
        # Check if the file exists before reading
        if salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS.exists():
            with open(salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification After NCBI FCS Bash File",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS.name} does not exist.")
        st.code("""
        dos2unix salmon_quantify_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.sh
        chmod +x salmon_quantify_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.sh
        nohup bash salmon_quantify_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.sh > salmon_quantify_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 &
        """, language="bash")
        # ----LOAD SALMON TRANSCRIPT QUANTIFICATION BASH AFTER NCBI FCS & TRANSRATE FILTERING FILE----
        # Check if the file exists before reading
        if salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.exists():
            with open(salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification After NCBI FCS & Transrate Filtering Bash File",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.name} does not exist.")
        st.write("✔️remove the tilde suffix inside all the quant.sf files produced by the salmon quantification run after NCBI FCS (otherwise you are going to encounter an error complaining that it cannot find certain transcripts ID due to mismatches between transcripts ID inside the both the gene_to_trans_map file and salmon quantification file)")
        st.code("""
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR11266554.sf # remember to copy the source code from the pycharm file and paste it to the terminal. You cannot copy and paste the command from the logbook website.
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR11266555.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR11266556.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR9690969.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR9690970.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR9690971.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR9690972.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR9690973.sf
        sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' quant.sf > quant_SRR9690974.sf
        """, language="bash")
        st.markdown("[Visit Salmon GitHub Page](https://github.com/COMBINE-lab/salmon)")
        st.markdown("[Visit Salmon User Manual Page](https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon)")
        st.markdown("[Visit Salmon Conda Installation Page](https://anaconda.org/bioconda/salmon)")

        st.write("###")

        st.write("**15. Compute the ExN50 statistics of the transcriptome assemblies (before deduplication, after deduplication, after NCBI FCS and after removing bad contigs via TransRate)**")
        st.write("✔️activate the 'trinity' conda environment")
        st.code("""
        conda activate trinity
        """, language="bash")
        st.write("✔️generate the gene to trans map file for the CPB transcriptome assembly for 4 times (before deduplication, after deduplication & after NCBI FCS & after removing bad contigs via TransRate)")
        st.code("""
        perl /home/cbr15/anaconda3/envs/trinity/bin/get_Trinity_gene_to_trans_map.pl /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_trimmomatic.fasta > /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_trimmomatic/trinity_assembly_after_trimmomatic.Trinity.fasta.gene_trans_map
        perl /home/cbr15/anaconda3/envs/trinity/bin/get_Trinity_gene_to_trans_map.pl /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta > /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity/trinity_assembly_after_cdhitest_trinity.Trinity.fasta.gene_trans_map
        perl /home/cbr15/anaconda3/envs/trinity/bin/get_Trinity_gene_to_trans_map.pl /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta > /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/trinity_assembly_after_cdhitest_trinity_NCBI_FCS.Trinity.fasta.gene_trans_map
        perl /home/cbr15/anaconda3/envs/trinity/bin/get_Trinity_gene_to_trans_map.pl /media/raid/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta > /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.Trinity.fasta.gene_trans_map
        """, language="bash")
        st.write("✔️compute and obtain the expression matrix after quantifying the abundance of transcripts via Salmon")
        st.code("""
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_trimmomatic/trinity_assembly_after_trimmomatic.Trinity.fasta.gene_trans_map --out_prefix salmon --name_sample_by_basedir --quant_files /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_trimmomatic/salmon_transcript_quantification_files_after_trimmomatic.list > trinity_abundance_to_estimates_after_trimmomatic_output.log 2>&1 &
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity/trinity_assembly_after_cdhitest_trinity.Trinity.fasta.gene_trans_map --out_prefix salmon --name_sample_by_basedir --quant_files /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity/salmon_transcript_quantification_files_after_cdhitest_trinity.list > trinity_abundance_to_estimates_after_cdhitest_trinity_output.log 2>&1 &
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/trinity_assembly_after_cdhitest_trinity_NCBI_FCS.Trinity.fasta.gene_trans_map --out_prefix salmon --name_sample_by_basedir --quant_files /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/salmon_transcript_quantification_files_after_cdhitest_trinity_NCBI_FCS.list > trinity_abundance_to_estimates_after_cdhitest_trinity_NCBI_FCS_output.log 2>&1 &
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.Trinity.fasta.gene_trans_map --out_prefix salmon --name_sample_by_basedir --quant_files /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/salmon_transcript_quantification_files_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.list > trinity_abundance_to_estimates_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 &
        """, language="bash")
        # ----LOAD SALMON TRANSCRIPT QUANTIFICATION LIST FILE BEFORE DEDUPLICATION----
        # Check if the file exists before reading
        if salmon_transcript_quantification_after_trimmomatic_list.exists():
            with open(salmon_transcript_quantification_after_trimmomatic_list, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification List File Before Deduplication",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_trimmomatic_list.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_trimmomatic_list.name} does not exist.")
            # ----LOAD SALMON TRANSCRIPT QUANTIFICATION LIST FILE AFTER DEDUPLICATION----
            # Check if the file exists before reading
        if salmon_transcript_quantification_after_cdhitest_trinity_list.exists():
            with open(salmon_transcript_quantification_after_cdhitest_trinity_list, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification List File After Deduplication",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_cdhitest_trinity_list.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_cdhitest_trinity_list.name} does not exist.")
            # ----LOAD SALMON TRANSCRIPT QUANTIFICATION LIST FILE AFTER NCBI FCS----
            # Check if the file exists before reading
        if salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_list.exists():
            with open(salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_list, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification List File After NCBI FCS",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_list.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_list.name} does not exist.")
            # ----LOAD SALMON TRANSCRIPT QUANTIFICATION LIST FILE AFTER REMOVING BAD CONTIGS VIA TRANSRATE----
            # Check if the file exists before reading
        if salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_list.exists():
            with open(salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_list, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Salmon Transcript Quantification List File After Removing Bad Contigs via TransRate",
                data=script_byte,
                file_name=salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_list.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{salmon_transcript_quantification_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_list.name} does not exist.")
        st.write("✔️compute the ExN50 statistics of the CPB transcriptome assembly at gene level (instead of transcript level) to remove bias")
        st.code("""
        perl /home/cbr15/anaconda3/envs/trinity/bin/contig_ExN50_statistic.pl /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_trimmomatic/salmon.isoform.TMM.EXPR.matrix /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_trimmomatic.fasta gene | tee ExN50_after_trimmomatic.gene.stats
        perl /home/cbr15/anaconda3/envs/trinity/bin/contig_ExN50_statistic.pl /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity/salmon.isoform.TMM.EXPR.matrix /media/raid/Wee/WeeYeZhi/output/get_longest_isoform_per_Trinity_gene_results/longest_isoform_per_gene_after_cdhitest_trinity.fasta gene | tee ExN50_after_cdhitest_trinity.gene.stats
        perl /home/cbr15/anaconda3/envs/trinity/bin/contig_ExN50_statistic.pl /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/salmon.isoform.TMM.EXPR.matrix /media/raid/Wee/WeeYeZhi/output/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta gene | tee ExN50_after_cdhitest_trinity_NCBI_FCS.gene.stats # run the command, "sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' clean.fasta > modified_clean.fasta" first to remove the ~ suffixes within the assembly file generated and cleaned via NCBI FCS pipeline, before using it to compute the ExN50 values
        perl /home/cbr15/anaconda3/envs/trinity/bin/contig_ExN50_statistic.pl /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/salmon.isoform.TMM.EXPR.matrix /media/raid/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta gene | tee ExN50_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.gene.stats # run the command, "sed -E 's/([^[:space:]]*)~[^[:space:]]*/\1/' clean.fasta > modified_clean.fasta" first to remove the ~ suffixes within the assembly file generated and cleaned via NCBI FCS pipeline, before using it to compute the ExN50 values
        """, language="bash")
        st.write("✔️remove the first line 'nohup: ignoring input' from each of the generated ExN50.stats file before plotting the ExN50 graph via sed to avoid encountering error when parsing the file")
        st.code("""
        sed -i '1d' ExN50_after_trimmomatic.gene.stats
        sed -i '1d' ExN50_after_cdhitest_trinity.gene.stats
        sed -i '1d' ExN50_after_cdhitest_trinity_NCBI_FCS.gene.stats
        sed -i '1d' ExN50_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.gene.stats
        """, language="bash")
        st.write("✔️install the ggplot2 conda package within the trinity conda environment via conda-forge")
        st.code("""
        conda install conda-forge::r-ggplot2
        """, language="bash")
        st.write("✔️download the 'plot_ExN50_statistic.Rscript' file from https://github.com/trinityrnaseq/trinityrnaseq/blob/master/util/misc/plot_ExN50_statistic.Rscript as the trinity conda package does not contain this Rscript (need to download separately)")
        st.write("✔️plot and visualize the ExN50 statistics results of the CPB transcriptome assembly")
        st.code("""
        nohup Rscript /media/raid/Wee/WeeYeZhi/output/trinity_ExN50_gene_level_statistics_calculation_results/trinity_assembly_after_trimmomatic/plot_ExN50_statistic.Rscript ExN50_after_trimmomatic.gene.stats > plot_ExN50_statistic_after_trimmomatic_output.log 2>&1 &
        nohup Rscript /media/raid/Wee/WeeYeZhi/output/trinity_ExN50_gene_level_statistics_calculation_results/trinity_assembly_after_cdhitest_trinity/plot_ExN50_statistic.Rscript ExN50_after_cdhitest_trinity.gene.stats > plot_ExN50_statistic_after_cdhitest_trinity_output.log 2>&1 &
        nohup Rscript /media/raid/Wee/WeeYeZhi/output/trinity_ExN50_gene_level_statistics_calculation_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS/plot_ExN50_statistic.Rscript ExN50_after_cdhitest_trinity_NCBI_FCS.gene.stats > plot_ExN50_statistic_after_cdhitest_trinity_NCBI_FCS_output.log 2>&1 &
        nohup Rscript /media/raid/Wee/WeeYeZhi/output/trinity_ExN50_gene_level_statistics_calculation_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/plot_ExN50_statistic.Rscript ExN50_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.gene.stats > plot_ExN50_statistic_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 &
        """, language="bash")

# Phase 2: Reference-Based Transcriptomics Analysis

if selected == "Phase 2: Transcriptomic and Structural-Based Analysis":
    with st.container():
        st.write("---")
        st.header("Phase 2: Transcriptomic and Structural-Based Analysis 🪢")
        st.write("###")

        st.write("**1. Identify coding regions within transcript sequences/transcriptome via TransDecoder**")
        st.write("✔️pull the TransDecoder docker image from DockerHub")
        st.code("""
        docker pull trinityrnaseq/transdecoder # you are actually pulling the 'latest' docker image & please record the digest ID 
        """, language="bash")
        st.write("✔️obtain the digest ID of the latest Trinotate docker image")
        st.code("""
        docker image inspect --format '{{.RepoDigests}}' trinityrnaseq/transdecoder:latest # you will obtain trinityrnaseq/transdecoder@sha256:80fa372fe94bc0ac4440fea7905c47e05187d13b1e837f8dca72c0ca49ab1bb4 where other scientists can run the command, docker pull trinityrnaseq/transdecoder@sha256:80fa372fe94bc0ac4440fea7905c47e05187d13b1e837f8dca72c0ca49ab1bb4 to use the same docker transdecoder image with the same software environment
        """, language="bash")
        st.write("✔️check whether the TransDecoder docker image has been successfully pulled")
        st.code("""
        docker images
        """, language="bash")
        st.write("✔️verify the pull of the TransDecoder docker image")
        st.code("""
        docker run --rm trinityrnaseq/transdecoder which makeblastdb
        docker run --rm trinityrnaseq/transdecoder which hmmpress
        docker run --rm trinityrnaseq/transdecoder which TransDecoder.LongOrfs
        docker run --rm trinityrnaseq/transdecoder which blastp
        docker run --rm trinityrnaseq/transdecoder which hmmsearch
        docker run --rm trinityrnaseq/transdecoder which TransDecoder.Predict
        """, language="bash")
        st.write("✔️download UniProt Swiss-Prot FASTA file externally from UniProt’s website (https://www.uniprot.org/downloads) under the Swiss-Prot section & create the UniProt database")
        st.code("""
        gunzip uniprot_sprot.fasta.gz # uncompress the file first before creating the database as makeblastdb does not accept gzipped file as input
        docker run --user 1000:1000 -d --name make_swissprot_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "makeblastdb -in /data/transdecoder_results/make_swissprot_database/uniprot_sprot.fasta -dbtype prot -out /data/transdecoder_results/make_swissprot_database/uniprot_sprot"
        """, language="bash")
        # ----LOAD UNIPROT SWISS-PROT FASTA FILE----
        # Check if the file exists before reading
        if swissprot_fasta_file.exists():
            with open(swissprot_fasta_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download UniProt Swiss-Prot Gzipped Fasta File",
                data=script_byte,
                file_name=swissprot_fasta_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{swissprot_fasta_file.name} does not exist.")
        st.write("✔️download and index the HMM profile database to run hmmsearch")
        st.code("""
        wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz # download the main HMM profile database for hmmsearch from the official EBI FTP server
        gunzip Pfam-A.hmm.gz # uncompress the file
        docker run --user 1000:1000 -d --name hmmpress -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "hmmpress /data/transdecoder_results/hmmsearch_results/Pfam-A.hmm"
        """, language="bash")
        st.write("✔️extract the long open reading frames (ORFs) from the CPB transcriptome assembly")
        st.code("""
        docker run --user 1000:1000 -d --name ORF_extraction -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "TransDecoder.LongOrfs -S -t /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --output_dir /data/transdecoder_results/transdecoder_ORF_extraction_results"
        """, language="bash")
        st.write("✔️run BLASTp to blast the protein sequences predicted from the CPB transcriptome assembly against the SwissProt protein database")
        st.code("""
        docker run --user 1000:1000 -d --name blastp_search -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "blastp -query /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep -db /data/transdecoder_results/make_swissprot_database/uniprot_sprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > /data/transdecoder_results/ncbi_blastp_results/blastp.outfmt6"
        """, language="bash")
        st.write("✔️run HMMSEARCH to identify the protein domains in the predicted ORFs using the Pfam database")
        st.code("""
        docker run --user 1000:1000 -d --name hmmsearch -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "hmmsearch --cpu 16 -E 1e-10 --domtblout /data/transdecoder_results/hmmsearch_results/pfam.domtblout /data/transdecoder_results/hmmsearch_results/Pfam-A.hmm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep"
        """, language="bash")
        st.write("✔️filter and retain predicted ORFs that have supported BLASTp hits and Pfam domains")
        st.code("""
        docker run --user 1000:1000 -d --name ORF_prediction -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "TransDecoder.Predict -t /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --retain_pfam_hits /data/transdecoder_results/hmmsearch_results/pfam.domtblout --retain_blastp_hits /data/transdecoder_results/ncbi_blastp_results/blastp.outfmt6 -O /data/transdecoder_results/transdecoder_ORF_extraction_results"
        """, language="bash")
        st.markdown("[Visit TransDecoder Docker Installation Page](https://github.com/TransDecoder/TransDecoder/wiki/TransDecoder_Docker_or_Singularity)")
        st.markdown("[Visit TransDecoder Docker Image on DockerHub](https://hub.docker.com/r/trinityrnaseq/transdecoder)")
        st.markdown("[Visit TransDecoder GitHub User Manual Page](https://github.com/TransDecoder/TransDecoder/wiki)")
        st.markdown("[Download the Pfam database](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/)")
        st.markdown("[Download UniProt Swiss-Prot Fasta File to create SwissProt Database](https://www.uniprot.org/downloads)")
        st.markdown("[Read how to perform BLASTp](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/)")

        st.write("###")

        st.write("**2. Annotate the CPB transcriptome assembly via Trinotate**")
        st.write("✔️pull the Trinotate docker image from DockerHub")
        st.code("""
        docker pull trinityrnaseq/trinotate
        """, language="bash")
        st.write("✔️check whether the Trinotate docker image has been successfully pulled")
        st.code("""
        docker images
        """, language="bash")
        st.write("✔️print the full path to the Trinotate executable file and verify the pull of the Trinotate docker image")
        st.code("""
        docker run --rm trinityrnaseq/trinotate find / -type f -name Trinotate 2>/dev/null
        docker run --rm trinityrnaseq/trinotate /usr/local/src/Trinotate/Trinotate -h
        """, language="bash")
        st.write("✔️start running the Trinotate docker container (with the container being named as 'trinotate') and keep it alive in the background, and set the environment variables")
        st.code("""
        docker run -dit --user 1000:1000 --name trinotate -e TRINOTATE_HOME=/usr/local/src/Trinotate/Trinotate -e TRINOTATE_DATA_DIR=/data/trinotate_results/trinotate_data_dir -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash # remember to run docker stop and docker rm after running trinotate, otherwise it is going to consume resources
        """, language="bash")
        st.write("✔️build the trinotate sqlite database")
        st.code("""
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate --create --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --trinotate_data_dir /data/trinotate_results/trinotate_data_dir --use_diamond" # bear in mind that the build process will help you to automatically download the two resources, namely 'uniprot_sprot.pep' and 'Pfam-A.hmm.gz'.
        """, language="bash")
        st.write("✔️prepare both the BLASTp and Pfam databases that have been previously created by the built process")
        st.code("""
        docker exec -d trinotate /bin/bash -c "makeblastdb -in \$TRINOTATE_DATA_DIR/uniprot_sprot.pep -dbtype prot"
        docker exec -d trinotate /bin/bash -c "gunzip \$TRINOTATE_DATA_DIR/Pfam-A.hmm.gz"
        docker exec -d trinotate /bin/bash -c "hmmpress \$TRINOTATE_DATA_DIR/Pfam-A.hmm"
        """, language="bash")
        st.write("✔️initialize the trinotate sqlite database")
        st.code("""
        docker exec -d trinotate /bin/bash -c "/$TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --init --gene_trans_map /data/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.Trinity.fasta.gene_trans_map --transcript_fasta /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --transdecoder_pep /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep"
        """, language="bash")
        st.write("✔️run diamond in BLASTx mode to blast translated nucleotide against protein")
        st.code("""
        docker exec -d trinotate /bin/bash -c "diamond blastx -d \$TRINOTATE_DATA_DIR/uniprot_sprot.pep -q /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -o /data/trinotate_results/diamond_blastx_results/blastx.outfmt6 -f 6 -p 16 -e 1e-5"
        """, language="bash")
        st.write("✔️run BLASTp & hmmsearch again using the files generated by trinotate sqlite database construction step")
        st.code("""
        docker exec -d trinotate /bin/bash -c "blastp -query /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep -db /data/trinotate_results/trinotate_data_dir/uniprot_sprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > /data/trinotate_results/ncbi_blastp_results/blastp.outfmt6"
        docker exec -d trinotate /bin/bash -c "hmmsearch --cpu 16 -E 1e-10 --domtblout /data/trinotate_results/hmmsearch_results/pfam.domtblout /data/trinotate_results/trinotate_data_dir/Pfam-A.hmm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep"
        """, language="bash")
        st.write("✔️create an isolated conda environment called signalp6, activate the conda environment, install the gzipped file of SignalP6 at https://services.healthtech.dtu.dk/services/SignalP-6.0/, verify the installation and run SignalP6 to predict the signal peptide molecules of the predicted CPB protein sequences")
        st.code("""
        conda create -n signalp6 python=3.8 -y
        conda activate signalp6
        tar -xzvf signalp-6.0h.fast.tar.gz # extract the signalp6 file
        cd /media/raid/Wee/WeeYeZhi/output/trinotate_results/signalp6_installation/signalp6_fast
        pip install torch==1.8.1+cpu torchvision==0.9.1+cpu torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
        pip install matplotlib tqdm
        pip install signalp-6-package/
        SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))") && mkdir -p "$SIGNALP_DIR/model_weights" && cp -r /media/raid/Wee/WeeYeZhi/output/trinotate_results/signalp6_installation/signalp6_fast/signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"
        which signalp6
        nohup signalp6 --fastafile /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep --organism eukarya --mode fast --format txt --output_dir /media/raid/Wee/WeeYeZhi/output/trinotate_results/signalp6_results > /media/raid/Wee/WeeYeZhi/output/trinotate_results/signalp6_results/signalp6_output.log 2>&1 &
        """, language="bash")
        st.write("✔️install the gzipped file of TMHMM at https://services.healthtech.dtu.dk/services/TMHMM-2.0/, verify the installation and run TMHMM to predict the transmembrane domains of the predicted CPB protein sequences")
        st.code("""
        tar -xzvf tmhmm-2.0c.Linux.tar.gz # extract the tmhmm file
        ./tmhmm -h # navigate into the working directory, /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/bin before you run the command
        tmhmm -m TMHMM2.0.model -f /media/raid/Wee/WeeYeZhi/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep
        """, language="bash")
        st.write("✔️install DeepTHMMM via Docker & run it to predict the transmembrane domains of the predicted CPB protein sequences")
        st.code("""
        docker pull gnick18/deepthmmm # pull the DeepTHMMM docker image from DockerHub
        docker images # verify whether the DeepTHMMM docker image has been successfully pulled
        docker run -d --user 1000:1000 --name deepthmmm -v /media/raid/Wee/WeeYeZhi/output:/data gnick18/deepthmmm /bin/bash -c "biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep > /data/trinotate_results/deeptmhmm_results/tmhmm.out" 
        docker ps # check whether the DeepTHMMM docker container has been successfully created and running
        """, language="bash")
        st.write("✔️install eggNOG-mapper via docker and run it to assign orthologs/evolutionary functional relationships to the predicted protein sequences of CPB")
        st.code("""
        docker pull nanozoo/eggnog-mapper:2.1.13--c16a7d2 # pull the eggNOGmapper Docker image from DockerHub
        docker images # verify whether the eggNOGmapper Docker image has been successfully pulled
        docker run -d --user 1000:1000 --name download_eggnogmapper_database -v /media/raid/Wee/WeeYeZhi/output:/data nanozoo/eggnog-mapper:2.1.13--c16a7d2 /bin/bash -c "download_eggnog_data.py --data_dir /data/trinotate_results/eggnogmapper_database" # download the eggNOGmapper database
        docker logs download_eggnogmapper_database # view the log content of the run
        docker run -d --user 1000:1000 --name eggnogmapper_search -v /media/raid/Wee/WeeYeZhi/output:/data nanozoo/eggnog-mapper:2.1.13--c16a7d2 /bin/bash -c "emapper.py -i /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep -m diamond --data_dir /data/trinotate_results/eggnogmapper_database --cpu 16 --itype protein -o CPB --output_dir /data/trinotate_results/eggnogmapper_results" # run eggNOGmapper to assign orthologs to the predicted protein sequences of CPB
        docker logs eggnogmapper_search # view the log content of the run
        docker ps # check whether the eggNOGmapper Docker container is still running in the background smoothly, just in case if the container exits due to errors 
        """, language="bash")
        st.write("✔️install infernal via docker and run it")
        st.code("""
        docker pull rfam/rfam-3d-seed-alignments # pull the infernal docker image from DockerHub
        docker images # verify whether the infernal docker image has been successfully pulled
        wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz # download the Rfam covariance model database via ftp (just in case if the docker image does not contain it)
        wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin # download the Rfam clan database via ftp (just in case if the docker image does not contain it)
        gunzip Rfam.cm.gz # unzip the downloaded Rfam CM database
        docker run --rm rfam/rfam-3d-seed-alignments cmpress -h # output the running options available in the cmpress program
        docker run --rm rfam/rfam-3d-seed-alignments cmscan -h # output the running options available in the cmscan program
        docker run -d --user 1000:1000 --name infernal_cmpress -v /media/raid/Wee/WeeYeZhi/output:/data rfam/rfam-3d-seed-alignments /bin/bash -c "cmpress Rfam.cm" # compress and index the downloaded covariance model (CM) database via cmpress. 
        docker run -d --user 1000:1000 --name infernal_cmscan -v /media/raid/Wee/WeeYeZhi/output:/data rfam/rfam-3d-seed-alignments /bin/bash -c "cmscan --cut_ga --rfam --nohmmonly --clanin Rfam.clanin --oskip --fmt 2 -o infernal_output.txt --tblout infernal.out Rfam.cm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.cds" # scan the predicted CPB coding regions of the transcript sequences against the CM models to check whether any of your predicted coding sequences match known RNA families from Rfam. Later, you need to filter and remove all the hits (detected non-coding RNAs) if there are any since you are only interested in the mRNA coding regions. Run infernal as a QC check to filter non-coding RNAs from the coding regions of the transcriptome
        """, language="bash")
        st.write("✔️Load the BLASTp, BLASTx, Pfam, SignalP6, TMHMM, Infernal & eggNOGmapper results into the created Trinotate SQLite Database")
        st.code("""
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_swissprot_blastp /data/transdecoder_results/ncbi_blastp_results/blastp.outfmt6" 
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_swissprot_blastx /data/trinotate_results/diamond_blastx_results/blastx.outfmt6"
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_pfam /data/transdecoder_results/hmmsearch_results/pfam.domtblout"
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_signalp /data/trinotate_results/signalp6_results/signalp6.out"
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_tmhmm /data/trinotate_results/deeptmhmm_results/tmhmm.out"
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_infernal /data/trinotate_results/infernal_results/infernal.out"
        docker exec -d trinotate /bin/bash -c "\$TRINOTATE_HOME/Trinotate /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite LOAD_eggnog /data/trinotate_results/eggnogmapper_results/eggnogmapper.out" # please modify the name of the output file after you have successfully run eggnogmapper
        """, language="bash")
        st.write("Generate the CPB transcriptome annotation report")
        st.code("""
        $TRINOTATE_HOME/Trinotate myTrinotate.sqlite --report > trinotate_annotation_report.tsv
        """, language="bash")
        st.markdown("[Visit Trinotate GitHub Wikipedia](https://github.com/Trinotate/Trinotate/wiki)")
        st.markdown("[Read how to load the results into the Trinotate SQLite Database](https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinotate-Functional-Annotation)")
        st.markdown("[Read how to cite Trinotate and all of the other tools](https://github.com/Trinotate/Trinotate/wiki/Lit-References)")
        st.markdown("[Visit Trinotate GitHub User Manual Page](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)")
        st.markdown("[Read how to run DeepTHMMM locally and via Docker](https://dtu.biolib.com/DeepTMHMM)")
        st.markdown("[Visit TMHMM GitHub User Manual Page](https://github.com/dansondergaard/tmhmm.py)")
        st.markdown("[Visit SignalP6 GitHub User Manual Page](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md)")
        st.markdown("[Visit TMHMM User Manual & Download Page](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)")
        st.markdown("[Visit SignalP6 User Manual & Download Page](https://services.healthtech.dtu.dk/services/SignalP-6.0/)")
        st.markdown("[Read Infernal Documentation](http://eddylab.org/infernal/)")
        st.markdown("[Download the latest covariance model files (.cm) & clan files (.clanin) manually from the Rfam database](https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/)")
        st.markdown("[How to solve the issue of getting an empty output results file/columns](https://github.com/Trinotate/Trinotate.github.io/issues/70)")

        st.write("###")

        st.write("**3. Perform differential expression gene (DEG) analysis using R**")
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
        """, language="r")
        st.write("✔️verify the installation of R packages by checking and recording its package version")
        st.code("""
        packageVersion("DESeq2")
        packageVersion("biomaRt")
        packageVersion("rtracklayer")
        packageVersion("GenomicFeatures")
        packageVersion("tidyverse")
        packageVersion("dplyr")
        packageVersion("stringr")
        packageVersion("ggplot2")
        packageVersion("ggrepel")
        """, language="r")
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
        """, language="r")
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
        st.write("✔️create heatmap plot and perform hierarchical clustering analysis using Morpheus (if you dont want to code)")
        st.markdown("[Visit Morpheus Broad Webpage](https://software.broadinstitute.org/morpheus/)")

        st.write("###")

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

        st.write("---")
        st.write("❗Bear in mind that the p-adjusted values provided by the DEG results will help to verify whether the particular gene is really differentially expressed or it's just a false positive result")
        st.write("❗Decide whether or not you need to accept null hypothesis or reject null hypothesis/accept alternative hypothesis after interpreting the DEG results")
        st.write("❗If there is a significance difference in the gene expression level between the two developmental stages where p-value < 0.05 (at 95% statistical significance level), then null hypothesis is rejected/ alternative hypothesis is accepted")
        st.write("❗If there is no significance difference in the gene expression level between the two developmental stages where p-value > 0.05 (at 95% statistical significance level), then null hypothesis is accepted/ alternative hypothesis is rejected")
        st.write("❗You can compare between the DEG results at p-value=0.05 and p-value=0.01 (most stringent gene filtering) and get the intersected results")
        st.write("❗Run these 3 analysis, 'drawing PCA plot', 'inspecting size factor', 'drawing dispersion plot' to double check the quality of the RNA-seq data")
        st.write("❗PCA is your QC/overview tool; heatmaps and volcanoes are for showing off the real DEGs.")
        st.write("❗PCA is typically performed on the expression data (counts or normalized counts) as shown in the 'analysisObject' dataframe, not on the results of differential gene expression as shown in the 'resOrdered_unique'.")
        st.write("---")

        st.write("**4. Install eggNOG-mapper to perform functional annotation (orthology-based functional annotation) of novel genome sequence of C. cramerella**")
        st.write("✔️create & activate the 'eggnogmapper' conda environment")
        st.code("""
               conda create -n eggnogmapper
               conda activate eggnogmapper
        """, language="bash")
        st.write("✔️install eggnog-mapper")
        st.code("conda install -c bioconda eggnog-mapper", language="bash")
        st.write("✔️install the eggnog database by running the below python script")
        st.code("download_eggnog_data.py", language="bash")
        st.write("✔️verify the installation of eggNOG-mapper")
        st.code("""
        emapper.py --version # emapper.py is the main command to invoke and run eggNOG-mapper
        emapper.py --help
        """, language="bash")
        st.write("✔️run eggNOG-mapper")
        st.code("""
        emapper.py -i braker.proteins.faa -o cramerella_annot --itype proteins --cpu 16 -m diamond # specify 'braker.proteins.faa' as input (datatype=protein), set the prefix of all the output files to be 'cramerella_annot', specify the mode to use DIAMOND to perform fast sequence alignment (if you dont specify -m, its also ok as eggNOG-mapper will use DIAMOND to perform sequene alignment by default)
        """, language="bash")
        st.markdown("[Visit eggNOG-mapper GitHub Page](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.0.0-v2.0.1)")
        st.markdown("[Read eggNOG-mapper Publication](https://academic.oup.com/mbe/article/34/8/2115/3782716?login=false)")
        st.markdown("[Read latest eggNOG-mapper Publication](https://pubmed.ncbi.nlm.nih.gov/30418610/)")

# Phase 3: Molecular Docking and Dynamics Simulation
if selected == "Phase 3: Molecular Docking and Dynamics Simulation":
    with st.container():
        st.write("---")
        st.header("Phase 3: Molecular Docking and Dynamics Simulation 🖥️🧪")
        st.write("###")
        st.write("Use both AlphaFold3 & ColabFold server to model the selected 3D developmental protein structures and superimpose the two structures with UCSF ChimeraX & compare their quality metrics in the form of table. To ensure compliance with AlphaFold3 server terms and usage, only use the results of AlphaFold3 server for structural analysis and structural comparison. Use ColabFold server for docking later")
        st.write("###")
        st.write("You can use logMD to visualize the trajectory of your protein-ligand complex easily (logMD functions the same as VMD)")
        st.write("generative AI drug design method, DrugHive")
        # ----LOAD GROMACS CODE----
        # Check if the file exists before reading
        if gromacs_file.exists():
            with open(gromacs_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Gromacs Code",
                data=script_byte,
                file_name=gromacs_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{gromacs_file.name} does not exist.")

        st.write("###")

        st.write("checkout visual molecular dynamics (VMD) 2 Alpha")

        st.write("###")

        st.write("Use the pymol movie script to make a movie for your protein")
        # ----LOAD PYMOL MOVIE SCRIPT----
        # Check if the file exists before reading
        if pymol_movie_script.exists():
            with open(pymol_movie_script, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Pymol Movie Script",
                data=script_byte,
                file_name=pymol_movie_script.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{pymol_movie_script.name} does not exist.")
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
        st.code("""
        kill 103839 100000 # This is the safest kill approach to kill processes running in the background. After running this command, you have to wait for 1-2 minutes. Only then, you try to run "ps aux | grep process_name" to check whether the process has been successfully terminated
        """, language="bash")
        st.code("pkill -9 -f spades",language="bash")  # force kill all 'spades' related processes without the need to specify the ID of the process # This is very dangerous to run as it might kill dwservice if the service happens to reference to the spades process in the launch command
        st.code("kill -9 103839", language="bash")  # force kill the main spades process with the process ID, '103839'. This is the safest command to kill the running bioinformatics process

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
        st.code("""
        conda env remove --name name_of_created_environment # or conda env remove -n name_of_created_environment
        """, language="bash")

        st.write("###")

        st.write("5. Delete and clean all the files within the working directory")
        st.code("rm -rf /path/to/your/directory/*", language="bash")

        st.write("###")

        st.write("6. Remember to stop running the Docker desktop if you dont use it to save up CPU resources. Double check whether the Docker desktop already stops running in the background. After checking the status, exit by running 'q'.")
        st.code("systemctl --user stop docker-desktop", language="bash")
        st.code("systemctl --user status docker-desktop", language="bash")
        st.code("q", language="bash")

        st.write("###")

        st.write("7. If you want to show/teach others how to code within Ubuntu/Linux terminal")
        st.code("PS1='$ '", language="bash") # run this command to change how the Linux terminal looks like so that it saves & continues displaying both the previous command lines and results, allowing you to continue coding
        st.code("#", language="bash") # specify this '#" to add a comment before executing any command line/code

        st.write("###")

        st.write("8. End the Ubuntu/Linux terminal session")
        st.code("exit", language="bash")

        st.write("###")

        st.write("9. Clear the terminal")
        st.code("clear", language="bash")

        st.write("###")

        st.write("10. To check how to cite R and R packages in RStudio")
        st.code("citation()", language="r")

        st.write("###")

        st.write("11. To run code in RStudio, you can use the keyboard shortcut 'Ctrl + Enter' to run the code instead of manually clicking the 'run' button")

        st.write("###")

        st.write("12. If you want to display all of the available bioconductor packages offered in RStudio")
        st.code("BiocManager::available()", language="r")

        st.write("###")

        st.write("13. Use the str() function to examine the structure of the file that you have stored in the variable") # examine how many rows of observations, how many columns of variables, examine whether it is a dataframe, list, vector, matrix being stored in the variable, examine the data type in RStudio
        st.code("str('meta')") # assume that the name of the variable used to store the sample metadata.csv file is meta

        st.write("###")

        st.write("14. To remove all the files within a single directory (by staying within the same directory)")
        st.code("rm *", language="bash")

        st.write("###")

        st.write("15. To remove all the subdirectories and files within a single directory (by staying within the same directory")
        st.code("rm -r *", language="bash")

        st.write("###")

        st.write("16. To decompress a file (.gz)")
        st.code("""
        gunzip filename.gz # decompress a single file without keeping the original gzipped file
        gunzip -k filename.gz # decompress a single file only one at a time while keeping the original gzipped file
        gunzip -k *.gz # decompress all the gzipped files with .gz extensions within the working directory all at once while keeping the original gzipped file
        """, language="bash")

        st.write("###")

        st.write("17. To check whether a program crashes/suddenly get killed, when you suspect a hardware/memory issue, when you debug weird behaviour in WSL or Linux")
        st.code("""
        dmesg | tail -n 20 # display the last 20 lines of kernel messages to view boot process details, device connection events, memory errors, processes being killed due to OOM conditions
        """, language="bash")

        st.write("###")

        st.write("18. Print the current working directory")
        st.code("pwd", language="bash")

        st.write("###")

        st.write("19. Rename a previously created conda environment")
        st.code("""
        conda deactivate / conda activate base # make sure you are not staying within the conda environment that you wish to rename
        conda rename -n falcon pb-assembly # apply this command, "conda rename -n old_env new_env". This will change the previously created 'falcon' environment to 'pb-assembly'. Make sure that your installed conda version is 4.1.0+. You can check the installed conda version via "conda --version"
        """)

        st.write("###")

        st.write("20. Verify the installation of the bioinformatic tool that you've just installed via conda")
        st.code("""
        conda list name_of_the_installed_bioinformatic_tool # list the installed bioinformatic tool's package environment
        which name_of_the_installed_bioinformatic_tool # display the filepath of the installed bioinformatic tool
        whereis name_of_the_installed_bioinformatic_tool # it does the same thing as using the 'which' command
        conda list | grep name_of_the_bioinformatic_tool_installed # this does the same thing as the first one
        name_of_the_bioinformatic_tool_installed --version # check the version of the installed bioinformatics tool
        name_of_the_bioinformatic_tool_installed --help # display the help (user manual) of the installed bioinformatic tool
        """, language="bash")

        st.write("###")

        st.write("21. Check disk space available inside the computer")
        st.code("""
        df -h # show total disk space, used space, available space, and percentage of use for each mounted filesystem
        """, language="bash")

        st.write("###")

        st.write("22. Check memory usage available inside the computer")
        st.code("""
        free -h # display the total amount of free and used physical & swap memory in the system in a human-readable format
        """, language="bash")

        st.write("###")

        st.write("23. To check back & record the previous commands that you have run on the terminal (this is useful for your logbook record)")
        st.code("""
        history | grep busco # view the list of history of running busco commands on the terminal
        """)

        st.write("###")

        st.write("24. Display all the processes running as a cbr15 user in the background")
        st.code("ps -U cbr15 -l -H | grep ' S '", language="bash")

        st.write("###")

        st.write("25. If strict channel priority is enabled in your Conda configuration which prevents you from installing a conda package (because another channel with a higher priority doesnt contain the required conda package), then you can add '--no-channel-priority' flag while installing the required conda package")
        st.code("""
        conda install bioconda::flash --no-channel-priority # be sure to run this command within a created conda environment as a non-root user to avoid affecting the system globally
        """, language="bash")

        st.write("###")

        st.write("26. If you want to send the output into an output.log file while printing/displaying the output in the terminal at the same time, run the 'tee' flag")
        st.code("""
        command | tee output.log
        """, language="bash")

        st.write("###")

        st.write("27. Clear the terminal inside the conda environment using the below command")
        st.code("""
        which clear
        /usr/bin/clear # run this command to clear the terminal within the created conda environment
        """, language="bash")

        st.write("###")

        st.write("28. Get your user ID & group ID")
        st.code("""
        id -u # get user ID
        id -g # get group ID
        """, language="bash")

        st.write("###")

        st.write("29. List the size of each folder within the /media/raid partition")
        st.code("""
        du -h --max-depth=1 /media/raid/ | sort -rh | head -n 10 # you can safely run this command within a created isolated conda environment as a non-root user
        """, language="bash")

        st.write("###")

        st.write("###")