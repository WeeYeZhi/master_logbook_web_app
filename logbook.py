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
deseq2rmd_file = current_dir / "assets" / "Multivariate_Statistical_Analysis_of_DEG_Results_Between_CPB_Larva_and_Adult.Rmd"
gromacs_file = current_dir / "assets" / "Gromacs_codes.txt"
generate_busco_plot_python_script = current_dir / "assets" / "generate_plot.py"
generate_busco_plot_Rscript = current_dir / "assets" / "busco_figure.R"
trimmomatic_file = current_dir / "assets" / "trimmomatic.sh"
trinity_samples_sequenced_by_LKM = current_dir / "assets" / "trinity_samples_sequenced_by_LKM.txt"
trinity_samples_sequenced_by_INBIOSIS = current_dir / "assets" / "trinity_samples_sequenced_by_INBIOSIS.txt"
trinity_samples_sequenced_by_LKM_and_INBIOSIS = current_dir / "assets" / "trinity_samples_sequenced_by_LKM_and_INBIOSIS.txt"
transrate_file = current_dir / "assets" / "transrate.sh"
swissprot_fasta_file = current_dir / "assets" / "uniprot_sprot.fasta.gz"
pymol_movie_script = current_dir / "assets" / "movie01_script.pml"
CPB_pic = current_dir / "assets" / "CPB.png"
bowtie2_mapping_after_trimmomatic = current_dir / "assets" / "bowtie2_mapping_after_trimmomatic.sh"
generate_mapping_statistics_after_trimmomatic = current_dir / "assets" / "generate_mapping_statistics_after_trimmomatic.sh"
bowtie2_mapping_after_cdhitest_trinity = current_dir / "assets" / "bowtie2_mapping_after_cdhitest_trinity.sh"
generate_mapping_statistics_after_cdhitest_trinity = current_dir / "assets" / "generate_mapping_statistics_after_cdhitest_trinity.sh"
bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS = current_dir / "assets" / "bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS.sh"
generate_mapping_statistics_after_cdhitest_trinity_NCBI_FCS = current_dir / "assets" /"generate_mapping_statistics_after_cdhitest_trinity_NCBI_FCS.sh"
bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate = current_dir / "assets" / "bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.sh"
generate_mapping_statistics_after_cdhitest_trinity_NCBI_FCS_transrate = current_dir / "assets" / "generate_mapping_statistics_after_cdhitest_trinity_NCBI_FCS_transrate.sh"
extract_TM_regions = current_dir / "assets" / "TM1_to_TM7_extraction.py"
original_OctB2R_phylogenetic_tree = current_dir / "assets" / "original_octopamine_receptor_phylogenetic_tree.nwk"
rscript_for_phylogenetic_tree_reconstruction = current_dir / "assets" / "phylogenetic_tree_reconstruction_after_MEGA11.Rmd"
test_alphafold3 = current_dir / "assets" / "fold_input.json"
CPB_OctB2R = current_dir / "assets" / "CPB_OctB2R.json"
ligands = current_dir / "assets" / "ligands.smi"
ligand_batch_preparation_script = current_dir / "assets" / "batch_preparation_of_ligands.sh"
autodock4_ligand_batch_preparation_script = current_dir / "assets" / "prepare_ligand4.py"

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
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM/trinity_samples_sequenced_by_LKM.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM/trinity_assembly_done_by_LKM_output.log 2>&1 & 
        """, language="bash")
        st.write("✔️run Trinity to perform transcriptome assembly to assemble all the 6 raw transcriptomic Illumina paired end reads (sequenced by INBIOSIS) after trimming the reads via Trimmomatic")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS/trinity_samples_sequenced_by_INBIOSIS.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_INBIOSIS/trinity_assembly_done_by_INBIOSIS_output.log 2>&1 & 
        """, language="bash")
        st.write("✔️run Trinity to perform comprehensive transcriptome assembly to assemble all the 9 raw transcriptomic Illumina paired end reads (sequenced by LKM and INBIOSIS)")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_samples_sequenced_by_LKM_and_INBIOSIS.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_CPB_transcriptome_assembly > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_CPB_transcriptome_assembly/trinity_output.log 2>&1 & 
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

        st.write("Run Trinity to assemble all the Illumina paired-end reads for each sample")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266554/SRR11266554.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266554 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266554/SRR11266554_output.log 2>&1 & 
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266555/SRR11266555.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266555 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266555/SRR11266555_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266556/SRR11266556.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266556 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR11266556/SRR11266556_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690969/SRR9690969.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690969 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690969/SRR9690969_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690970/SRR9690970.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690970 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690970/SRR9690970_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690971/SRR9690971.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690971 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690971/SRR9690971_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690972/SRR9690972.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690972 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690972/SRR9690972_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690973/SRR9690973.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690973 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690973/SRR9690973_output.log 2>&1 &
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690974/SRR9690974.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690974 > /media/raid/Wee/WeeYeZhi/output/trinity_assembly_for_each_sample/trinity_SRR9690974/SRR9690974_output.log 2>&1 &
        """, language="bash")

        st.write("###")

        st.write("Build a customized diamond-indexed primates and rodent database to be used for running BLASTx for each sample of Illumina paired-end reads to check which sample has high percentage of primates and rodents")
        st.write("✔️download the protein fasta files of primate and rodents")
        st.code("""
        wget -O primates.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=(taxonomy_id:9443)"
        wget -O rodents.fasta.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=(taxonomy_id:9989)"
        """, language="bash")
        st.write("✔️unzip the protein fasta files, tag the fasta headers with PRIMATE and RODENT, and merge the two protein fasta files together")
        st.code("""
        gunzip primates.fasta.gz
        gunzip rodents.fasta.gz
        sed 's/^>/\>PRIMATE|/' primates.fasta > primates.tagged.fasta
        sed 's/^>/\>RODENT|/' rodents.fasta  > rodents.tagged.fasta
        cat primates.tagged.fasta rodents.tagged.fasta > primate_rodent.dmnd.fasta
        """, language="bash")
        st.write("✔️build the diamond database")
        st.code("""
        diamond makedb --in primate_rodent.dmnd.fasta --db primate_rodent_db --threads 16 # name the diamond database as primate_rodent_db.dmnd
        """, language="bash")
        st.write("✔️run DIAMOND BLASTx against the customized protein database (containing only primates and rodents protein sequences) for each sample")
        st.code("""
        docker run --user 1000:1000 -d --name blastx_SRR11266554_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR11266554/blastx_SRR11266554_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR11266554_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR11266555_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR11266555/blastx_SRR11266555_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR11266555_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR11266556_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR11266556/blastx_SRR11266556_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR11266556_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR9690969_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR9690969/blastx_SRR9690969_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR9690969_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR9690970_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR9690970/blastx_SRR9690970_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR9690970_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR9690971_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR9690971/blastx_SRR9690971_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR9690971_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR9690972_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR9690972/blastx_SRR9690972_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR9690972_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR9690973_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR9690973/blastx_SRR9690973_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR9690973_against_primate_rodent_results" 
        docker run --user 1000:1000 -d --name blastx_SRR9690974_tabular_results_against_primate_rodent_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/blastx_against_primate_rodent_db_results/primate_rodent_db.dmnd -q /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz -o /data/blastx_against_primate_rodent_db_results/SRR9690974/blastx_SRR9690974_against_primate_rodent.outfmt6 --outfmt 6 --fast --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_SRR9690974_against_primate_rodent_results" 
        """, language="bash")
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
        st.write("✔️remove the suffix from the assembly fasta file to keep only the transcript ID")
        st.code("""
        sed 's/ .*$//' modified_clean.fasta > latest_modified_clean.fasta
        mv latest_modified_clean.fasta modified_clean.fasta
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
        nohup busco -m transcriptome -i /media/cbr14/Two/Wee/WeeYeZhi/output/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta -c 16 -l lepidoptera_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after NCBI FCS
        nohup busco -m transcriptome -i /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -c 16 -l lepidoptera_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after removing bad contigs via TransRate (by using lepidoptera_odb12 dataset)
        nohup busco -m transcriptome -i /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -c 16 -l insecta_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after removing bad contigs via TransRate (by using insecta_odb12 dataset)
        nohup busco -m transcriptome -i /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -c 16 -l endopterygota_odb12 -o CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering > CPB_transcriptome_assembly_busco_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 & # compute the BUSCO completeness score of the CPB transcriptome assembly after removing bad contigs via TransRate (by using endopterygota_odb12 dataset)
        """, language="bash")
        st.write("✔️draw BUSCO plot for comparison using python script (built inside the busco conda package))")
        st.code("""
        which busco # to find the location of BUSCO executable 
        head -n 5 /home/cbr15/anaconda3/envs/busco/bin/generate_plot.py # display the first 5 lines to see whether there is a shebang line, "#!/usr/bin/env python3" at the first line of the file. If have, then no need to specify the flag, "python3" while running the command to draw the BUSCO plot
        nohup python3 /home/cbr14/anaconda3/envs/busco/bin/generate_plot.py -wd /media/cbr14/Two/Wee/WeeYeZhi/output/Buscoresults/BUSCO_summaries_transcriptome -rt specific --no_r > busco_plot_output.log 2>&1 & # output only the R code to be run inside the RStudio later to draw the BUSCO plot # need to specify the flag, "python3" if the file doesn't have the shebang line
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
        perl /home/cbr14/anaconda3/envs/trinity/bin/TrinityStats.pl /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta > stats_trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.log
        """, language="bash")

        st.write("###")

        st.write("**13. Perform reference-free quality assessment of de novo CPB transcriptome assembly for 4 times via TransRate (before deduplication, after deduplication, & after NCBI-FCS & after removing bad contigs via TransRate)**")
        st.write("✔️install transrate via docker container")
        st.code("""
        docker pull genevia/transrate:v1.0.3_orp # pull the transrate docker image from the DockerHub registry
        docker images # verify that the transrate docker image has been successfully pulled
        id -u && id -g # check the current user ID and group ID
        """, language="bash")
        st.write("✔️run transrate to evaluate the CPB transcriptome assembly (after NCBI-FCS) via Docker")
        st.code("""
        docker run --user 1000:1000 -d --name orp_transrate_after_NCBI_FCS -v /media/cbr14/Two/Wee/WeeYeZhi/output:/data genevia/transrate:v1.0.3_orp /bin/bash -c "transrate --assembly=/data/ncbi_fcs_results/after_tsa_submission_filtering_results/modified_clean.fasta --left=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz --right=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz --threads 16 --output=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS" 
        """, language="bash")
        st.write("✔️run transrate to evaluate the CPB transcriptome assembly (after filtering & removing bad contigs via Transrate) via Docker")
        st.code("""
        docker run --user 1000:1000 -d --name orp_transrate_after_NCBI_FCS_transrate_filtering -v /media/cbr14/Two/Wee/WeeYeZhi/output:/data genevia/transrate:v1.0.3_orp /bin/bash -c "transrate --assembly=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --left=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_1_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_1_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_1_paired.fastq.gz --right=/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266554_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266555_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR11266556_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690969_2_paired.fastq.gz, /data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690970_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690971_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690972_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690973_2_paired.fastq.gz,/data/trimmomatic_results/trimmomatic_raw_RNA_seq_from_NCBI_SRA_results/SRR9690974_2_paired.fastq.gz --threads 16 --output=/data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering" 
        """, language="bash")
        st.write("✔️show the log file content of the transrate docker that is running in the background")
        st.code("""
        docker logs orp_transrate_after_NCBI_FCS # show the log file content of the run by using the designated name of the docker container after deduplication & NCBI-FCS
        docker logs orp_transrate_after_NCBI_FCS_transrate_filtering # show the log file content of the run by using the designated name of the docker container after deduplication & NCBI-FCS & removing bad contigs via Transrate      
        docker ps -a # list all the current actively running and stopped docker containers. Those that exited with code 0 is the successful run
        """, language="bash")
        st.markdown("[Visit TransRate GitHub Page](https://github.com/blahah/transrate)")
        st.markdown("[Visit TransRate Conda Installation Page](https://anaconda.org/bioconda/transrate)")
        st.markdown("[Visit TransRate User Manual Page](https://hibberdlab.com/transrate/)")
        st.markdown("[Visit TransRate Step by Step User Manual Page](https://hibberdlab.com/transrate/getting_started.html)")

        st.write("###")

        st.write("**14. Check the total number of full length coding transcripts**")
        st.write("✔️download the UniProt SwissProt database")
        st.code("""
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        """, language="bash")
        st.write("✔️build a blastable database")
        st.code("""
        makeblastdb -in uniprot_sprot.fasta -dbtype prot
        """, language="bash")
        st.write("✔️perform the blast search, reporting only the top alignment")
        st.code("""
        blastx -help
        nohup blastx -query /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta -db /media/cbr14/Two/Wee/WeeYeZhi/output/count_full_length_transcript_results/uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 16 -max_target_seqs 1 -outfmt 6 2> error_output.log &
        """, language="bash")
        st.write("✔️examine the percentage of the target being aligned to by the best matching Trinity transcript")
        st.code("""
        nohup perl /home/cbr14/anaconda3/envs/trinity/bin/analyze_blastPlus_topHit_coverage.pl /media/cbr14/Two/Wee/WeeYeZhi/output/count_full_length_transcript_results/blastx_search_results/blastx.outfmt6 /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta /media/cbr14/Two/Wee/WeeYeZhi/output/count_full_length_transcript_results/uniprot_sprot.fasta 2> error_output.log &
        """, language="bash")
        st.write("✔️group multiple HSPs per transcript/database_match pairing to collect all the pieces of alignment between one transcript and one protein")
        st.code("""
        perl -I /home/cbr14/anaconda3/envs/trinity/bin/PerlLib /home/cbr14/anaconda3/envs/trinity/bin/blast_outfmt6_group_segments.pl /media/cbr14/Two/Wee/WeeYeZhi/output/count_full_length_transcript_results/blastx_search_results/blastx.outfmt6 /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta /media/cbr14/Two/Wee/WeeYeZhi/output/count_full_length_transcript_results/uniprot_sprot.fasta > blastx.outfmt6.grouped
        """, language="bash")
        st.write("✔️compute the coverage percentage by length to check the percentage of the protein that is covered by the transcript")
        st.code("""
        perl -I /home/cbr14/anaconda3/envs/trinity/bin/PerlLib /home/cbr14/anaconda3/envs/trinity/bin/blast_outfmt6_group_segments.tophit_coverage.pl /media/cbr14/Two/Wee/WeeYeZhi/output/count_full_length_transcript_results/blast_outfmt6_group_segments.pl_results/blastx.outfmt6.grouped > blast_outfmt6_group_segments.tophit_coverage.tsv
        """, language="bash")

        st.write("###")

        st.write("**15. Generate the gene to trans map file for the CPB transcriptome assembly**")
        st.code("""
        perl /home/cbr14/anaconda3/envs/trinity/bin/get_Trinity_gene_to_trans_map.pl /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta > /media/cbr14/Two/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/good.modified_clean.fasta.gene_trans_map
        """, language="bash")

        st.write("###")

        st.write("**16. Align the trimmed reads back to the CPB transcriptome assembly via Bowtie2 (alignment-based approach**)")
        st.write("✔️create a bowtie2 conda environment, activate the environment, and install bowtie2 via bioconda")
        st.code("""
        conda create -n bowtie2
        conda activate bowtie2
        conda install bioconda::bowtie2
        """, language="bash")
        st.write("✔️index the CPB transcriptome assembly via Bowtie2")
        st.code("""
        nohup bowtie2-build /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta /media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_build_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate > /media/cbr14/Two/Wee/WeeYeZhi/output/bowtie2_results/bowtie2_build_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate/bowtie2_build_output.log 2>&1 &
        """, language="bash")
        st.write("✔️map the trimmed reads back to the CPB transcriptome assembly via Bowtie2 by using a bash script")
        st.code("""
        dos2unix bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.sh
        chmod +x bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.sh
        nohup bash bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.sh > bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate_output.log 2>&1 &
        """, language="bash")
        # ----LOAD THE BOWTIE2 MAPPING AFTER CDHITEST TRINITY NCBI FCS TRANSRATE SCRIPT----
        # Check if the file exists before reading
        if bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.exists():
            with open(bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Bowtie2 Mapping After CDHITEST Trinity NCBI FCS TransRate Script",
                data=script_byte,
                file_name=bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{bowtie2_mapping_after_cdhitest_trinity_NCBI_FCS_transrate.name} does not exist.")
        st.markdown("[Visit Bioconda Installation of Bowtie2](https://anaconda.org/)")

        st.write("###")

        st.write("17. Run the transcript quantification step (follow the Trinity GitHub manual)")
        st.write("✔️activate the 'trinity' conda environment")
        st.code("""
        conda activate trinity
        """, language="bash")
        st.write("✔️run the align_and_estimate_abundance.pl script by running alignment-based approach via Bowtie2 and transcript quantification via RSEM")
        st.code("""
        nohup perl /home/cbr14/anaconda3/envs/trinity/bin/align_and_estimate_abundance.pl --transcripts /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --seqType fq --samples_file /media/cbr14/Two/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS/trinity_samples_sequenced_by_LKM_and_INBIOSIS.txt --est_method RSEM --aln_method bowtie2 --gene_trans_map /media/cbr14/Two/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/good.modified_clean.fasta.gene_trans_map --prep_reference --SS_lib_type RF --thread_count 16 --output_dir /media/cbr14/Two/Wee/WeeYeZhi/output/align_and_estimate_abundance.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate > trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_output.log 2>&1 &
        """, language="bash")

        st.write("###")

        st.write("**18. Compute the ExN50 statistics of the transcriptome assemblies (before deduplication, after deduplication, after NCBI FCS and after removing bad contigs via TransRate)**")
        st.write("✔️compute and obtain the expression matrix after quantifying the abundance of transcripts via RSEM (Trinity_trans.gene.TMM.EXPR.matrix will be used to run DEG analysis later & Trinity_trans.isoform.TMM.EXPR.matrix will be used to compute the EX90N50 value of the CPB transcriptome assembly later")
        st.code("""
        nohup perl /home/cbr14/anaconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map /media/cbr14/Two/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/good.modified_clean.fasta.gene_trans_map --out_prefix Trinity_trans --name_sample_by_basedir --quant_files /media/cbr14/Two/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/RSEM_transcript_quantification_files_after_cdhitest_trinity_NCBI_FCS_transrate.list > /media/cbr14/Two/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/abundance_estimates_to_matrix_after_cdhitest_trinity_NCBI_FCS_transrate_output.log 2>&1 &
        """, language="bash")
        st.write("✔️compute the ExN50 statistics of the CPB transcriptome assembly at gene level (instead of transcript level) to remove bias")
        st.code("""
        perl /home/cbr14/anaconda3/envs/trinity/bin/contig_ExN50_statistic.pl /media/cbr14/Two/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/Trinity_trans.isoform.TMM.EXPR.matrix /media/cbr14/Two/Wee/WeeYeZhi/output/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta gene | tee ExN50_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.gene.stats 
        """, language="bash")
        st.write("✔️install the ggplot2 conda package within the trinity conda environment via conda-forge")
        st.code("""
        conda install conda-forge::r-ggplot2
        """, language="bash")
        st.write("✔️download the 'plot_ExN50_statistic.Rscript' file from https://github.com/trinityrnaseq/trinityrnaseq/blob/master/util/misc/plot_ExN50_statistic.Rscript as the trinity conda package does not contain this Rscript (need to download separately)")
        st.write("✔️plot and visualize the ExN50 statistics results of the CPB transcriptome assembly")
        st.code("""
        nohup Rscript /media/cbr14/Two/Wee/WeeYeZhi/output/trinity_ExN50_gene_level_statistics_calculation_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/plot_ExN50_statistic.Rscript ExN50_after_cdhitest_trinity_NCBI_FCS_transrate_filtering.gene.stats > plot_ExN50_statistic_after_cdhitest_trinity_NCBI_FCS_transrate_filtering_output.log 2>&1 &
        """, language="bash")
        st.markdown("[Bear in mind that running abundance_estimates_to_matrix.pl script automatically generates both the transcript-level and gene-level expression matrices for you already](https://github.com/trinityrnaseq/trinityrnaseq/issues/356)")

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
        st.write("✔️download and index the HMM profile database to run hmmsearch")
        st.code("""
        wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz # download the main HMM profile database for hmmsearch from the official EBI FTP server
        gunzip Pfam-A.hmm.gz # uncompress the file
        docker run --user 1000:1000 -d --name hmmpress -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "hmmpress /data/transdecoder_results/hmmsearch_results/Pfam-A.hmm"
        """, language="bash")
        st.write("✔️extract the long open reading frames (ORFs) from the CPB transcriptome assembly")
        st.code("""
        docker run --user 1000:1000 -d --name ORF_extraction -v /media/cbr14/Two/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "TransDecoder.LongOrfs -S -t /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --output_dir /data/transdecoder_results/transdecoder_ORF_extraction_results"
        """, language="bash")
        st.write("✔️download the already preformatted BLAST NCBI NR database file (NCBI non-redundant protein database) and decompress each preformatted file by using the update_blastdb.pl script")
        st.code("""
        conda create -n blast
        conda activate blast
        conda install bioconda::blast
        which update_blastdb.pl
        update_blastdb.pl --help
        nohup perl /home/cbr14/anaconda3/envs/blast/bin/update_blastdb.pl --passive --decompress nr 2> error_output.log &
        """, language="bash")
        st.write("✔️download the NCBI NR database FASTA file from https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/, create a DIAMOND-formatted reference NR database, and run DIAMOND BLASTp")
        st.code("""
        conda create -n diamond && conda activate diamond && conda install bioconda::diamond
        nohup wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz 2> error_output.log & # download the FASTA format of the NCBI NR database file
        nohup diamond makedb --in /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_indexed_ncbi_nr_database/nr.gz --db nr --threads 16 2> /media/raid/Wee/WeeYeZhi/output/transdecoder_results/ncbi_nr_database/error_output.log &
        nohup diamond blastp --db /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_indexed_ncbi_nr_database/nr --query /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep --out /media/raid/Wee/WeeYeZhi/output/transdecoder_results/ncbi_blastp_against_nr_results/blastp.outfmt5 --outfmt 5 --ultra-sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_blastp_tmp_dir/diamond_blastp_XML 2> /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_blastp_tmp_dir/diamond_blastp_XML/error_output.log & # output BLASTp results in XML format
        nohup diamond blastp --db /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_indexed_ncbi_nr_database/nr --query /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep --out /media/raid/Wee/WeeYeZhi/output/transdecoder_results/ncbi_blastp_against_nr_results/blastp.outfmt6 --outfmt 6 --ultra-sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_blastp_tmp_dir/diamond_blastp_tabular 2> /media/raid/Wee/WeeYeZhi/output/transdecoder_results/diamond_blastp_tmp_dir/diamond_blastp_tabular/error_output.log & # output BLASTp results in tabular format
        """, language="bash")
        st.write("✔️alternatively, if you want to use a smaller protein database file like uniprot_sprot.fasta, you can download UniProt Swiss-Prot FASTA file externally from UniProt’s website (https://www.uniprot.org/downloads) under the Swiss-Prot section & create the UniProt database")
        st.code("""
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
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
        st.write("✔️run NCBI BLASTp to blast the protein sequences predicted from the CPB transcriptome assembly against the SwissProt protein database")
        st.code("""
        docker run --user 1000:1000 -d --name blastp_search_against_swissprot_output_tabular -v /media/cbr14/Two/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "blastp -query /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep -db /data/transdecoder_results/make_swissprot_database/uniprot_sprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > /data/transdecoder_results/ncbi_blastp_against_swissprot_results/blastp.outfmt6"
        """, language="bash")
        st.write("✔️run NCBI BLASTp to blast the protein sequences predicted from the CPB transcriptome assembly against the NCBI NR database")
        st.code("""
        docker run --user 1000:1000 -d --name blastp_search_against_nr_output_tabular -v /media/cbr14/Two/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "blastp -query /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep -db /data/transdecoder_results/ncbi_nr_database/nr -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > /data/transdecoder_results/ncbi_blastp_against_nr_results/blastp.outfmt6" # output the tabular format of BLASTp results
        docker run --user 1000:1000 -d --name blastp_search_against_nr_output_XML -v /media/cbr14/Two/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "blastp -query /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep -db /data/transdecoder_results/ncbi_nr_database/nr -max_target_seqs 1 -outfmt 5 -evalue 1e-5 -num_threads 16 > /data/transdecoder_results/ncbi_blastp_against_nr_results/blastp.outfmt5" # output the XML format of the BLASTp results
        """, language="bash")
        st.write("✔️run HMMSEARCH to identify the protein domains in the predicted ORFs using the Pfam database")
        st.code("""
        docker run --user 1000:1000 -d --name hmmsearch -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "hmmsearch --cpu 16 -E 1e-5 --domtblout /data/transdecoder_results/hmmsearch_results/pfam.domtblout /data/transdecoder_results/hmmsearch_results/Pfam-A.hmm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder_dir/longest_orfs.pep"
        """, language="bash")
        st.write("✔️filter and retain predicted ORFs that have supported BLASTp hits and Pfam domains")
        st.code("""
        docker run --user 1000:1000 -d --name ORF_prediction -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/transdecoder /bin/bash -c "TransDecoder.Predict -t /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --retain_pfam_hits /data/transdecoder_results/hmmsearch_results/pfam.domtblout --retain_blastp_hits /data/transdecoder_results/ncbi_blastp_against_nr_results/blastp.outfmt6 -O /data/transdecoder_results/transdecoder_ORF_extraction_results"
        """, language="bash")
        st.markdown("[Read how to download NCBI NR, swissprot, and nt database files before running BLAST](https://www.ncbi.nlm.nih.gov/books/NBK62345/)")
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
        docker run --rm trinityrnaseq/trinotate find / -type f -name Trinotate 2>/dev/null # search for the full path to the Trinotate executable file
        docker run --rm trinityrnaseq/trinotate /usr/local/src/Trinotate/Trinotate -h # print the command-line options available to run Trinotate correctly
        """, language="bash")
        st.write("✔️start running the Trinotate docker container (with the container being named as 'trinotate') and keep it alive in the background, and set the environment variables")
        st.code("""
        docker run -it --user 0:0 --name trinotate -v /media/raid/Wee/WeeYeZhi/output:/data -e TRINOTATE_HOME=/usr/local/src/Trinotate -e TRINOTATE_DATA_DIR=/data/trinotate_results/trinotate_data_dir trinityrnaseq/trinotate /bin/bash # remember to run docker stop and docker rm after running trinotate, otherwise it is going to consume resources
        """, language="bash")
        st.write("✔️connect back to the running Trinotate docker container (run this command everytime after you exit the terminal to connect back to the running Trinotate docker container)")
        st.code("""
        docker exec -it trinotate /bin/bash
        """, language="bash")
        st.write("✔️build the trinotate sqlite database")
        st.code("""
        $TRINOTATE_HOME/Trinotate --create --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --trinotate_data_dir /data/trinotate_results/trinotate_data_dir --use_diamond # bear in mind that the build process will help you to automatically download the resources, namely 'uniprot_sprot.pep', 'Pfam-A.hmm.gz', 'Rfam.cm', 'Rfam.clanin',and 'eggnog.db.gz' 
        """, language="bash")
        st.write("✔️prepare both the BLASTp and Pfam databases that have been previously created by the built process")
        st.code("""
        makeblastdb -in \$TRINOTATE_DATA_DIR/uniprot_sprot.pep -dbtype prot # no need to run this command already as the previous Trinotate sqlite database creation step already ran this command
        gunzip \$TRINOTATE_DATA_DIR/Pfam-A.hmm.gz  # no need to run this command already as the previous Trinotate sqlite database creation step already ran this command
        hmmpress \$TRINOTATE_DATA_DIR/Pfam-A.hmm  # no need to run this command already as the previous Trinotate sqlite database creation step already ran this command
        """, language="bash")
        st.write("✔️import CPB transcriptome and protein data into the Trinotate sqlite database to initialize the trinotate sqlite database")
        st.code("""
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --init --gene_trans_map /data/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/good.modified_clean.fasta.gene_trans_map --transcript_fasta /data/transrate_results/transrate_CPB_transcriptome_assembly_after_cdhitest_trinity_NCBI_FCS/modified_clean/good.modified_clean.fasta --transdecoder_pep /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep
        """, language="bash")
        st.write("✔️run diamond in BLASTx mode to blast translated nucleotide against protein (perform local blastx)")
        st.code("""
        sed -E 's/>(TRINITY[^ ]+_i[0-9]+)\.p[0-9]+.*/>\1/' /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.cds > /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.cds # trim off the long coding FASTA sequence headers to retain only the transcript IDs
        diamond blastx -d $TRINOTATE_DATA_DIR/uniprot_sprot.dmnd -q /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.cds -o /data/trinotate_results/diamond_blastx_results/blastx.outfmt6 -f 6 -p 16 -e 1e-5 # run BLASTx against the SwissProt database
        docker run --user 1000:1000 -d --name blastx_tabular_results_against_nr_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/transdecoder_results/diamond_indexed_ncbi_nr_database/nr.dmnd -q /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.cds -o /data/trinotate_results/diamond_blastx_results/blastx.outfmt6 --outfmt 6 --ultra-sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_tabular_results_of_trinotate" # run BLASTX against the NR database to output the tabular format of the BLASTX results
        docker run --user 1000:1000 -d --name blastx_XML_results_against_nr_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastx -d /data/transdecoder_results/diamond_indexed_ncbi_nr_database/nr.dmnd -q /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.cds -o /data/trinotate_results/diamond_blastx_results/blastx.outfmt5 --outfmt 5 --ultra-sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastx_XML_results_of_trinotate" # run BLASTX against the NR database to output the XML format of the BLASTX results
        """, language="bash")
        st.write("✔️run BLASTp & hmmsearch again using the files generated by trinotate sqlite database construction step (perform local blastp)")
        st.code("""
        sed -E 's/^(>[^ ]*).*/\1/' /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep # trim off the long protein FASTA sequence headers to retain only protein IDs in order to avoid encountering "OSError: [Errno 36] File name too long"
        docker run --user 1000:1000 -d --name blastp_tabular_results_against_nr_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastp -d /data/transdecoder_results/diamond_indexed_ncbi_nr_database/nr.dmnd -q /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep -o /data/trinotate_results/diamond_blastp_results/blastp.outfmt6 --outfmt 6 --ultra-sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastp_tabular_results_of_trinotate" # output the tabular format of the BLASTp results
        docker run --user 1000:1000 -d --name blastp_XML_results_against_nr_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastp -d /data/transdecoder_results/diamond_indexed_ncbi_nr_database/nr.dmnd -q /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep -o /data/trinotate_results/diamond_blastp_results/blastp.outfmt5 --outfmt 5 --ultra-sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastp_XML_results_of_trinotate" # output the XML format of the BLASTp results
        docker run --user 1000:1000 -d --name hmmsearch -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "hmmsearch --cpu 16 -E 1e-5 --domtblout /data/trinotate_results/hmmsearch_results/pfam.domtblout /data/trinotate_results/trinotate_data_dir/Pfam-A.hmm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep"
        """, language="bash")
        st.write("✔️run interproscan to scan the predicted CPB protein sequences against the protein family, domain, function, and binding sites classification database to try annotate the predicted CPB protein sequences")
        st.code("""
        docker pull interpro/interproscan:5.76-107.0 # pull the interproscan docker image with a specific tag of 5.76-107.0 from DockerHub 
        docker images # verify whether the interproscan docker image has been successfully pulled
        nohup curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/alt/interproscan-data-5.76-107.0.tar.gz 2> download_interproscan_database_error_output.log & # download the gzipped file of the interproscan database
        curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/alt/interproscan-data-5.76-107.0.tar.gz.md5 # download the .md5 extension of the interproscan database 
        md5sum -c interproscan-data-5.76-107.0.tar.gz.md5 # verify whether the download is complete to make sure the download inteproscan database file is not corrupted
        tar -pxzf interproscan-data-5.76-107.0.tar.gz # extract the downloaded interproscan database to interproscan-5.76-107.0/data
        sed '/^>/! s/*//g' /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed_for_interproscan.fasta.transdecoder.pep # remove all the asterisk signs (which represent stop codons inside the fasta sequence lines) before running interproscan since interproscan expects the input protein file to only include known amino acid letters
        docker run --user 1000:1000 -d --name interproscan -v /media/raid/Wee/WeeYeZhi/output:/data -v /media/raid/Wee/WeeYeZhi/output/trinotate_results/interproscan_results/interproscan_database/interproscan-5.76-107.0/data:/opt/interproscan/data interpro/interproscan:5.76-107.0 --input /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed_for_interproscan.fasta.transdecoder.pep --seqtype p --formats TSV,XML --cpu 16 --goterms --iprlookup --pathways --tempdir /data/trinotate_results/interproscan_results/interproscan_temp_dir --output-dir /data/trinotate_results/interproscan_results
        """, language="bash")
        st.write("✔️create an isolated conda environment called signalp6, activate the conda environment, install the gzipped file of SignalP6 at https://services.healthtech.dtu.dk/services/SignalP-6.0/, verify the installation and run SignalP6 to predict the signal peptide molecules of the predicted CPB protein sequences")
        st.code("""
        conda create -n signalp6 python=3.8 -y
        conda activate signalp6
        sed -E 's/^(>[^ ]*).*/\1/' /media/cbr14/Two/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.fasta.transdecoder.pep > /media/cbr14/Two/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep # trim off the long protein FASTA sequence headers to retain only protein IDs in order to avoid encountering "OSError: [Errno 36] File name too long"
        tar -xzvf signalp-6.0h.fast.tar.gz # extract the signalp6 file
        cd /media/cbr14/Two/Wee/WeeYeZhi/output/trinotate_results/signalp6_results/signalp6_installation/signalp6_fast
        pip install torch==1.8.1+cpu torchvision==0.9.1+cpu torchaudio==0.8.1 -f https://download.pytorch.org/whl/torch_stable.html
        pip install matplotlib tqdm
        pip install signalp-6-package/
        SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))") && mkdir -p "$SIGNALP_DIR/model_weights" && cp -r /media/cbr14/Two/Wee/WeeYeZhi/output/trinotate_results/signalp6_results/signalp6_installation/signalp6_fast/signalp-6-package/models/* "$SIGNALP_DIR/model_weights/"
        which signalp6
        nohup signalp6 --fastafile /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep --organism eukarya --mode fast --format txt --output_dir /media/raid/Wee/WeeYeZhi/output/trinotate_results/signalp6_results > /media/raid/Wee/WeeYeZhi/output/trinotate_results/signalp6_results/signalp6_output.log 2>&1 & # After running this command, you will automatically obtain a final output file, known as prediction_results.txt file
        """, language="bash")
        st.write("✔️install the gzipped file of TMHMM at https://services.healthtech.dtu.dk/services/TMHMM-2.0/, verify the installation and run TMHMM to predict the transmembrane domains of the predicted CPB protein sequences")
        st.code("""
        tar -xzvf tmhmm-2.0c.Linux.tar.gz # extract the tmhmm file
        cd /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin # navigate into this working directory before you run tmhmm
        which perl # check the working directory to perl for your system
        nano tmhmm # change the shebang line of the tmhmm program perl script to /home/cbr14/opt/Colabfold/localcolabfold/colabfold-conda/bin/perl
        nano tmhmmformat.pl # change the shebang line of the tmhmmformat program perl script to /home/cbr14/opt/Colabfold/localcolabfold/colabfold-conda/bin/perl
        export PATH="/media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin:$PATH" # note that the export of the PATH environment variable to the tmhmm executable file is just temporary and will only take effect for the current terminal session
        which tmhmm # double check the system can access the tmhmm executable file correctly after exporting the PATH variable 
        tmhmm --short /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_results/tmhmm.out
        """, language="bash")
        st.write("✔️use the eggNOG-mapper program that comes together with the Trinotate Docker image and run it to assign orthologs/evolutionary functional relationships to the predicted protein sequences of CPB")
        st.code("""
        docker run -it --user 0:0 --name eggnogmapper -v /media/raid/Wee/WeeYeZhi/output:/data -e EGGNOGMAPPER=/usr/local/src/eggnog-mapper-2.1.9 trinityrnaseq/trinotate /bin/bash # Note that the previous Trinotate docker image pull process already came with the pre-built eggnogmapper package. You no longer need to install eggnogmapper externally. Remember to run docker stop and docker rm after running eggnogmapper, otherwise it is going to consume resources
        $EGGNOGMAPPER/emapper.py -i /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep -m diamond --ultra-sensitive --data_dir /data/trinotate_results/trinotate_data_dir/EGGNOG_DATA_DIR --cpu 16 --itype proteins -o CPB --output_dir /data/trinotate_results/eggnogmapper_results # run eggNOGmapper to assign orthologs to the predicted protein sequences of CPB. Note that the eggnogmapper database was already previously downloaded and saved in the EGGNOG_DATA_DIR during the Trinotate sqlite database creation step
        docker run --user 0:0 -d --name eggnogmapper -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "/usr/local/src/eggnog-mapper-2.1.9/emapper.py -i /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep -m diamond --sensmode ultra-sensitive --data_dir /data/trinotate_results/trinotate_data_dir/EGGNOG_DATA_DIR --cpu 16 --itype proteins -o CPB --output_dir /data/trinotate_results/eggnogmapper_results" # Alternatively, you can run this command at one go to run the docker container in the background with the most sensitive parameter of DIAMOND search
        """, language="bash")
        st.write("✔️use the infernal program that comes together with the Trinotate Docker image and run it")
        st.code("""
        wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz # download the Rfam covariance model database via ftp (just in case if the docker image does not contain it). Note that this command was already run by the previous trinotate sqlite database creation step
        wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin # download the Rfam clan database via ftp (just in case if the docker image does not contain it).  Note that this command was already run by the previous trinotate sqlite database creation step
        gunzip Rfam.cm.gz # unzip the downloaded Rfam CM database. Note that this command was already run by the previous trinotate sqlite database creation step
        docker run -it --user 0:0 --name infernal -v /media/raid/Wee/WeeYeZhi/output:/data -e INFERNAL_CMSCAN=/usr/local/src/infernal-1.1.2/src trinityrnaseq/trinotate /bin/bash # Note that the previous Trinotate docker image pull process already came with the pre-built infernal package, hence you can just use its cmscan programme. You no longer need to install infernal externally. Remember to run docker stop and docker rm after running infernal, otherwise it is going to consume resources
        $INFERNAL_CMSCAN/cmscan --cut_ga --rfam --nohmmonly --clanin /data/trinotate_results/trinotate_data_dir/Rfam.clanin --oskip --fmt 2 -o /data/trinotate_results/infernal_results/infernal_output.txt --tblout /data/trinotate_results/infernal_results/infernal.out /data/trinotate_results/trinotate_data_dir/Rfam.cm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.cds # scan the predicted CPB coding regions of the transcript sequences against the CM models to check whether any of your predicted coding sequences match known RNA families from Rfam. Later, you need to filter and remove all the hits (detected non-coding RNAs) if there are any since you are only interested in the mRNA coding regions. Run infernal as a QC check to filter non-coding RNAs from the coding regions of the transcriptome
        docker run --user 0:0 -d --name infernal -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "/usr/local/src/infernal-1.1.2/src/cmscan --cut_ga --rfam --nohmmonly --clanin /data/trinotate_results/trinotate_data_dir/Rfam.clanin --oskip --fmt 2 -o /data/trinotate_results/infernal_results/infernal_output.txt --tblout /data/trinotate_results/infernal_results/infernal.out /data/trinotate_results/trinotate_data_dir/Rfam.cm /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.cds" # Alternatively, you can run this infernal command at one go in the background
        """, language="bash")
        st.write("✔️Load the BLASTp, BLASTx, Pfam, SignalP6, TMHMM, Infernal & eggNOGmapper results into the created Trinotate SQLite Database")
        st.code("""
        $TRINOTATE_HOME/Trinotate -h # display all the command-line options available to run Trinotate correctly
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_swissprot_blastp /data/trinotate_results/diamond_blastp_results/blastp.outfmt6
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_swissprot_blastx /data/trinotate_results/diamond_blastx_results/blastx.outfmt6
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_pfam /data/trinotate_results/hmmsearch_results/pfam.domtblout
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_signalp /data/trinotate_results/signalp6_results/prediction_results.txt
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_tmhmmv2 /data/trinotate_results/tmhmm_results/tmhmm.out
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_infernal /data/trinotate_results/infernal_results/infernal.out
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --LOAD_EggnogMapper /data/trinotate_results/eggnogmapper_results/CPB.emapper.annotations
        """, language="bash")
        st.write("✔️Generate the CPB transcriptome annotation report")
        st.code("""
        $TRINOTATE_HOME/Trinotate --db /data/trinotate_results/trinotate_sqlite_database/myTrinotate.sqlite --report -E 1e-5 > /data/trinotate_results/trinotate_annotation_report/trinotate_annotation_report.xlsx
        """, language="bash")
        st.write("✔️Extract GO assignments for the predicted CPB protein sequences")
        st.code("""
        docker run -it --user 0:0 --name GO_extraction -v /media/raid/Wee/WeeYeZhi/output:/data -e GO_EXTRACTION_SCRIPT=/usr/local/src/Trinotate/util trinityrnaseq/trinotate /bin/bash # remember to run docker stop and docker rm after running the GO extraction process, otherwise it is going to consume resources
        $GO_EXTRACTION_SCRIPT/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls /data/trinotate_results/trinotate_annotation_report/trinotate_annotation_report.xls -G --include_ancestral_terms > /data/trinotate_results/GO_annotation_results/GO_annotation.txt
        """, language="bash")
        st.markdown("[Read how to run DIAMOND BLASTx and BLASTp](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options)")
        st.markdown("[Visit Trinotate GitHub Wikipedia](https://github.com/Trinotate/Trinotate/wiki)")
        st.markdown("[Read how to load the results into the Trinotate SQLite Database](https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinotate-Functional-Annotation)")
        st.markdown("[Read how to cite Trinotate and all of the other tools](https://github.com/Trinotate/Trinotate/wiki/Lit-References)")
        st.markdown("[Visit Trinotate GitHub User Manual Page](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)")
        st.markdown("[Read how to run DeepTHMMM locally and via Docker](https://dtu.biolib.com/DeepTMHMM)")
        st.markdown("[Visit TMHMM GitHub User Manual Page](https://github.com/dansondergaard/tmhmm.py)")
        st.markdown("[Visit SignalP6 GitHub User Manual Page](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md)")
        st.markdown("[Visit TMHMM User Manual & Download Page](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)")
        st.markdown("[Visit SignalP6 User Manual & Download Page](https://services.healthtech.dtu.dk/services/SignalP-6.0/)")
        st.markdown("[Learn how to run SignalP6 and TMHMM](https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Secreted_Protein_Prediction_with_SignalP_and_TMHMM.html#gsc.tab=0)")
        st.markdown("[Read Infernal Documentation](http://eddylab.org/infernal/)")
        st.markdown("[Download the latest covariance model files (.cm) & clan files (.clanin) manually from the Rfam database](https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/)")
        st.markdown("[How to solve the issue of getting an empty output results file/columns](https://github.com/Trinotate/Trinotate.github.io/issues/70)")
        st.markdown("[Learn how to run Trinotate to annotate transcripts](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/annotating-transcripts.html#gsc.tab=0)")
        st.markdown("[Read how to run TMHMM using the TMHMM Server](https://medium.com/@snippetsbio/tmhmm-server-usage-and-result-analysis-in-3-simple-steps-4a37edfead39)")
        st.markdown("[Read how to run the whole DESeq2 pipeline](https://docs.hpc.ufs.ac.za/training/transcriptomics_tutorial/transcriptomics_tutorial_1/)")
        st.markdown("[Read how to run interproscan to annotate protein sequences](https://interproscan-docs.readthedocs.io/en/v5/HowToRun.html)")
        st.markdown("[Read how to pull interproscan docker image from DockerHub](https://hub.docker.com/r/interpro/interproscan)")

        st.write("###")

        st.write("**3. Perform differential expression gene (DEG) analysis by running the run_DE_analysis.pl (trinity utility perl script) (considering all the 4 larva replicates and 4 adult replicates)**")
        st.write("✔️convert abundance estimates into matrices by using the abundance_estimates_to_matrix.pl script")
        st.code("""
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map /media/raid/Wee/WeeYeZhi/output/get_Trinity_gene_to_trans_map_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/good.modified_clean.fasta.gene_trans_map --out_prefix Trinity --name_sample_by_basedir --quant_files /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/RSEM_transcript_quantification_files_after_cdhitest_trinity_NCBI_FCS_transrate.list > /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/abundance_estimates_to_matrix_after_cdhitest_trinity_NCBI_FCS_transrate_output.log 2>&1 &
        """, language="bash")
        st.write("✔️run DEG by using the trinity per script, run_DE_analysis.pl")
        st.code("""
        R --version # ensure R is installed within the Linux environment first before you run the script
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/run_DE_analysis.pl --matrix /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/Trinity.gene.counts.matrix --method DESeq2 --samples_file /media/raid/Wee/WeeYeZhi/output/run_DE_analysis.pl_results/samples_for_DEG_analysis.txt --reference_sample larva --output /media/raid/Wee/WeeYeZhi/output/run_DE_analysis.pl_results 2> error_output_run_DE_analysis.log &
        """, language="bash")
        st.write("✔️extract and cluster differentially expressed transcripts")
        st.code("""
        nohup perl /home/cbr15/anaconda3/envs/trinity/bin/analyze_diff_expr.pl --matrix /media/raid/Wee/WeeYeZhi/output/trinity_abundance_estimates_to_matrix.pl_results/trinity_assembly_after_cdhitest_trinity_NCBI_FCS_transrate_filtering/Trinity.gene.TMM.EXPR.matrix -P 1e-3 -C 2 --output DEG --samples /media/raid/Wee/WeeYeZhi/output/run_DE_analysis.pl_results/samples_for_DEG_analysis.txt 2> error_output_analyze_diff_expr.log & # make sure to run this command within the same directory containing the results of run_DE_analysis.pl script
        """, language="bash")
        st.write("✔️create a R markdown script in RStudio to draw all the 4 plots, namely, volcano plot, hierarchical clustering heatmap, PCA plot and MA plot automatically")
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
        st.write("✔️Alternatively, you can also create heatmap plot and perform hierarchical clustering analysis using Morpheus (if you dont want to code)")
        st.markdown("[Visit Morpheus Broad Webpage](https://software.broadinstitute.org/morpheus/)")
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

# Phase 3: Molecular Docking and Dynamics Simulation
if selected == "Phase 3: Molecular Docking and Dynamics Simulation":
    with st.container():
        st.write("---")
        st.header("Phase 3: Molecular Docking and Dynamics Simulation 🖥️🧪")
        st.write("###")

        st.write("**1. Build a customized diamond-indexed Plutella xylostella octopamine receptor database to be used for running BLASTp to blastp the CPB transcriptome assembly against the customized protein database to detect & identify the presence of octopamine receptor amino acid sequences harbouring within the CPB transcriptome**")
        st.write("✔️download the protein fasta files of OAMB, OctB1, OctB2, OctB3, and TAR-OctR of Plutella xylostella from NCBI protein database")
        st.code("""
        wget -O Px_octopamine_receptors.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=XP_048486473.1,XP_011557485.1,XP_011568733.1,XP_011556023.1,XP_048485348.1&rettype=fasta&retmode=text"
        """, language="bash")
        st.write("✔️check the downloaded protein fasta file content")
        st.code("""
        grep ">" Px_octopamine_receptors.fasta 
        head -n 100 Px_octopamine_receptors.fasta
        """, language="bash")
        st.write("✔️build the diamond database")
        st.code("""
        diamond makedb --in Px_octopamine_receptors.fasta --db Px_octopamine_receptor_db --threads 16 
        """, language="bash")
        st.write("✔️check the version of diamond and blastp used in the trinityrnaseq/trinotate docker container")
        st.code("""
        docker run --rm trinityrnaseq/trinotate diamond --version
        docker run --rm trinityrnaseq/trinotate blastp --version
        """, language="bash")
        st.write("✔️run DIAMOND BLASTp of the CPB transcriptome assembly against the customized protein database containing only Plutella xylostella octopamine receptor sequences")
        st.code("""
        docker run --user 1000:1000 -d --name blastp_CPB_transcriptome_tabular_results_against_PxOAR_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastp -d /data/blastp_against_PxOAR_db_results/Px_octopamine_receptor_db.dmnd -q /data/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep -o /data/blastp_against_PxOAR_db_results/blastp_CPB_transcriptome_tabular_results_against_PxOAR_database.outfmt6 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --threads 16 --max-target-seqs 1 --evalue 1e-5 --tmpdir /data/tmp_dir/blastp_CPB_transcriptome_tabular_results_against_PxOAR_database" 
       """, language="bash")


        st.write("###")

        st.write("**2. Extract the filtered CPB octopamine receptor sequences from the predicted CPB peptide sequences via SeqKit**")
        st.write("✔️create an ID list of the filtered CPB octopamine receptor sequences")
        st.code("""
        echo -e "TRINITY_DN108405_c0_g1_i3.p1\nTRINITY_DN114314_c1_g1_i1.p1\nTRINITY_DN219751_c0_g1_i1.p1\nTRINITY_DN22425_c0_g1_i3.p1" > CPB_octopamine_receptor_id.txt # navigate into the working directory, /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results, and run the command to save the file inside this working directory
        """, language="bash")
        st.write("✔️extract the filtered CPB octopamine receptor sequences")
        st.code("""
        seqkit grep -f /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results/CPB_octopamine_receptor_id.txt /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results/filtered_CPB_octopamine_receptors.fasta
        """, language="bash")

        st.write("###")

        st.write("**3. Run TMHMM to predict the transmembrane domains of the filtered CPB octopamine receptor sequences to check whether the CPB octopamine receptor sequences have 7 transmembrane domains of GPCR. Retain only those with 7 transmembrane helices (TMHs) for downstream analysis**")
        st.code("""
        cd /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin # navigate into this working directory before you run tmhmm
        export PATH="/media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin:$PATH" # note that the export of the PATH environment variable to the tmhmm executable file is just temporary and will only take effect for the current terminal session
        which tmhmm # double check the system can access the tmhmm executable file correctly after exporting the PATH variable 
        tmhmm /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results/filtered_CPB_octopamine_receptors.fasta > /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results/tmhmm_results/tmhmm_CPB_octopamine_receptors.out
        """, language="bash")

        st.write("###")

        st.write("**4. Detect the presence of DRY and NPxxY motifs present within the filtered CPB octopamine receptor**")
        st.code("""
        seqkit locate -p "DRY" /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results/filtered_CPB_octopamine_receptors.fasta > DRY_motif_list_CPB_octopamine_receptors.txt # search for the DRY (Asp–Arg–Tyr) amino acid portion that is often conserved among Class A (Rhodopsin-like) GPCRs & check whether the DRY motif is located at the end of the TM3 region
        seqkit locate -r -p "NP..Y" /media/raid/Wee/WeeYeZhi/output/blastp_against_PxOAR_db_results/filtered_CPB_octopamine_receptors.fasta > NPxxY_motif_list_CPB_octopamine_receptors.txt # search for the NPxxY (Asn-Pro-x-x-Tyr) amino acid portion that is often conserved among Class A (Rhodopsin-like) GPCRs & check whether the NPxxy motif is located within the TM7 region of the receptor. x within the "NPxxY" can be any amino acids
        """, language="bash")

        st.write("###")

        st.write("**5. Build a customized diamond-indexed octopamine beta 2 receptor database to be used for running BLASTp to blastp the CPB octopamine beta 2 receptor against the customized protein database to compute the similarity percentage**")
        st.write("✔️download the protein fasta files of PxOctβ2R, VdOctβ2R, DmOctβ2R, AmOctβ2R, TcOctβ2R, PrOctβ2R, CsOctβ2R, BmOctβ2R, and NlOctβ2R from the NCBI protein database")
        st.code("""
        wget -O octopamine_beta2_receptors.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=XP_011568733.1,XP_022664697.1,NP_001034049.1,XP_396348.4,NP_001280501.1,XP_022116353.1,AEO89318.1,NP_001171666.1,ASA47149.1&rettype=fasta&retmode=text"
        """, language="bash")
        st.write("✔️check the downloaded protein fasta file content")
        st.code("""
        grep ">" octopamine_beta2_receptors.fasta 
        head -n 100 octopamine_beta2_receptors.fasta
        """, language="bash")
        st.write("✔️extract the CPB octopamine beta 2 receptor")
        st.code("""
        seqkit grep -p "TRINITY_DN22425_c0_g1_i3.p1" /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/CPB_OctB2R_TRINITY_DN22425_c0_g1_i3_p1.fasta
        """, language="bash")
        st.write("✔️build the diamond database")
        st.code("""
        diamond makedb --in octopamine_beta2_receptors.fasta --db octopamine_beta2_receptor_db --threads 16 
        """, language="bash")
        st.write("✔️run DIAMOND BLASTp of the CPB octopamine beta 2 receptor against the customized protein database containing PxOctβ2R, VdOctβ2R, DmOctβ2R, AmOctβ2R, TcOctβ2R, PrOctβ2R, CsOctβ2R, BmOctβ2R, and NlOctβ2R")
        st.code("""
        docker run --user 1000:1000 -d --name blastp_CPB_octopamine_beta2_receptor_tabular_results_against_beta2_receptor_database -v /media/raid/Wee/WeeYeZhi/output:/data trinityrnaseq/trinotate /bin/bash -c "diamond blastp -d /data/blastp_CPB_beta2_receptor_against_beta2_receptor_database_results/octopamine_beta2_receptor_db.dmnd -q /data/MSA_octopamine_beta2_receptors/CPB_OctB2R_TRINITY_DN22425_c0_g1_i3_p1.fasta -o /data/blastp_CPB_beta2_receptor_against_beta2_receptor_database_results/blastp_CPB_beta2_receptor_against_beta2_receptor_database.outfmt6 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --threads 16 --evalue 1e-5 --tmpdir /data/tmp_dir/blastp_CPB_beta2_receptor_against_beta2_receptor_database" 
        """, language="bash")

        st.write("###")

        st.write("**6. Perform multiple sequence alignment (MSA) to align the CPB octopamine beta 2 receptor sequence with the octopamine beta 2 receptor sequences of Apis mellifera, Drosophila melanogaster, Tribolium castaneum, Varroa destructor, Pieris rapae, Chilo suppressalis, Bombyx mori, and Nilaparvata lugens via MUSCLE**")
        st.write("✔️create a conda environment named 'muscle', activate the conda environment, and install muscle via Bioconda")
        st.code("""
        conda create -n muscle
        conda activate muscle
        conda install bioconda::muscle
        muscle --version
        """, language="bash")
        st.write("✔️download the protein fasta files of PxOctβ2R, VdOctβ2R, DmOctβ2R, AmOctβ2R, TcOctβ2R, PrOctβ2R, CsOctβ2R, BmOctβ2R, and NlOctβ2R from the NCBI protein database")
        st.code("""
        wget -O octopamine_beta2_receptors.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=XP_011568733.1,XP_022664697.1,NP_001034049.1,XP_396348.4,NP_001280501.1,XP_022116353.1,AEO89318.1,NP_001171666.1,ASA47149.1&rettype=fasta&retmode=text"
        """, language="bash")
        st.write("✔️extract the CPB octopamine beta 2 receptor sequence from the CPB peptide sequences (previously predicted by TransDecoder) via SeqKit")
        st.code("""
        seqkit grep -p "TRINITY_DN22425_c0_g1_i3.p1" /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/CPB_OctB2R_TRINITY_DN22425_c0_g1_i3_p1.fasta
        """, language="bash")
        st.write("✔️concatenate the CPB octopamine beta 2 receptor sequence with the other insect OctB2R sequences")
        st.code("""
        cat /media/raid/Wee/WeeYeZhi/output/blastp_CPB_beta2_receptor_against_beta2_receptor_database_results/octopamine_beta2_receptors.fasta /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/CPB_OctB2R_TRINITY_DN22425_c0_g1_i3_p1.fasta > /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/combined_octopamine_beta2_receptors.fasta
        """, language="bash")
        st.write("✔️rename the FASTA sequence headers")
        st.code("""
        sed 's/^>XP_396348.4.*/>AmOct\u03B22R/; s/^>NP_001034049.1.*/>DmOct\u03B22R/; s/^>NP_001280501.1.*/>TcOct\u03B22R/; s/^>XP_011568733.1.*/>PxOct\u03B22R/; s/^>XP_022664697.1.*/>VdOct\u03B22R/; s/^>TRINITY_DN22425_c0_g1_i3.p1.*/>CcOct\u03B22R/; s/^>XP_022116353.1.*/>PrOct\u03B22R/; s/^>AEO89318.1.*/>CsOct\u03B22R/; s/^>NP_001171666.1.*/>BmOct\u03B22R/; s/^>ASA47149.1.*/>NlOct\u03B22R/' /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/combined_octopamine_beta2_receptors.fasta > /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/renamed_combined_octopamine_beta2_receptors.fasta
        """, language="bash")
        st.write("✔️run MSA via MUSCLE")
        st.code("""
        nohup muscle -align /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/renamed_combined_octopamine_beta2_receptors.fasta -output /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/MSA_between_CPB_and_other_insects_OctB2R_results.aln -threads 16 -consiters 2 -refineiters 100 > /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/MSA_between_CPB_and_other_insects_OctB2R_output.log 2>&1 &
        """, language="bash")
        st.write("✔️run TMHMM to predict the transmembrane domains of the all the octopamine beta 2 receptor sequences to check whether they have 7 transmembrane domains of GPCR. Retain only those with 7 transmembrane helices (TMHs) for downstream analysis**")
        st.code("""
        cd /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin # navigate into this working directory before you run tmhmm
        export PATH="/media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin:$PATH" # note that the export of the PATH environment variable to the tmhmm executable file is just temporary and will only take effect for the current terminal session
        which tmhmm # double check the system can access the tmhmm executable file correctly after exporting the PATH variable 
        tmhmm /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/renamed_combined_octopamine_beta2_receptors.fasta > /media/raid/Wee/WeeYeZhi/output/MSA_octopamine_beta2_receptors/tmhmm_octopamine_beta2_receptors_results/tmhmm_octopamine_beta2_receptors_results.out
        """, language="bash")
        st.markdown("[Visit MUSCLE GitHub Page](https://github.com/rcedgar/muscle?tab=readme-ov-file)")
        st.markdown("[Visit MUSCLE GitHub Manual Page](https://drive5.com/muscle5/manual/commands.html)")

        st.write("###")

        st.write("**7. Construct a phylogenetic tree to cluster and verify the identity of the detected CPB octopamine receptor subtypes**")
        st.write("✔️download all the 41 representative members of the different octopamine receptors (OARs) subtypes from Lepidoptera, Diptera, Coleoptera, Hymenoptera, Hemiptera, Primates, and Mesostigmata orders")
        st.code("""
        wget -O representative_octopamine_receptors.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=AFX62896.1,XP_022116353.1,XP_022123444.1,AEQ33589.1,AEO89318.1,AGV79326.1,AIC75370.1,AFG26689.1,BAD11157.1,NP_001171666.1,XP_004922133.2,NP_001091748.1,XP_048485348.1,XP_011568733.1,XP_048486473.1,NP_001034049.1,NP_651057.1,NP_001034043.2,NP_650754.2,NP_732541.1,NP_001280501.1,NP_001280505.1,XP_015839170.1,NP_001280520.1,XP_396348.4,XP_397139.3,XP_006557730.1,NP_001011565.1,CCO13925.1,XP_001122075.3,ASA47149.1,XP_022664697.1,NP_000675.1,NP_000015.2,NP_000016.1,NP_000671.2,NP_000670.1,NP_000669.1,NP_000672.3,NP_000674.2,NP_000673.2&rettype=fasta&retmode=text"
        grep -c ">" representative_octopamine_receptors.fasta # count and make sure that you've downloaded 33 OARs from the NCBI protein database
        """, language="bash")
        st.write("✔️extract the TyrOctR and OctB2R receptor sequences of CPB")
        st.code("""
        seqkit grep -p "TRINITY_DN219751_c0_g1_i1.p1" -p "TRINITY_DN22425_c0_g1_i3.p1" /media/raid/Wee/WeeYeZhi/output/transdecoder_results/transdecoder_ORF_extraction_results/good.modified_clean.trimmed.fasta.transdecoder.pep > /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/extracted_CPB_octopamine_receptor_sequences.fasta
        """, language="bash")
        st.write("✔️concatenate both the CPB OARs and its representative OARs together")
        st.code("""
        cat representative_octopamine_receptors.fasta extracted_CPB_octopamine_receptor_sequences.fasta > combined_CPB_and_representative_OARs_sequences.fasta
        """, language="bash")
        st.write("✔️rename the FASTA sequence headers")
        st.code("""
        sed 's/^>AFX62896.1.*/>AFX62896.1/; s/^>XP_022116353.1.*/>XP_022116353.1/; s/^>XP_022123444.1.*/>XP_022123444.1/; s/^>AEQ33589.1.*/>AEQ33589.1/; s/^>AEO89318.1.*/>AEO89318.1/; s/^>AGV79326.1.*/>AGV79326.1/; s/^>AIC75370.1.*/>AIC75370.1/; s/^>AFG26689.1.*/>AFG26689.1/; s/^>BAD11157.1.*/>BAD11157.1/; s/^>NP_001171666.1.*/>NP_001171666.1/; s/^>XP_004922133.2.*/>XP_004922133.2/; s/^>NP_001091748.1.*/>NP_001091748.1/; s/^>XP_048485348.1.*/>XP_048485348.1/; s/^>XP_011568733.1.*/>XP_011568733.1/; s/^>XP_048486473.1.*/>XP_048486473.1/; s/^>NP_001034049.1.*/>NP_001034049.1/; s/^>NP_651057.1.*/>NP_651057.1/; s/^>NP_001034043.2.*/>NP_001034043.2/; s/^>NP_650754.2.*/>NP_650754.2/; s/^>NP_732541.1.*/>NP_732541.1/; s/^>NP_001280501.1.*/>NP_001280501.1/; s/^>NP_001280505.1.*/>NP_001280505.1/; s/^>XP_015839170.1.*/>XP_015839170.1/; s/^>NP_001280520.1.*/>NP_001280520.1/; s/^>XP_396348.4.*/>XP_396348.4/; s/^>XP_397139.3.*/>XP_397139.3/; s/^>XP_006557730.1.*/>XP_006557730.1/; s/^>NP_001011565.1.*/>NP_001011565.1/; s/^>CCO13925.1.*/>CCO13925.1/; s/^>XP_001122075.3.*/>XP_001122075.3/; s/^>ASA47149.1.*/>ASA47149.1/; s/^>XP_022664697.1.*/>XP_022664697.1/; s/^>NP_000671.2.*/>NP_000671.2/; s/^>NP_000670.1.*/>NP_000670.1/; s/^>NP_000669.1.*/>NP_000669.1/; s/^>NP_000672.3.*/>NP_000672.3/; s/^>NP_000674.2.*/>NP_000674.2/; s/^>NP_000673.2.*/>NP_000673.2/; s/^>NP_000675.1.*/>NP_000675.1/; s/^>NP_000015.2.*/>NP_000015.2/; s/^>NP_000016.1.*/>NP_000016.1/' /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/combined_CPB_and_representative_OARs_sequences.fasta > /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/renamed_combined_CPB_and_representative_OARs_sequences.fasta
        """, language="bash")
        st.write("✔️upload the protein fasta file, 'renamed_combined_CPB_and_representative_OARs_sequences.fasta' into the DeepTMHMM 1.0 Web Server to predict the transmembrane domains of the all the 43 octopamine receptor sequences, to check whether they have 7 transmembrane domains of GPCR. Retain only those with 7 transmembrane helices (TMHs) for downstream analysis. Download the resulting .gff3 file from the DeepTMHMM Web Server and use the python script below to extract all the 7 TMs.**")
        # ----LOAD TM EXTRACTION PYTHON SCRIPT----
        # Check if the file exists before reading
        if extract_TM_regions.exists():
            with open(extract_TM_regions, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download TM Extraction Python Script",
                data=script_byte,
                file_name=extract_TM_regions.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{extract_TM_regions.name} does not exist.")
        st.write("✔️rename the FASTA sequence headers of the extracted TM1-TM7 regions")
        st.code("""
        sed 's/^>AFX62896.1_TM1-7.*/>PrTyrOctR(AFX62896.1)/; s/^>XP_022116353.1_TM1-7.*/>PrOctbeta2R(XP_022116353.1)/; s/^>XP_022123444.1_TM1-7.*/>PrOctbeta1R(XP_022123444.1)/; s/^>AEQ33589.1_TM1-7.*/>CsOctalpha1R(AEQ33589.1)/; s/^>AEO89318.1_TM1-7.*/>CsOctbeta2R(AEO89318.1)/; s/^>AGV79326.1_TM1-7.*/>CsOctbeta1R(AGV79326.1)/; s/^>AIC75370.1_TM1-7.*/>CsOctalpha2R(AIC75370.1)/; s/^>AFG26689.1_TM1-7.*/>CsTyrOctR(AFG26689.1)/; s/^>BAD11157.1_TM1-7.*/>BmTyrOctR(BAD11157.1)/; s/^>NP_001171666.1_TM1-7.*/>BmOctbeta2R(NP_001171666.1)/; s/^>XP_004922133.2_TM1-7.*/>BmOctbeta1R(XP_004922133.2)/; s/^>NP_001091748.1_TM1-7.*/>BmOctalpha1R(NP_001091748.1)/; s/^>XP_048485348.1_TM1-7.*/>PxTyrOctR(XP_048485348.1)/; s/^>XP_011568733.1_TM1-7.*/>PxOctbeta2R(XP_011568733.1)/; s/^>XP_048486473.1_TM1-7.*/>PxOAMB(XP_048486473.1)/; s/^>NP_001034049.1_TM1-7.*/>DmOctbeta2R(NP_001034049.1)/; s/^>NP_651057.1_TM1-7.*/>DmOctbeta1R(NP_651057.1)/; s/^>NP_001034043.2_TM1-7.*/>DmOctbeta3R(NP_001034043.2)/; s/^>NP_650754.2_TM1-7.*/>DmOctalpha2R(NP_650754.2)/; s/^>NP_732541.1_TM1-7.*/>DmOctalpha1R(NP_732541.1)/; s/^>NP_001280501.1_TM1-7.*/>TcOctbeta2R(NP_001280501.1)/; s/^>NP_001280505.1_TM1-7.*/>TcOctbeta3R(NP_001280505.1)/; s/^>XP_015839170.1_TM1-7.*/>TcOctalpha2R(XP_015839170.1)/; s/^>NP_001280520.1_TM1-7.*/>TcOctalpha1R(NP_001280520.1)/; s/^>XP_396348.4_TM1-7.*/>AmOctbeta2R(XP_396348.4)/; s/^>XP_397139.3_TM1-7.*/>AmOctbeta1R(XP_397139.3)/; s/^>XP_006557730.1_TM1-7.*/>AmOctbeta3R(XP_006557730.1)/; s/^>NP_001011565.1_TM1-7.*/>AmOctalpha1R(NP_001011565.1)/; s/^>CCO13925.1_TM1-7.*/>AmOctbeta4R(CCO13925.1)/; s/^>XP_001122075.3_TM1-7.*/>AmOctalpha2R(XP_001122075.3)/; s/^>ASA47149.1_TM1-7.*/>NiOctbeta2R(ASA47149.1)/; s/^>TRINITY_DN219751_c0_g1_i1.p1_TM1-7.*/>CcTyrOctR(TRINITY_DN219751_c0_g1_i1.p1)/; s/^>TRINITY_DN22425_c0_g1_i3.p1_TM1-7.*/>CcOctbeta2R(TRINITY_DN22425_c0_g1_i3.p1)/; s/^>XP_022664697.1_TM1-7.*/>VdOctbeta2R(XP_022664697.1)/; s/^>NP_000675.1_TM1-7.*/>HsOctbeta1AR(NP_000675.1)/; s/^>NP_000015.2_TM1-7.*/>HsOctbeta2AR(NP_000015.2)/; s/^>NP_000016.1_TM1-7.*/>HsOctbeta3AR(NP_000016.1)/; s/^>NP_000671.2_TM1-7.*/>HsOctalpha1AAR(NP_000671.2)/; s/^>NP_000670.1_TM1-7.*/>HsOctalpha1BAR(NP_000670.1)/; s/^>NP_000669.1_TM1-7.*/>HsOctalpha1DAR(NP_000669.1)/; s/^>NP_000672.3_TM1-7.*/>HsOctalpha2AAR(NP_000672.3)/; s/^>NP_000674.2_TM1-7.*/>HsOctalpha2CAR(NP_000674.2)/; s/^>NP_000673.2_TM1-7.*/>HsOctalpha2BAR(NP_000673.2)/;' /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/extracted_TM1-7_output.fasta > /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/renamed_extracted_TM1-7_output.fasta
        """, language="bash")
        st.write("✔️run MSA again to align all the extracted TM1-TM7 regions across the species to get the MSA alignment fasta file (.fasta). Input all the 7 TMs (MSA alignment fasta file) into the MEGA 11 software to construct the phylogenetic tree")
        st.code("""
        nohup muscle -align /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/renamed_extracted_TM1-7_output.fasta -output /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/MSA_extracted_TM1-7_output.fasta -threads 16 -consiters 2 -refineiters 100 > /media/raid/Wee/WeeYeZhi/output/phylogenetic_tree_construction_results/MSA_between_CPB_and_representative_OARs_results/MSA_extracted_TM1-7_output.log 2>&1 &
        """, language="bash")
        st.write("Below is kindly attached with the phylogenetic tree newick file (using the original tree, not bootstrap consensus tree) which is used as an input to reconstruct the phylogenetic tree using ape v5.8.1 and ggtree v3.16.3 via R")
        # ----LOAD THE ORIGINAL PHYLOGENETIC TREE NEWICK FILE----
        # Check if the file exists before reading
        if original_OctB2R_phylogenetic_tree.exists():
            with open(original_OctB2R_phylogenetic_tree, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Original Phylogenetic Tree Newick File",
                data=script_byte,
                file_name=original_OctB2R_phylogenetic_tree.name,  # Extract just the file name
                mime="text/plain",  # MIME type for plain text
            )
        else:
            st.error(f"{original_OctB2R_phylogenetic_tree.name} does not exist.")
        # ----LOAD THE PHYLOGENETIC TREE RECONSTRUCTION RSCRIPT----
        # Check if the file exists before reading
        if rscript_for_phylogenetic_tree_reconstruction.exists():
            with open(rscript_for_phylogenetic_tree_reconstruction, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Octopamine Receptors Phylogenetic Tree Reconstruction Rscript",
                data=script_byte,
                file_name=rscript_for_phylogenetic_tree_reconstruction.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{rscript_for_phylogenetic_tree_reconstruction.name} does not exist.")

        st.write("###")

        st.write("**8. Predict the 3D structure of the CPB octopamine beta2 receptor via AlphaFold3**")
        st.write("✔️install git via conda-forge")
        st.code("""
        conda create -n git
        conda install conda-forge::git
        conda activate git
        """, language="bash")
        st.write("✔️git clone the AlphaFold3 GitHub repository into your HPC system to obtain the AlphaFold3 source code and install AlphaFold3 locally")
        st.code("""
        git clone https://github.com/google-deepmind/alphafold3.git
        """, language="bash")
        st.write("✔️ensure that your HPC system has already installed wget and zstd for you to download & compress the required AlphaFold3 databases")
        st.code("""
        which wget zstd
        """, language="bash")
        st.write("✔️download all the required AlphaFold3 databases namely, BFD small, MGnify, PDB, PDB seqres, UniProt, UniRef90, NT, RFam, and RNACentral")
        st.code("""
        cd AlphaFold3
        chmod +x fetch_databases.sh # remember to make the script executable otherwise you are going to encounter the permission denied error
        ls -l fetch_databases.sh # double check whether the script is executable
        nohup ./fetch_databases.sh /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases > /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases/downloaded_alphafold3_public_databases.log 2>&1 &
        """, language = "bash")
        st.write("✔️build the AlphaFold3 Docker container with all the right python dependencies")
        st.code("""
        docker build -t alphafold3 -f docker/Dockerfile . # navigate to the working directory, /media/raid/Wee/WeeYeZhi/output/cloned_alphafold3_github_repository/alphafold3 first before you run this docker build command to build the AlphaFold3 Docker Image using the instructions written inside the Dockerfile
        """, language="bash")
        st.write("✔️display all the available flags of run_alphafold.py")
        st.code("""
        docker run alphafold3 python run_alphafold.py --helpfull
        """, language="bash")
        st.write("✔️ask your research institution/lab representative to help apply for the AlphaFold3 model parameters by filling in the application form attached below. Download the model parameters into the desired model_dir after you have successfully received it from Google Deepmind AlphaFold3")
        st.markdown("[Visit and Fill Up the AlphaFold3's Form to Apply for Model Parameters](https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform)")
        st.write("✔️unzip the alphafold3 model parameters file after you received it from Google")
        st.code("""
        zstd -d --keep af3.bin.zst
        """)
        st.write("✔️test your setup using the provided AlphaFold3 JSON file to check whether you have successfully set up AlphaFold3 on your HPC system")
        st.code("""
        docker run --user 0:0 -d --name alphafold3_testing -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/testing/fold_input.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output/testing
        """, language="bash")
        # ----LOAD ALPHAFOLD3 TESTING JSON FILE----
        # Check if the file exists before reading
        if test_alphafold3.exists():
            with open(test_alphafold3, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download AlphaFold3 Testing JSON File",
                data=script_byte,
                file_name=test_alphafold3.name,  # Extract just the file name
                mime="application/json",  # MIME type for shell scripts
            )
        else:
            st.error(f"{test_alphafold3.name} does not exist.")
        st.write("✔️prepare the input JSON file for the CPB octopamine beta 2 receptor and run AlphaFold3 to predict the 3D receptor structure")
        st.code("""
        docker run --user 0:0 -d --name alphafold3_prediction_of_CPB_octopamine_beta2_receptor -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/CPB_octopamine_beta2_receptor/CPB_OctB2R.json --model_dir=/root/models --output_dir=/root/af_output/CPB_octopamine_beta2_receptor
        """, language="bash")
        # ----LOAD THE CPB OCTOPAMINE BETA 2 RECEPTOR JSON FILE----
        # Check if the file exists before reading
        if CPB_OctB2R.exists():
            with open(CPB_OctB2R, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download CPB_OctB2R JSON File",
                data=script_byte,
                file_name=CPB_OctB2R.name,  # Extract just the file name
                mime="application/json",  # MIME type for shell scripts
            )
        else:
            st.error(f"{CPB_OctB2R.name} does not exist.")
        st.markdown("[Visit the AlphaFold3's Official GitHub Page](https://github.com/google-deepmind/alphafold3)")
        st.markdown("[Read how to install and run your first AlphaFold3 prediction (tutorial)](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)")
        st.markdown("[Visit the AlphaFold3's Introduction Page](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/using-the-alphafold-3-source-code/)")
        st.markdown("[Visit the AlphaFold3's Model Parameters Page](https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#obtaining-model-parameters)")
        st.markdown("[Test your setup using the AlphaFold JSON file to check whether you have successfully set up AlphaFold3 on your HPC system](https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#installation-and-running-your-first-prediction)")
        st.markdown("[Check how to create an input JSON file to run AlphaFold3 3D protein structure prediction](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)")

        st.write("###")

        st.write("**9. Prepare the list of ligands using Gypsum-DL**")
        st.write("✔️install the latest version of Gypsum-DL via the Python Package Index (PyPI) website using the following command")
        st.code("""
        conda create -n gypsum-dl
        conda activate gypsum-dl
        conda install conda-forge::pip
        pip install gypsum-dl # pillow v12.2.0, dimorphite-dl v2.0.2, gypsum-dl v1.3.0, loguru v0.7.3, molvs v0.1.1, mpi4py v4.1.1, numpy v2.4.4, rdkit v2025.9.6, scipy v1.17.1, & six v1.17.0
        gypsum-dl --help
        """, language="bash")
        st.write("✔️install the latest version of OpenBabel via CMake compilation (compile OpenBabel on your own)")
        st.code("""
        conda create -n obabel
        conda activate obabel
        conda install conda-forge::cmake # install cmake
        cd /media/raid/Wee/WeeYeZhi/output/cloned_obabel_github_source_code
        git clone https://github.com/openbabel/openbabel.git # clone the Open Babel source code from GitHub
        nano /media/raid/Wee/WeeYeZhi/output/cloned_obabel_github_source_code/openbabel/include/openbabel/obutil.h # edit the obutil.h file by adding the line, #include <ctime>, below the line, #include <math.h>, and save it. This is important to ensure the compilation process run smoothly and successfully
        cd openbabel # navigate to the Open Babel source code directory
        git checkout tags/openbabel-3-1-1 # checkout the version 3.1.1 release tag
        mkdir build # create a build directory
        cd build
        cmake .. -DCMAKE_INSTALL_PREFIX=/media/raid/Wee/WeeYeZhi/output/obabel_installation -DCMAKE_BUILD_TYPE=Release # configure Open Babel for building in your build directory with your chosen install prefix
        make -j16 # compile the source code with 16 HPC cores
        make install # install open babel
        export PATH=/media/raid/Wee/WeeYeZhi/output/obabel_installation/bin:$PATH # add openbabel to your PATH to run it easily after succesfully installing OpenBabel. If you do not export into PATH, you have run the command, /media/raid/Wee/WeeYeZhi/output/obabel_installation/bin/obabel every single time when you want to run OpenBabel
        obabel -H # output the usage information to test whether open babel has been successfully installed
        obabel -V # output its version number
        rm -rf openbabel # you can delete all the open babel source code files after you have successfully installed Open Babel to your designated directory to save space
        """)
        st.write("✔️prepare the input file containing all the 18 SMILES strings of ligands with the PubChem ID, 26451, 18526103, 76145148, 19606232, 1480785, 26752, 577782, 255273, 1001, 36324, 2726, 2803, 36326, 439570, 4184, 4436, 5775, and 8969 from the PubChem database")
        # ----LOAD THE LIGAND SMILES STRING FILE----
        # Check if the file exists before reading
        if ligands.exists():
            with open(ligands, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Ligand SMILES String File (.smi)",
                data=script_byte,
                file_name=ligands.name,  # Extract just the file name
                mime="text/plain",  # MIME type for plain text
            )
        else:
            st.error(f"{ligands.name} does not exist.")
        st.write("✔️run Gypsum-DL to generate the 3D molecular structures of ligands from SMILES format to SDF and PDB. Check all potential ionization, tautomeric, chiral, cis/trans isomeric, and ring-conformational forms at a pH range between 6.5 and 7.5 (The whole Gypsum-DL pipeline desalt all molecules, ionize all molecules, generate tautomers for all molecules, apply Durrant-lab filters to all molecules, enumerate all possible enantiomers for all molecules, enumerate all possible cis-trans isomers for all molecules, convert all molecules to 3D structures, generate several conformers of molecules with non-aromatic rings (boat vs chair), minimize all 3D molecular structures, make PDB output files, and save all molecules)")
        st.code("""
        nohup gypsum-dl --source /media/raid/Wee/WeeYeZhi/output/ligand_preparation_results/ligands.smi --output_folder /media/raid/Wee/WeeYeZhi/output/ligand_preparation_results/gypsumdl_results --job_manager multiprocessing --num_processors 16 --max_variants_per_compound 5 --thoroughness 3 --min_ph 6.5 --max_ph 7.5 --pka_precision 1 --use_durrant_lab_filters --separate_output_files >  /media/raid/Wee/WeeYeZhi/output/ligand_preparation_results/gypsumdl_results/gypsumdl_results.log 2>&1 & # dont use the flag, --add_pdb_output, because the bond order information of ligands get lost after converting the ligands to PDB format, so you better omit this flag and use the generated default SDF files and later use OpenBabel to split and convert the SDF files into MOL2 format. 
        """, language="python")
        st.write("✔️split each SDF file of ligand (generated by Gypsum-DL) into multiple separate files of ligand variants (due to having different tautomers, ionization states, & cis-trans isomers) via OpenBabel")
        st.code("""
        obabel ../gypsumdl_results/1-_4-chlorophenyl_imidazolidin-2-one__input1.sdf -O PubChem_26451_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-bromophenyl_imidazolidin-2-one__input2.sdf -O PubChem_18526103_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-ethoxyphenyl_imidazolidin-2-one__input3.sdf -O PubChem_76145148_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-ethylphenyl_imidazolidin-2-one__input4.sdf -O PubChem_19606232_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-fluorophenyl_imidazolidin-2-one__input5.sdf -O PubChem_1480785_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-methoxyphenyl_imidazolidin-2-one__input6.sdf -O PubChem_26752_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-methylphenyl_imidazolidin-2-one__input7.sdf -O PubChem_577782_variant.sdf -m 
        obabel ../gypsumdl_results/1-phenylimidazolidin-2-one__input8.sdf -O PubChem_255273_variant.sdf -m 
        obabel ../gypsumdl_results/2-phenylethylamine__input9.sdf -O PubChem_1001_variant.sdf -m 
        obabel ../gypsumdl_results/amitraz__input10.sdf -O PubChem_36324_variant.sdf -m 
        obabel ../gypsumdl_results/chlorpromazine__input11.sdf -O PubChem_2726_variant.sdf -m 
        obabel ../gypsumdl_results/clonidine__input12.sdf -O PubChem_2803_variant.sdf -m 
        obabel ../gypsumdl_results/DPMF__input13.sdf -O PubChem_36326_variant.sdf -m 
        obabel ../gypsumdl_results/L-carvone__input14.sdf -O PubChem_439570_variant.sdf -m 
        obabel ../gypsumdl_results/mianserin__input15.sdf -O PubChem_4184_variant.sdf -m 
        obabel ../gypsumdl_results/naphazoline__input16.sdf -O PubChem_4436_variant.sdf -m 
        obabel ../gypsumdl_results/phentolamine__input17.sdf -O PubChem_5775_variant.sdf -m 
        obabel ../gypsumdl_results/yohimbine__input18.sdf -O PubChem_8969_variant.sdf -m 
        """, language="bash")
        st.write("✔️convert all the 84 SDF files of ligands from the .sdf format to .mol2 format via OpenBabel")
        st.code("""
        obabel *.sdf -omol2 -m
        """, language="bash")
        st.write("✔️perform batch preparation of ligands using the prepare_ligand4.py package of AutoDock MGL Tools to convert all the 84 ligands from PDB format to PDBQT format")
        st.code("""
        tar -xvf mgltools_x86_64Linux2_1.5.7p1.tar.gz # download the mgltools_x86_64Linux2_1.5.7p1.tar.gz file from https://ccsb.scripps.edu/mgltools/downloads/ and unzip it to install AutoDock MGL Tools via Linux
        tar -xvf MGLToolsPckgs.tar.gz # look for MGLToolsPckgs.tar.gz and unzip it
        navigate to autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py
        navigate to autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py # just in case if you want to use it to prepare protein in batch later
        """, language="bash")
        # ----LOAD AUTODOCK4 prepare_ligand4.py SCRIPT----
        # Check if the file exists before reading
        if autodock4_ligand_batch_preparation_script.exists():
            with open(autodock4_ligand_batch_preparation_script, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download AutoDock4 prepare_ligand4.py script",
                data=script_byte,
                file_name=autodock4_ligand_batch_preparation_script.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{autodock4_ligand_batch_preparation_script.name} does not exist.")
        # ----LOAD THE LIGAND BATCH PREPARATION BASH SCRIPT----
        # Check if the file exists before reading
        if ligand_batch_preparation_script.exists():
            with open(ligand_batch_preparation_script, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Ligand Batch Preparation Bash Script",
                data=script_byte,
                file_name=ligand_batch_preparation_script.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{ligand_batch_preparation_script.name} does not exist.")
        st.markdown("[Visit Gypsum-DL official PyPi Installation Page](https://pypi.org/project/gypsum-dl/)")
        st.markdown("[Visit Gypsum-DL official GitHub Page](https://github.com/durrantlab/gypsum_dl)")
        st.markdown("[Visit OpenBabel official Tutorial Page](https://openbabel.org/docs/Command-line_tools/babel.html)")
        st.markdown("[Read how to install OpenBabel in Linux Ubuntu](https://stackoverflow.com/questions/75814835/install-open-babel-in-ubuntu#:~:text=1%20Answer,sudo%20make%20install)")
        st.markdown("[Read how to perform batch preparation of ligands for docking via AutoDock Vina](https://www.researchgate.net/post/Batch_Ligand_Preparation_on_Autodock_Vina)")
        st.markdown("[Visit the official AutoDock website](https://autodock.scripps.edu/)")
        st.markdown("[Read the useful resources of AutoDock4 and AutoDock Vina](https://autodock.scripps.edu/resources/)")
        st.markdown("[Watch YouTube video to learn how to perform batch preparation of ligands](https://www.youtube.com/watch?v=_Blz2DxSAtQ&t=44s)")

        st.write("###")

        st.markdown("[Read how to install and run AutoDock Vina](https://autodock-vina.readthedocs.io/en/latest/installation.html)")

        st.write("Gibbs free energy of a chemical reaction is normally calculated to check whether the chemical reaction can take place naturally (whether the reaction is considered spontaneous). If the Gibbs free energy is negative, then it means that the chemical reaction is considered spontaneous and the chemical reaction is considered exergonic (release free energy). The second law of thermodynamics states that the entropy of the universe is always positive. Entropy is the measure of disorder of the system, where increasing the volume, temperature (Kelvin), and number of moles of particles will increase the disorderness/randomness (entropy) of the system. Entropy is the nature's way to spread/disperse energy throughout the space to achieve an energitically favourable system")
        st.write("Try out the GROMICA web app to analyze the RMSD, RMSF, Hbond, SASA, and Rg of MD simulation with zero coding to produce publication-quality images")
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