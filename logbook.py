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
trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa = current_dir / "assets" / "trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa.txt"
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
CPBOctB2R_octopamine_complex = current_dir / "assets" / "CPBOctB2R_octopamine_complex.json"
ligands = current_dir / "assets" / "ligands.smi"
ligand_batch_preparation_script = current_dir / "assets" / "batch_preparation_of_ligands.sh"
autodock4_ligand_batch_preparation_script = current_dir / "assets" / "prepare_ligand4.py"
galaxyrefine_usage_information_file = current_dir / "assets" / "galaxyrefine_README.md"
qmeanbrane_manual = current_dir / "assets" / "qmeanbrane_interpretation_manual.txt"
prosa2003_manual = current_dir / "assets" / "prosa2003_manual.pdf"
gpcritasser_cscore = current_dir / "assets" / "cscore_of_GPCR-I-TASSER-predicted_CcOctB2R.txt"
alphafold2_jupyternotebook = current_dir / "assets" / "AlphaFold2.ipynb"
ENA_manifest = current_dir / "assets" / "manifest.txt"
CPB_transcriptome_assembly = current_dir / "assets" / "CPB_transcriptome_assembly_v2.fasta.gz"

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
        options=["Phase 1: Sequence-Based Analysis", "Phase 2: Transcriptomic and Structural-Based Analysis", "Phase 3: Molecular Docking and Dynamics Simulation", "Submit Transcriptome Assembly to ENA", "Additional Notes"],
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
        st.write("✔️run Trinity to perform comprehensive transcriptome assembly to assemble all the 9 raw transcriptomic Illumina paired end reads (sequenced by LKM and INBIOSIS) except the pupal sample")
        st.code("""
        nohup Trinity --seqType fq --max_memory 450G --samples_file /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS_except_pupa/trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa.txt --CPU 16 --normalize_by_read_set --SS_lib_type RF --full_cleanup --output /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS_except_pupa > /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS_except_pupa/trinity_output.log 2>&1 & 
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
        # ----LOAD TRINITY SAMPLE FILE sequenced BY LKM & INBIOSIS except pupa----
        # Check if the file exists before reading
        if trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa.exists():
            with open(trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Trinity Sample File Sequenced by LKM & INBIOSIS except pupa",
                data=script_byte,
                file_name=trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{trinity_samples_sequenced_by_LKM_and_INBIOSIS_except_pupa.name} does not exist.")
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
        nohup cd-hit-est -i /media/raid/Wee/WeeYeZhi/output/trinity_results/trinity_assembly_done_by_LKM_and_INBIOSIS_except_pupa/trinity_assembly_done_by_LKM_and_INBIOSIS_except_pupa.Trinity.fasta -o /media/raid/Wee/WeeYeZhi/output/cd_hit_est_results/cd_hit_est_after_trimmomatic_without_pupa -c 0.95 -n 5 -T 4 > cd_hit_est_after_trimmomatic_without_pupa_output.log 2>&1 & # specify fewer number of threads (use -T 4 in this case) if possible to avoid encountering OutofMemory error
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
                mime="application/gzip",  # MIME type for shell scripts
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
        st.write("✔️run MSA via MUSCLE and view the MSA results (.aln) via JalView")
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

        st.write("**8. Predict the 3D structure of the CPB octopamine beta2 receptor via AlphaFold3 (run AlphaFold3 locally)**")
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
        st.write("✔️check the installed AlphaFold3 version")
        st.code("""
        docker run --rm -it alphafold3 bash
        git describe --tags # The version of the AlphaFold3 that you are using is AlphaFold3 v3.0.1-127-g608edb6
        """, language="bash")
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
        st.write("✔️prepare the input JSON file for the CPB octopamine beta 2 receptor and other insect octopamine beta 2 receptors (Plutella xylostella, Pieris rapae, Bombyx mori). Run AlphaFold3 to predict the 3D receptor structure (use 5 different random seeds and perform a multi-seed prediction to predict 5 independent AlphaFold3 structures and choose the best AlphaFold3 structure with the highest ranking score for downstream docking or perform an ensemble docking to dock a ligand to five independent AlphaFold3-predicted structures. Compute the RMSD values between the best AlphaFold3-predicted structure with the highest ranking score and other 4 best model seed structures, if their RMSD values are very similar to each other, you can just choose the best AlphaFold3 structure with the highest ranking score for downstream docking. Check whether the key binding residues are aligned and whether the side chains are stable)")
        st.code("""
        docker run --user 0:0 -d --name alphafold3_prediction_of_CPB_octopamine_beta2_receptor -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/CPB_octopamine_beta2_receptor/CPB_OctB2R.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output/CPB_octopamine_beta2_receptor
        docker run --user 0:0 -d --name alphafold3_prediction_of_Px_octopamine_beta2_receptor -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/Px_octopamine_beta2_receptor/Px_OctB2R.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output/Px_octopamine_beta2_receptor
        docker run --user 0:0 -d --name alphafold3_prediction_of_Pr_octopamine_beta2_receptor -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/Pr_octopamine_beta2_receptor/Pr_OctB2R.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output/Pr_octopamine_beta2_receptor
        docker run --user 0:0 -d --name alphafold3_prediction_of_Bm_octopamine_beta2_receptor -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/Bm_octopamine_beta2_receptor/Bm_OctB2R.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output/Bm_octopamine_beta2_receptor
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

        st.write("**9. Predict the 3D structure of the CPB octopamine beta 2 receptor via Swiss Model webserver**")
        st.write("✔️use G3M4F8.1.A as the template to perform homology modelling of the CPB octopamine beta 2 receptor since Swiss Model webserver itself identified G3M4F8.1.A has the highest coverage of 0.9795, GMQE score of 0.8396, sequence identity percentage of 86.91")
        st.markdown("[Visit Swiss Model webserver](https://swissmodel.expasy.org/interactive)")

        st.write("###")

        st.write("**10. Predict the 3D structure of the CPB octopamine beta 2 receptor via GPCR-I-TASSER webserver**")
        # ----LOAD GCPR-I-TASSER results interpretation manual----
        # Check if the file exists before reading
        if gpcritasser_cscore.exists():
            with open(gpcritasser_cscore, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download GPCR-I-TASSER results interpretation manual to interpret GPCR-I-TASSER of CcOctB2R results",
                data=script_byte,
                file_name=gpcritasser_cscore.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{gpcritasser_cscore.name} does not exist.")
        st.markdown("[Visit GPCR-I-TASSER webserver](https://aideepmed.com/GPCR-I-TASSER/)")
        st.markdown("[Read how to interpret the GPCR-I-TASSER results to know how to explain the C-score, TM-score, RMSD, No of decoys, and cluster density of the predicted structure](https://aideepmed.com/GPCR-I-TASSER/output/G1856/cscore.txt)")

        st.write("###")

        st.write("**11. Predict the 3D structure of the CPB octopamine beta 2 receptor via AlphaFold2**")
        # ----LOAD AlphaFold2 results interpretation manual----
        # Check if the file exists before reading
        if alphafold2_jupyternotebook.exists():
            with open(alphafold2_jupyternotebook, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download AlphaFold2 Jupyter Notebook to run AlphaFold2 prediction of CPB OctB2R",
                data=script_byte,
                file_name=alphafold2_jupyternotebook.name,  # Extract just the file name
                mime="application/x-ipynb+json",  # MIME type for shell scripts
            )
        else:
            st.error(f"{alphafold2_jupyternotebook.name} does not exist.")

        st.write("###")

        st.write("**12. Color the best ranked AlphaFold3-predicted CPB octopamine beta 2 receptor by pLDDT score and save it using UCSF ChimeraX**")
        st.code("""
        cd C:\\Users\\USER\\Documents\\Master_of_Science_in_Systems_Biology_GRA\\Experimental_results\\alphafold3_prediction_results\\output\\CPB_octopamine_beta2_receptor\\five_top_ranked_AlphaFold3_models # change to this working directory first in UCSF ChimeraX
        open CPB_OctB2R_seed-1_sample-2_model.cif # load the AlphaFold3-predicted CPBOctB2R structure in .cif format. Remember to load the structure in .cif format, not the .pdb format that you converted using OpenBabel, otherwise the command won't work
        color bfactor palette alphafold # color the AlphaFold3-predicted CPBOctB2R structure based on the AlphaFold3 pLDDT score
        """, language="bash")

        st.write("###")

        st.write("**13. Prepare the list of ligands using Gypsum-DL**")
        st.write("✔️filter the ligands based on their insecticide-likeness properties by submitting SDF files of ligands to the APPi web server via http://pesticides.cau.edu.cn/APPi.")
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
        obabel -V # output its version number (you are using OpenBabel v3.1.0)
        rm -rf openbabel # you can delete all the open babel source code files after you have successfully installed Open Babel to your designated directory to save space
        """)
        st.write("✔️prepare the input file containing all the 18 SMILES strings of ligands with the PubChem ID, 4581, 36324, 36326, 26451, 18526103, 76145148, 19606232, 1480785, 26752, 577782, 255273, 2726, 4184, 5775, and 8969 from the PubChem database")
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
        nohup gypsum-dl --source /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/ligands.smi --output_folder /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/gypsumdl_results --job_manager multiprocessing --num_processors 16 --max_variants_per_compound 5 --thoroughness 3 --min_ph 6.5 --max_ph 7.5 --pka_precision 1 --use_durrant_lab_filters --separate_output_files >  /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/gypsumdl_results/gypsumdl_results.log 2>&1 & # dont use the flag, --add_pdb_output, because the bond order information of ligands get lost after converting the ligands to PDB format, so you better omit this flag and use the generated default SDF files and later use OpenBabel to split and convert the SDF files into MOL2 format. 
        """, language="python")
        st.write("✔️split each SDF file of ligand (generated by Gypsum-DL) into multiple separate files of ligand variants (due to having different tautomers, ionization states, & cis-trans isomers) via OpenBabel")
        st.code("""
        obabel ../gypsumdl_results/1-_4-chlorophenyl_imidazolidin-2-one__input1.sdf -O PubChem_26451_variant.sdf -m # navigate the working directory, /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/openbabel_results, before you run all these list of obabel commands for splitting.
        obabel ../gypsumdl_results/1-_4-bromophenyl_imidazolidin-2-one__input2.sdf -O PubChem_18526103_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-ethoxyphenyl_imidazolidin-2-one__input3.sdf -O PubChem_76145148_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-ethylphenyl_imidazolidin-2-one__input4.sdf -O PubChem_19606232_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-fluorophenyl_imidazolidin-2-one__input5.sdf -O PubChem_1480785_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-methoxyphenyl_imidazolidin-2-one__input6.sdf -O PubChem_26752_variant.sdf -m 
        obabel ../gypsumdl_results/1-_4-methylphenyl_imidazolidin-2-one__input7.sdf -O PubChem_577782_variant.sdf -m 
        obabel ../gypsumdl_results/1-phenylimidazolidin-2-one__input8.sdf -O PubChem_255273_variant.sdf -m 
        obabel ../gypsumdl_results/amitraz__input9.sdf -O PubChem_36324_variant.sdf -m 
        obabel ../gypsumdl_results/chlorpromazine__input10.sdf -O PubChem_2726_variant.sdf -m 
        obabel ../gypsumdl_results/DPMF__input11.sdf -O PubChem_36326_variant.sdf -m 
        obabel ../gypsumdl_results/mianserin__input12.sdf -O PubChem_4184_variant.sdf -m 
        obabel ../gypsumdl_results/phentolamine__input13.sdf -O PubChem_5775_variant.sdf -m 
        obabel ../gypsumdl_results/yohimbine__input14.sdf -O PubChem_8969_variant.sdf -m 
        obabel ../gypsumdl_results/octopamine__input15.sdf -O PubChem_4581_variant.sdf -m 
        """, language="bash")
        st.write("✔️convert all the 73 SDF files of ligands from the .sdf format to .mol2 format via OpenBabel")
        st.code("""
        ls -l | wc -l # count the total number of SDF files (variants) generated after splitting 
        obabel *.sdf -omol2 -m
        """, language="bash")
        st.write("✔️perform batch preparation of ligands using the prepare_ligand4.py package of AutoDock MGL Tools to convert all the 84 ligands from PDB format to PDBQT format")
        st.code("""
        tar -xvf mgltools_x86_64Linux2_1.5.7p1.tar.gz # download the mgltools_x86_64Linux2_1.5.7p1.tar.gz file from https://ccsb.scripps.edu/mgltools/downloads/ and unzip it to install AutoDock MGL Tools via Linux
        tar -xvf MGLToolsPckgs.tar.gz # look for MGLToolsPckgs.tar.gz and unzip it
        navigate to autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py
        navigate to autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py # just in case if you want to use it to prepare protein in batch later
        /media/raid/Wee/WeeYeZhi/output/autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7/bin/pythonsh /media/raid/Wee/WeeYeZhi/output/autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -h # output the command-line usage. remember to use the MGLTools internal Python interpreter (pythonsh) instead of using the HPC system's python when you run the prepare_ligand4.py python script
        bash /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/ligand_batch_preparation_script/batch_preparation_of_ligands.sh # run this bash script to perform batch preparation of ligands by adding polar hydrogens,adding Gasteiger charges, and removing non-polar hydrogens from the ligands. Remember to run this command in the directory containing a list of Avogadro2 energy-minimised ligands (.mol2) files (you have to change your working directory to the directory containing a list of Avogadro2 energy-minimised ligands (.mol2) files first)
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
        st.markdown("[Read the command-line usage manual of the prepare_ligand4.py script](https://github.com/lmdu/AutoDockTools_py3/blob/master/AutoDockTools/Utilities24/prepare_ligand4.py)")
        st.markdown("[Visit the official AutoDock website](https://autodock.scripps.edu/)")
        st.markdown("[Read the useful resources of AutoDock4 and AutoDock Vina](https://autodock.scripps.edu/resources/)")
        st.markdown("[Watch YouTube video to learn how to perform batch preparation of ligands](https://www.youtube.com/watch?v=_Blz2DxSAtQ&t=44s)")
        st.write("Note: Unlike Open Babel, Gypsum-DL often generates multiple variants of each compound with differing ionization, tautomeric, chiral, cis/trans isomeric, and ring-conformational states.")

        st.write("###")

        st.write("**14. Filter the list of ligands based on ZINC22 filters, 3D compounds, molecular weight <= 450Da, logP<=5, only clean compound, only in stock purchase, ref and mid pH, and default charge value**")
        st.code("""
        curl -X GET https://cartblanche22.docking.org/smiles.txt -F smiles=@/media/raid/Wee/WeeYeZhi/output/ligand_preparation_results/ZINC22_filter_results/one_ligand_smiles.txt -F dist=0 -F adist=0 -o /media/raid/Wee/WeeYeZhi/output/ligand_preparation_results/ZINC22_filter_results/ZINC22_filter_results.json
        curl -X GET https://cartblanche22.docking.org/search/saveResult/22059c83-86db-4f0b-89af-b4acb89afb65.txt -o zinc_results.txt
        curl -L "https://cartblanche.docking.org/search/saveResult/916a2d21-1286-49de-9d7d-29100dfc1f8c.txt" -o zinc_results.txt
        """, language="bash")
        st.markdown("[Visit ZINC22 SMILES string search database](https://cartblanche22.docking.org/search/smiles)")
        st.markdown("[Read how to search the ZINC22 database by ZINCIDs, SMILES string, supplier code IDs, and download molecules in bulk using the curl command](https://wiki.docking.org/index.php/Zinc22:Searching)")

        st.write("###")

        st.write("**15. Run TMHMM to predict the transmembrane domains of the inactivate state of the human B2-adrenergic GPCR sequence with PDB ID 2RH1**")
        st.code("""
        cd /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin # navigate into this working directory before you run tmhmm
        export PATH="/media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin:$PATH" # note that the export of the PATH environment variable to the tmhmm executable file is just temporary and will only take effect for the current terminal session
        which tmhmm # double check the system can access the tmhmm executable file correctly after exporting the PATH variable 
        tmhmm /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/rcsb_pdb_2RH1.fasta > /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/tmhmm_results/tmhmm_2RH1.out
        """, language="bash")

        st.write("###")

        st.write("**16. Perform structural superimposition between the human B2-adrenergic receptor (2RH1) and the AlphaFold3-predicted CPB octopamine receptor (seed1sample2) using UCSF ChimeraX**")
        st.write("✔️install gpcrmining python package via pip and obtain the class A GPCR human beta 2 adrenergic receptor sequence to get the mapping between the GPCR generic residues and UniProt numbering system residues")
        st.code("""
        conda create -n gpcrmining
        conda activate gpcrmining
        conda install conda-forge::pip
        pip install gpcrmining # successfully installed certifi-2026.2.25 charset_normalizer-3.4.7 click-8.3.2 gpcrmining-0.2.0 idna-3.11 numpy-2.4.4 pandas-3.0.2 python-dateutil-2.9.0.post0 requests-2.33.1 six-1.17.0 urllib3-2.6.3
        python -m gpcrmining.gpcrdb -n adrb2_human -d /media/raid/Wee/WeeYeZhi/output/gpcrmining_results/2RH1
        """, language="bash")
        st.write("✔️perform structural superimposition between ADRB2_HUMAN (P07550) and 2RH1 to identify the superimposed ligand binding residues of 2RH1 since 2RH1 is a type of ADRB2_HUMAN receptor (human beta 2 adrenergic receptor)")
        st.code("""
        pwd # print the current working directory of UCSF ChimeraX to check your file location now
        cd C:\\Users\\USER\\Documents\\Master_of_Science_in_Systems_Biology_GRA\\Experimental_results\\ligand_binding_site_identification_results # navigate into the desired working directory containing the predicted 3D protein structures
        open AF-P07550-F1-model_v6.pdb # load the ADRB2_HUMAN (P07550) AlphaFold structure which uses the UniProt numbering system that follow the GPCRdb UniProt numbering system converted by gpcrmining python package
        open 2RH1.pdb # load the 2RH1 structure into the UCSF ChimeraX interface and it will become your first loaded structure (labelled as #1)
        select #2; hide sel target a # select your second loaded protein structure (#1) and hide all of its atom and stick representation
        select clear # clear your current selection
        matchmaker #2 to #1 # perform structural superimposition to structurally superimpose 2RH1 (labelled as #2) to the human beta 2 adrenergic receptor (ADRB2_HUMAN) (P07550) as the reference structure (labelled as #1) to identify the ligand binding residues of 2RH1 that superimpose with the residues of P07550 (W109,D113,V114,V117,F193,Y199,S203,S204,S207,W286,F289,F290,N293,Y308,N312,Y316). In this case, #1 is the reference structure.
        select #1/A:109,113,114,117,193,199,203,204,207,286,289,290,293,308,312,316; show sel target ab # select the ligand binding residues of ADRB2_HUMAN (P07550) and show them in atom and bond (stick) representation
        manually hover your mouse over each of the ligand binding residue of ADRB2_HUMAN (P07550) and record down the superimposed ligand binding residue of 2RH1
        select #2/A:116,120,121,124,200,206,210,211,214,421,424,425,428,443,447,451; show sel target ab # select the ligand binding residues of 2RH1 and show them in atom and bond (stick) representation
        select #1/A; cartoon hide sel # hide the cartoon representation of structure #1
        select #2/A; cartoon hide sel # hide the cartoon representation of structure #2
        select #1/A:109,113,114,117,193,199,203,204,207,286,289,290,293,308,312,316; label sel text "{0.name} {0.number}{0.insertion_code}" # label the amino acid residue with three-letter codes
        select #2/A:116,120,121,124,200,206,210,211,214,421,424,425,428,443,447,451; label sel text "{0.name} {0.number}{0.insertion_code}"
        ui mousemode right "move label" # change the mouse setting so that you can right click the amino acid residue label and move them to arrange properly
        save "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/ligand_binding_site_identification_results/superimposition between ligand binding residues of ADRB2_HUMAN and 2RH1.cxs"
        """, language="bash")
        st.write("✔️perform structural superimposition between 2RH1 and AlphaFold3-predicted CPB octopamine beta 2 receptor to identify the ligand binding residues of the AlphaFold3-predicted CPB octopamine beta 2 receptor")
        st.code("""
        pwd # print the current working directory of UCSF ChimeraX to check your file location now
        cd C:\\Users\\USER\\Documents\\Master_of_Science_in_Systems_Biology_GRA\\Experimental_results\\ligand_binding_site_identification_results # navigate into the desired working directory containing the predicted 3D protein structures
        open 2RH1.pdb (This will become your #1 structure)
        open CPB_OctB2R_seed-1_sample-2_model.pdb (This will become your #2 structure)
        select #1; hide sel target a # select your first loaded protein structure (#1) and hide all of its atom and stick representation
        matchmaker #2 to #1 # perform structural superimposition to structurally superimpose the AlphaFold3-predicted CPB octopamine beta 2 receptor (labelled as #2) to 2RH1 as the reference structure (labelled as #1) to identify the ligand binding residues of AlphaFold3-predicted CPB octopamine beta 2 receptor that superimpose with the residues of 2RH1 (W116,D120,D121,V124,F200,Y206,S210,S211,S214,W421,F424,F425,N428,Y443,N447,Y451). In this case, #1 is the reference structure.
        select #1/A:116,120,121,124,200,206,210,211,214,421,424,425,428,443,447,451; show sel target ab # select the ligand binding residues of 2RH1 and show them in atom and bond (stick) representation
        manually hover your mouse over each of the ligand binding residue of 2RH1 and record down the superimposed ligand binding residue of the AlphaFold3-predicted CPB octopamine beta 2 receptor
        select #2/A:105,109,110,113,188,194,198,199,202,290,293,294,297,313,317,321; show sel target ab # select the ligand binding residues of CPB OctB2R and show them in atom and bond (stick) representation
        select #1/A; cartoon hide sel # hide the cartoon representation of structure #1
        select #2/A; cartoon hide sel # hide the cartoon representation of structure #2
        select #1/A:116,120,121,124,200,206,210,211,214,421,424,425,428,443,447,451; label sel text "{0.name} {0.number}{0.insertion_code}"
        select #2/A:105,109,110,113,188,194,198,199,202,290,293,294,297,313,317,321; label sel text "{0.name} {0.number}{0.insertion_code}"
        ui mousemode right "move label" # change the mouse setting so that you can right click the amino acid residue label and move them to arrange properly
        save "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/ligand_binding_site_identification_results/superimposition between ligand binding residues of 2RH1 and CPB OctB2R.cxs"
        """, language="bash")
        st.write("✔️perform pairwise structural superimposition between the ligand binding site of refined, energy-minimized AlphaFold3-predicted CPB octopamine beta 2 receptor (seed1sample2) and the ligand binding site of other refined, energy-minimized AlphaFold3-predicted PxOctB2R (Px_OctB2R_seed-3_sample-3_model.pdb), AlphaFold3-predicted PrOctB2R (Pr_OctB2R_seed-1_sample-3_model.pdb), AlphaFold3-predicted BmOctB2R (Bm_OctB2R_seed-1_sample-2_model.pdb), and CsOctB2R (AF-G3M4F8-F1-model_v6.pdb) (already deposited in UniProt database) respectively to investigate whether the AlphaFold3-predicted CPB octopamine beta 2 receptor structure is similar in other insects with the same taxonomy Lepidopteran order")
        st.code("""
        pwd # print the current working directory of UCSF ChimeraX to check your file location now
        cd C:\\Users\\USER\\Documents\\Master_of_Science_in_Systems_Biology_GRA\\Experimental_results\\yasara_energy_minimization_results # navigate into the desired working directory containing the predicted 3D protein structures
        open alphafold3_prediction\\energy_minimized_alphafold3_predicted_CPB_octopamine_beta2_receptor.pdb # this will become your #1 structure
        open PxOctB2R\\em_PxOctB2R_model.pdb # this will become your #2 structure
        open PrOctB2R\\em_PrOctB2R_model.pdb # this will become your #3 structure
        open BmOctB2R\\em_BmOctB2R_model.pdb # this will become your #4 structure
        open CsOctB2R\\em_CsOctB2R_model.pdb # this will become your #5 structure
        matchmaker #2-5 to #1 # superimpose all the #2, #3, #4, and #5 structures with the CPB octopamine beta 2 receptor (reference structure #1)
        select #1/A:105,109,110,113,188,194,198,199,202,290,293,294,297,313,317,321; show sel target ab # select the ligand binding residues of CPB OctB2R and show them in atom and bond (stick) representation
        select #2/A:100,104,105,108,183,189,193,194,197,282,285,286,289,306,310,314; show sel target ab # select the ligand binding residues of PxOctB2R and show them in atom and bond (stick) representation
        select #3/A:104,108,109,112,185,191,195,196,199,286,289,290,293,309,313,317; show sel target ab # select the ligand binding residues of PrOctB2R and show them in atom and bond (stick) representation
        select #4/A:111,115,116,119,192,198,202,203,206,293,296,297,300,316,320,324; show sel target ab # select the ligand binding residues of BmOctB2R and show them in atom and bond (stick) representation
        select #5/A:113,117,118,121,194,200,204,205,208,295,298,299,302,318,322,326; show sel target ab # select the ligand binding residues of CsOctB2R and show them in atom and bond (stick) representation
        select #1/A; cartoon hide sel # hide the cartoon representation of structure #1
        select #2/A; cartoon hide sel # hide the cartoon representation of structure #2
        select #3/A; cartoon hide sel # hide the cartoon representation of structure #3
        select #4/A; cartoon hide sel # hide the cartoon representation of structure #4
        select #5/A; cartoon hide sel # hide the cartoon representation of structure #5
        delete H # delete all the hydrogen atoms attached the residues to simplify the ligand binding site representation otherwise the superposed ligand binding site is going to look very crowded and complex
        select #1/A:105 | #2/A:100 | #3/A:104 | #4/A:111 | #5/A:113; 2dlabel create gpcr328 text "3.28x28" # collectively label the residues across 5 insect octopamine beta 2 receptors as 3.28x28
        2dlabel change gpcr328 size 12 # change the size of the label to 12. After changing the size of the label, you can change the 'right mouse' setting in UCSF ChimeraX GUI Interface to use the right click function to move and arrange the label around
        select #1/A:109 | #2/A:104 | #3/A:108 | #4/A:115 | #5/A:117; 2dlabel create gpcr332 text "3.32x32"
        2dlabel change gpcr332 size 12
        select #1/A:110 | #2/A:105 | #3/A:109 | #4/A:116 | #5/A:118; 2dlabel create gpcr333 text "3.33x33"
        2dlabel change gpcr333 size 12
        select #1/A:113 | #2/A:108 | #3/A:112 | #4/A:119 | #5/A:121; 2dlabel create gpcr336 text "3.36x36"
        2dlabel change gpcr336 size 12
        select #1/A:188 | #2/A:183 | #3/A:185 | #4/A:192 | #5/A:194; 2dlabel create gpcr4552 text "45.52x52"
        2dlabel change gpcr4552 size 12
        select #1/A:194 | #2/A:189 | #3/A:191 | #4/A:198 | #5/A:200; 2dlabel create gpcr539 text "5.38x39"
        2dlabel change gpcr539 size 12
        select #1/A:198 | #2/A:193 | #3/A:195 | #4/A:202 | #5/A:204; 2dlabel create gpcr543 text "5.42x43"
        2dlabel change gpcr543 size 12
        select #1/A:199 | #2/A:194 | #3/A:196 | #4/A:203 | #5/A:205; 2dlabel create gpcr544 text "5.43x44"
        2dlabel change gpcr544 size 12
        select #1/A:202 | #2/A:197 | #3/A:199 | #4/A:206 | #5/A:208; 2dlabel create gpcr546 text "5.46x461"
        2dlabel change gpcr546 size 12
        select #1/A:290 | #2/A:282 | #3/A:286 | #4/A:293 | #5/A:295; 2dlabel create gpcr648 text "6.48x48"
        2dlabel change gpcr648 size 12
        select #1/A:293 | #2/A:285 | #3/A:289 | #4/A:296 | #5/A:298; 2dlabel create gpcr651 text "6.51x51"
        2dlabel change gpcr651 size 12
        select #1/A:294 | #2/A:286 | #3/A:290 | #4/A:297 | #5/A:299; 2dlabel create gpcr652 text "6.52x52"
        2dlabel change gpcr652 size 12
        select #1/A:297 | #2/A:289 | #3/A:293 | #4/A:300 | #5/A:302; 2dlabel create gpcr655 text "6.55x55"
        2dlabel change gpcr655 size 12
        select #1/A:313 | #2/A:306 | #3/A:309 | #4/A:316 | #5/A:318; 2dlabel create gpcr734 text "7.35x34"
        2dlabel change gpcr734 size 12
        select #1/A:317 | #2/A:310 | #3/A:313 | #4/A:320 | #5/A:322; 2dlabel create gpcr738 text "7.39x38"
        2dlabel change gpcr738 size 12
        select #1/A:321 | #2/A:314 | #3/A:317 | #4/A:324 | #5/A:326; 2dlabel create gpcr742 text "7.43x42"
        2dlabel change gpcr742 size 12
        name site1 #1/A:105,109,110,113,188,194,198,199,202,290,293,294,297,313,317,321 # name the ligand binding site of CPB OctB2R as site1
        name site2 #2/A:100,104,105,108,183,189,193,194,197,282,285,286,289,306,310,314 # name the ligand binding site of PxOctB2R as site2
        name site3 #3/A:104,108,109,112,185,191,195,196,199,286,289,290,293,309,313,317 # name the ligand binding site of PrOctB2R as site3
        name site4 #4/A:111,115,116,119,192,198,202,203,206,293,296,297,300,316,320,324 # name the ligand binding site of BmOctB2R as site4
        name site5 #5/A:113,117,118,121,194,200,204,205,208,295,298,299,302,318,322,326 # name the ligand binding site of CsOctB2R as site5
        matchmaker site2 to site1 # perform pairwise structural superposition between the ligand binding residues of CPB OctB2R and other insect species OctB2R (superpose only on the ligand binding site itself, not on the whole structure itself)
        matchmaker site3 to site1
        matchmaker site4 to site1
        matchmaker site5 to site1
        preset "overall look" "publication 1 (silhouettes)"
        select clear # remember to clear selection before rendering and saving the image
        save "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/ligand_binding_site_identification_results/superimposed ligand binding site of CPB octopamine beta 2 receptor with other insect OctB2R.cxs"
        """, language="bash")
        st.markdown("[Visit the GPCRdb generic residue numbering (PDB) section to assign the GPCRdb generic numbers, and Ballesteros-Weinstein numbering to your input PDB GPCR structure. Remember to download the pymol script (pymol_view_generic_numbers.pml) to view the assigned generic numbers via PyMol to double check your identified ligand binding residues](https://gpcrdb.org/structure/generic_numbering_index)")
        st.markdown("[Visit the GitHub page of gpcrmining python package](https://github.com/drorlab/GPCR-mining)")
        st.markdown("[Visit the GPCRdb website to view the list of ligand binding residues of 2RH1](https://gpcrdb.org/interaction/2RH1)")
        st.markdown("[Visit the GPCRdb website to view the diagram of 7TM of ADRB2_HUMAN receptor by hovering your mouse cursor over each residue on the diagram to display the residue's GPCR BW numbers](https://gpcrdb.org/protein/adrb2_human/)")

        st.write("###")

        st.write("**17. Perform pairwise sequence alignment between the CPB octopamine beta 2 receptor and 2RH1 via MUSCLE to identify whether the conserved motifs and seven transmembrane helices are conserved and aligned well to each other (to convince reviewer why you use 2RH1 (a type of ADRB2_HUMAN receptor) as the template when it comes to predicting the ligand binding site of CcOctB2R). Perform structural superimposition between the CPB octopamine beta 2 receptor and 2RH1 via UCSF ChimeraX to compute the refined RMSD value**")
        st.write("✔️concatenate the CPB octopamine beta 2 receptor FASTA sequence with the 2RH1 FASTA sequence")
        st.code("""
        cat CPB_OctB2R_TRINITY_DN22425_c0_g1_i3_p1.fasta rcsb_pdb_2RH1.fasta > combined_CCOctB2R_and_2RH1_sequence.fasta # navigate to the working directory, /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template before you run this command
        """, language="bash")
        st.write("✔️rename the FASTA sequence headers")
        st.code("""
        sed 's/^>TRINITY_DN22425_c0_g1_i3.p1.*/>CcOct\u03B22R/; s/^>2RH1.*/>2RH1/' /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/combined_CCOctB2R_and_2RH1_sequence.fasta > /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/renamed_combined_CCOctB2R_and_2RH1_sequence.fasta
        """, language="bash")
        st.write("✔️run MSA via MUSCLE and view the MSA results (.aln) via JalView")
        st.code("""
        nohup muscle -align /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/renamed_combined_CCOctB2R_and_2RH1_sequence.fasta -output /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template.aln -threads 16 -consiters 2 -refineiters 100 > /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template_output.log 2>&1 &
        """, language="bash")
        st.write("✔️run TMHMM to predict the transmembrane domains of the all the octopamine beta 2 receptor sequences to check whether they have 7 transmembrane domains of GPCR. Retain only those with 7 transmembrane helices (TMHs) for downstream analysis**")
        st.code("""
        cd /media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin # navigate into this working directory before you run tmhmm
        export PATH="/media/raid/Wee/WeeYeZhi/output/trinotate_results/tmhmm_installation/tmhmm-2.0c/bin:$PATH" # note that the export of the PATH environment variable to the tmhmm executable file is just temporary and will only take effect for the current terminal session
        which tmhmm # double check the system can access the tmhmm executable file correctly after exporting the PATH variable 
        tmhmm /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/renamed_combined_CCOctB2R_and_2RH1_sequence.fasta > /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/MSA_between_CPB_octopamine_beta2_receptor_and_2RH1_template/tmhmm_results/tmhmm_CPB_octopamine_beta2_receptors_and_2RH1_results.out
        """, language="bash")
        st.code("""
        pwd # print the current working directory of UCSF ChimeraX to check your file location now
        cd C:\\Users\\USER\\Documents\\Master_of_Science_in_Systems_Biology_GRA\\Experimental_results\\ligand_binding_site_identification_results\\superimposition_between_energy_minimized_AlphaFold3_predicted_CPB_octopamine_beta2_receptor_and_2RH1 # navigate into the desired working directory containing the predicted 3D protein structures
        open 2RH1.pdb (This will become your #1 structure)
        open energy_minimized_.pdb (This will become your #2 structure)
        select #1; hide sel target a # select your first loaded protein structure (#1) and hide all of its atom and stick representation
        matchmaker #2 to #1 # perform structural superimposition to structurally superimpose the AlphaFold3-predicted CPB octopamine beta 2 receptor (labelled as #2) to 2RH1 as the reference structure (labelled as #1) to identify the ligand binding residues of AlphaFold3-predicted CPB octopamine beta 2 receptor that superimpose with the residues of 2RH1 (W116,D120,D121,V124,F200,Y206,S210,S211,S214,W421,F424,F425,N428,Y443,N447,Y451). In this case, #1 is the reference structure.
        select #1/A:116,120,121,124,200,206,210,211,214,421,424,425,428,443,447,451; show sel target ab # select the ligand binding residues of 2RH1 and show them in atom and bond (stick) representation
        manually hover your mouse over each of the ligand binding residue of 2RH1 and record down the superimposed ligand binding residue of the AlphaFold3-predicted CPB octopamine beta 2 receptor
        select #2/A:105,109,110,113,188,194,198,199,202,290,293,294,297,313,317,321; show sel target ab # select the ligand binding residues of CPB OctB2R and show them in atom and bond (stick) representation
        select #1/A; cartoon hide sel # hide the cartoon representation of structure #1
        select #2/A; cartoon hide sel # hide the cartoon representation of structure #2
        select #1/A:116,120,121,124,200,206,210,211,214,421,424,425,428,443,447,451; label sel text "{0.name} {0.number}{0.insertion_code}"
        select #2/A:105,109,110,113,188,194,198,199,202,290,293,294,297,313,317,321; label sel text "{0.name} {0.number}{0.insertion_code}"
        ui mousemode right "move label" # change the mouse setting so that you can right click the amino acid residue label and move them to arrange properly
        save "C:/Users/USER/Documents/Master_of_Science_in_Systems_Biology_GRA/Experimental_results/ligand_binding_site_identification_results/superimposition between ligand binding residues of 2RH1 and CPB OctB2R.cxs"
        """, language="bash")

        st.write("###")

        st.write("**18. install GalaxyRefine and run it locally to refine the structure of the AlphaFold3-predicted, SwissModel-predicted, and GPCR-I-TASSER-predicted CPB octopamine beta 2 receptor respectively. After finished running GalaxyRefine, view the results in UCSF ChimeraX and extract the first structure (delete the other 4 structures) for downstream protein structure evaluation.**")
        st.code("""
        navigate to /media/raid/Wee/WeeYeZhi/output/galaxyrefine_installation
        fill up the galaxyweb software request form at https://galaxy.seoklab.org/request_softwares.html to request for the GalaxyRefine software installation files & download the software installation files inside the installation directory, /media/raid/Wee/WeeYeZhi/output/galaxyrefine_installation
        unzip GalaxyRefine.zip
        cd GalaxyRefine-master/
        export GALAXY_HOME=/media/raid/Wee/WeeYeZhi/output/galaxyrefine_installation/GalaxyRefine-master # set the environment variable for GALAXY_HOME 
        export PATH=$GALAXY_HOME/bin:$PATH
        chmod +x $GALAXY_HOME/bin/*
        $GALAXY_HOME/bin/GalaxyRefine -h # output the usage information of GalaxyRefine
        mkdir /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output/CPB_octopamine_beta2_receptor/CPB_OctB2R/five_top_ranked_AlphaFold3_model/CPB_OctB2R_seed-1_sample-2_model.pdb -t GalaxyRefine_CPBOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/galaxyrefine_CPBOctB2R_results.log 2>&1 & # refine the AlphaFold3-predicted CPB octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/swissmodel_prediction_results/swissmodel_predicted_CPB_octopamine_beta2_receptor.pdb -t GalaxyRefine_SwissModel_Predicted_CPBOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/galaxyrefine_swissmodel_predicted_CPBOctB2R_results.log 2>&1 & # refine the SwissModel-predicted CPB octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/gpcritasser_prediction_results/gpcritasser_predicted_CPB_octopamine_beta2_receptor.pdb -t GalaxyRefine_GPCRITASSER_Predicted_CPBOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/galaxyrefine_GPCRITASSER_predicted_CPBOctB2R_results.log 2>&1 & # refine the GPCRItasser-predicted CPB octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/alphafold2_prediction_results/AlphaFold2_predicted_CcOctB2R_2c236_relaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb -t GalaxyRefine_AlphaFold2_Predicted_CPBOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/galaxyrefine_AlphaFold2_predicted_CPBOctB2R_results.log 2>&1 & # refine the AlphaFold2-predicted CPB octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/AF-G3M4F8-F1-model_v6.pdb -t GalaxyRefine_CsOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/GalaxyRefine_CsOctB2R/galaxyrefine_CsOctB2R_results.log 2>&1 & # refine the Chilo suppressalis octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/Bm_OctB2R_seed-1_sample-2_model.pdb -t GalaxyRefine_BmOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/GalaxyRefine_BmOctB2R/galaxyrefine_BmOctB2R_results.log 2>&1 & # refine the Bombyx mori octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/Pr_OctB2R_seed-1_sample-3_model.pdb -t GalaxyRefine_PrOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/GalaxyRefine_PrOctB2R/galaxyrefine_PrOctB2R_results.log 2>&1 & # refine the Pieris rapae octopamine beta 2 receptor
        nohup $GALAXY_HOME/bin/GalaxyRefine -p /media/raid/Wee/WeeYeZhi/output/ligand_binding_site_identification_results/Px_OctB2R_seed-3_sample-3_model.pdb -t GalaxyRefine_PxOctB2R -s 16 -o 5 > /media/raid/Wee/WeeYeZhi/output/galaxyrefine_results/GalaxyRefine_PxOctB2R/galaxyrefine_PxOctB2R_results.log 2>&1 & # refine the Plutella xylostella octopamine beta 2 receptor
        """, language="bash")
        # ----LOAD THE GALAXYREFINE USAGE INFORMATION FILE----
        # Check if the file exists before reading
        if galaxyrefine_usage_information_file.exists():
            with open(galaxyrefine_usage_information_file, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download The GalaxyRefine Usage Information File",
                data=script_byte,
                file_name=galaxyrefine_usage_information_file.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{galaxyrefine_usage_information_file.name} does not exist.")

        st.write("###")

        st.write("**19. run Yasara energy minimization of the Galaxy-refined AlphaFold3-predicted, SwissModel-predicted, and GPCR-I-TASSER-predicted CPB octopamine beta 2 receptor respectively via YASARA webserver and download Yasara View to view and save the YASARA-energy minimized CPB octopamine beta 2 receptor model in .PDB format**")
        st.markdown("[Visit YASARA energy minimization webserver](https://www.yasara.org/minimizationserver.htm)")
        st.markdown("[Visit YASARA tutorial](https://ceberndsen.github.io/YASARA-guide/file-types-and-how-to-work-with-them-in-yasara.html)")
        st.markdown("[Visit YASARA View webpage to download YASARA View](https://www.yasara.org/viewdl.htm)")
        st.markdown("[Read to understand the PDB file format](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)")

        st.write("###")

        st.write("**20. perform 3D protein structure evaluation to evaluate the 3D structure of the CPB octopamine beta 2 receptor predicted by AlphaFold3, Swiss Model, and GPCR-I-TASSER respectively via SAVESserver, ProSA2003, MolProbity, and QMEANBrane**")
        st.write("✔️install ProSa2003 (since the ProSA-Web server is closed temporarily due to attack by hackers)")
        st.code("""
        install ProSa2003 via Window # choose the ProSa2003 (Win XP/2000)
        check your email and copy the academic license key into the LicenseDB file and save the LicenseDB file (without any file extension) in the .\ProSaData\LicenseDB\. working directory
        unzip prosa2003.zip # download the prosa2003.zip file and unzip it
        double click the Setup.exe to download the Window version of ProSA2003
        """, language="bash")
        st.write("✔️compute the zscore of the energy minimized AlphaFold3-predicted CPB octopamine beta 2 receptor via ProSa2003")
        st.code("""
        init zscore # after you open the ProSA2003 application, run this command to begin computing the Z-score of the predicted CPB octopamine beta 2 receptor
        read pdb energy_minimized_alphafold3_predicted_CPB_octopamine_beta2_receptor.pdb AlphaFold3_predicted_CPB_octopamine_beta2_receptor # remember to move the file, energy_minimized_alphafold3_predicted_CPB_octopamine_beta2_receptor.pdb into the ProSA installation PDB directory, /mnt/c/prosa_installation/ProSa2003/PDB. In this case, name the receptor as AlphaFold3_predicted_CPB_octopamine_beta2_receptor
        zscore AlphaFold3_predicted_CPB_octopamine_beta2_receptor AlphaFold3_predicted_CPB_octopamine_beta2_receptor_zscore # compute the Z-score of the CPB octopamine beta 2 receptor and save it. All the files will be saved in the working directory, /mnt/c/prosa_installation/ProSa2003
        analyse energy AlphaFold3_predicted_CPB_octopamine_beta2_receptor # calculate the energy of AlphaFold3_predicted_CPB_octopamine_beta2_receptor
        plot # display the energy graph of AlphaFold3_predicted_CPB_octopamine_beta2_receptor
        export plot AlphaFold3_predicted_CPB_octopamine_beta2_receptor_plot # export a postscript file called obj1_plot.ps and save it. All the files will be saved in the working directory, /mnt/c/prosa_installation/ProSa2003
        sudo apt install ghostscript # install ghostscript to handle the postscript file
        ps2pdf AlphaFold3_predicted_CPB_octopamine_beta2_receptor_plot.ps AlphaFold3_predicted_CPB_octopamine_beta2_receptor_plot.pdf # convert the postscript file to pdf version to view the plot
        """, language="bash")
        st.write("✔️compute the zscore of the energy minimized SwissModel-predicted CPB octopamine beta 2 receptor via ProSa2003")
        st.write("✔️compute the zscore of the energy minimized GPCR-I-TASSER-predicted CPB octopamine beta 2 receptor via ProSa2003")
        st.write("✔️use the Microsoft Excel advanced filter function to filter the DALI web server results table by the Z-score (from largest to smallest) first, followed by percentage of sequence identity (%id) (from largest to smallest), and rmsd columns (from smallest to largest)")
        st.markdown("[Visit DALI web server](http://ekhidna2.biocenter.helsinki.fi/dali/)")
        st.markdown("[Visit SAVESserver v6.1 webpage](https://saves.mbi.ucla.edu/)")
        st.markdown("[Visit MolProbity](http://molprobity.biochem.duke.edu/)")
        st.markdown("[Visit Swiss Model 3D Protein Structure Assessment](https://swissmodel.expasy.org/assess)")
        st.markdown("[Visit Swiss Model QMEAN (Qualitative Model Energy Analysis) to run QMEANBrane](https://swissmodel.expasy.org/qmean/)")
        st.markdown("[Visit ProSA2003 download page](https://www.came.sbg.ac.at/prosa.php)")
        # ----LOAD ProSa2003 manual----
        # Check if the file exists before reading
        if prosa2003_manual.exists():
            with open(prosa2003_manual, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download ProSa2003 manual to know how to run the standalone version of ProSa2003 using command lines (since ProSa web server has been temporarily shut down)",
                data=script_byte,
                file_name=prosa2003_manual.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{prosa2003_manual.name} does not exist.")
        # ----LOAD QMEANBrane interpretation manual----
        # Check if the file exists before reading
        if qmeanbrane_manual.exists():
            with open(qmeanbrane_manual, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download QMEANBrane manual to interpret QMEANBrane results",
                data=script_byte,
                file_name=qmeanbrane_manual.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{qmeanbrane_manual.name} does not exist.")

        st.write("###")

        st.write("**21. Run GROMACS**")
        st.write("✔️check whether the NVIDIA GPU driver is set up and properly configured")
        st.code("""
        gmx mdrun -version # check whether GROMACS build was compiled with GPU support by checking the line, GPU support:
        nvidia-smi # check the NVIDIA GPU driver version and CUDA version. Check whether the NVIDIA Driver is installed, check whether GPU is detected, and check whether the CUDA runtime exists. The driver version of the NVIDIA GPU is 575.57.08 and CUDA version is 12.9. 
        nvcc --version # check the version of the CUDA compiler to check whether the CUDA toolkit is fully installed. The version of the NVIDIA CUDA compiler driver is 12.9.r12.9. CUDA is the backend of NVIDIA GPU driver that allows GROMACS to communicate with and utilize the NVIDIA GPU driver
        modinfo nvidia | grep -w version # check the version of the NVIDIA GPU driver. The driver version is 575.57.08
        """, language="bash")
        st.write("✔️build GPU-enabled GROMACS from source code and configure GROMACS to utilize GPU while running")
        st.code("""
        wget ftp://ftp.gromacs.org/gromacs/gromacs-2026.2.tar.gz # navigate to the working directory, /media/raid/Wee/WeeYeZhi/output/GPU_enabled_GROMACS_installation, followed by downloading the source code of GROMACS
        sudo apt-get install cmake build-essential snap gcc-11 g++-11 libopenmpi-dev # download the dependencies required for building GROMACS to run GROMACS using GPU
        tar -xvzf gromacs-2026.2.tar.gz # extract the downloaded GROMACS folder
        cd gromacs-2026.2
        mkdir build
        cd build
        cmake .. -DGMX_GPU=CUDA -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11 -DGMX_BUILD_OWN_FFTW=ON -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.9 # tell GROMACS that it is a CUDA configured GPU & build GROMACS on GPU. Remember to mention the GPU configuration using the -DGMX_GPU tag, mention the C and C++ compilers used by the GROMACS to build the C files using the -DCMAKE-C_COMPILER & -DCMAKE_CXX_COMPILER tag, mention the build protocol that should be used by GROMACS to build the GROMACS on GPU using the -DGMX_BUILD_OWN=ON tag, mention the CUDA toolkit path using the -DCUDA_TOOLKIT_ROOT_DIR tag. Run the command, which nvcc, to check the path of your system's CUDA toolkit. The -DGMX_GPU=CUDA tag builds GROMACS with NVIDIA GPU support enabled
        grep GMX_GPU CMakeCache.txt # look for the line, GMX_GPU:STRING-CUDA to verify that CMake has successfully detected CUDA and configured GROMACS for NVIDIA GPU support
        nohup make -j$(nproc) > build_GPU_enabled_GROMACS.log 2>&1 & # compile all the C files of GROMACS by using all the CPU cores simultaneously during compilation
        sudo make install # install the newly compiled GROMACS 
        source /usr/local/gromacs/bin/GMXRC # activate the newly compiled gromacs. You need to activate the newly compiled GROMACS everytime you open the new terminal session
        which gmx # check which GROMACS installation is active
        gmx mdrun -version # verify whether the GPU support has been successfully enabled. check whether GROMACS successfully detected and compiled with NVIDIA GPU CUDA support
        nvidia-smi # During MD simulation, monitor GPU utilization. Check the utilization percentage of GPU using another Ubuntu terminal session while you are running the GROMACS mdrun program at the background (to run energy minimization, MD simulation, and calculation) to see whether GROMACS is run using GPU. If you see there's a certain utilization percentage of GPU, then it means that you have successfully run GROMACS using GPU instead of using the default CPU
        """, language="bash")
        st.markdown("[Visit Official GROMACS Download Page](https://manual.gromacs.org/current/download.html)")
        st.markdown("[Visit Official GROMACS Installation Guide Page to build a GPU-enabled GROMACS](https://manual.gromacs.org/current/install-guide/index.html)")
        st.markdown("[Watch YouTube video to learn how to build a GPU-enabled GROMACS. Note that the -DGMX_BUILD_OWN_FFTW=ON is the correct tag instead of specifying the outdated tag, -DGMX_BUILD_OWN=ON](https://www.youtube.com/watch?v=If1D2BIYgFg&t=1s)")
        st.write("Note: Read the two publications, 'A beginner's guide to molecular dynamics simulations and the identification of cross-correlation networks for enzyme engineering' and 'Introductory tutorials for simulating protein dynamics with GROMACS' to understand how to run MD simulation")

        st.write("###")

        st.write("**22. run MD simulation of the CPB OctB2R-lipid bilayer complex via GROMACS for 100 ns first, followed by extending the MD simulation to 500 ns and to 1000 ns**")
        st.write("✔️invoke the grompp and mdrun programme to run energy minimization of the receptor-lipid bilayer complex system to lower the overall potential energy of the system by removing overlapping atoms, steric clashes, suboptimal hydrogen bonding, and bad contacts")
        st.code("""
        nohup gmx grompp -f step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p topol.top -o em.tpr > em_tpr.log 2>&1 & # navigate to the working directory, /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs, before you run the command
        nohup gmx mdrun -deffnm em -nb gpu > em_mdrun.log 2>&1 & # run the calculation for non-bonded atomic interactions, PME electrostatic interactions, bonded interactions, and coordinate update/integration on GPU to accelerate the overall MD simulation process
        gmx energy -f em.edr -o potential.xvg # execute the energy command to extract the potential energy (as a function of energy minimization step) of the system from the generated binary energy file (.edr). Select 13 (Potential) and 0 (to terminate input). Based on the potential energy plot, you will notice that the potential energy of the system drops dramatically in the first few steps as water molecules rearrange to optimize hydrogen bonding
        """, language="bash")
        st.write("Note: The receptor-lipid bilayer complex system was energy-minimized using the steepest descent method under positional and dihedral restraints. Positional restraint force constants of 4000, 2000, and 1000 kJ mol⁻¹ nm⁻² were applied to protein backbone atoms, protein side-chain atoms, and lipid molecules, respectively, while dihedral restraints were applied with a force constant of 1000 kJ mol⁻¹. This is to prevent both the protein and lipid bilayer membrane from distorting too much while bad steric clashes are removed. The maximum energy minimization step was set as 5000 and the energy minimization process stops when the maximum force reaches 1000 kJ mol⁻¹ nm−1. The neighbour list, van der Waals, and electrostatic cutoff were all set as 1.2 nm. Under periodic boundary conditions, the long-range electrostatic interactions were calculated using the Particle Mesh Ewald (PME) method with a real-space cutoff of 1.2nm. The LINCS method was used to constrain hydrogen-containing bonds since hydrogen atoms tend to vibrate very rapidly.")

        st.write("###")

        st.write("✔️invoke the grompp and mdrun programme to run NVT and NPT equilibration of the receptor-lipid bilayer complex system to bring the whole system to a stable temperature and pressure")
        st.code("""
        nohup gmx grompp -f step6.1_equilibration.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o step6.1.tpr > step6.1_tpr.log 2>&1 & # run first NVT equilibration
        nohup gmx mdrun -s step6.1.tpr -deffnm step6.1 -nb gpu -pme gpu -bonded gpu -update gpu > step6.1_mdrun.log 2>&1 &
        nohup gmx grompp -f step6.2_equilibration.mdp -c step6.1.gro -r step6.1.gro -t step6.1.cpt -p topol.top -n index.ndx -o step6.2.tpr > step6.2_tpr.log 2>&1 & # run second NVT equilibration
        nohup gmx mdrun -s step6.2.tpr -deffnm step6.2 -nb gpu -pme gpu -bonded gpu -update gpu > step6.2_mdrun.log 2>&1 &
        gmx energy -f step6.2.edr -o temperature_step6.2.xvg # extract the temperature of the system during NVT equilibration to plot a graph of system temperature against time. Select 17 (temperature) and 0 (to terminate input) 
        nohup gmx grompp -f step6.3_equilibration.mdp -c step6.2.gro -r step6.2.gro -t step6.2.cpt -p topol.top -n index.ndx -o step6.3.tpr > step6.3_tpr.log 2>&1 & # run first NPT equilibration
        nohup gmx mdrun -s step6.3.tpr -deffnm step6.3 -nb gpu -pme gpu -bonded gpu -update gpu > step6.3_mdrun.log 2>&1 &
        nohup gmx grompp -f step6.4_equilibration.mdp -c step6.3.gro -r step6.3.gro -t step6.3.cpt -p topol.top -n index.ndx -o step6.4.tpr > step6.4_tpr.log 2>&1 & # run second NPT equilibration
        nohup gmx mdrun -s step6.4.tpr -deffnm step6.4 -nb gpu -pme gpu -bonded gpu -update gpu > step6.4_mdrun.log 2>&1 &
        nohup gmx grompp -f step6.5_equilibration.mdp -c step6.4.gro -r step6.4.gro -t step6.4.cpt -p topol.top -n index.ndx -o step6.5.tpr > step6.5_tpr.log 2>&1 & # run third NPT equilibration
        nohup gmx mdrun -s step6.5.tpr -deffnm step6.5 -nb gpu -pme gpu -bonded gpu -update gpu > step6.5_mdrun.log 2>&1 &
        nohup gmx grompp -f step6.6_equilibration.mdp -c step6.5.gro -r step6.5.gro -t step6.5.cpt -p topol.top -n index.ndx -o step6.6.tpr > step6.6_tpr.log 2>&1 & # run fourth NPT equilibration
        nohup gmx mdrun -s step6.6.tpr -deffnm step6.6 -nb gpu -pme gpu -bonded gpu -update gpu > step6.6_mdrun.log 2>&1 &
        gmx energy -f step6.6.edr -o pressure_step6.6.xvg # extract the pressure of the system during NPT equilibration to plot a graph of system pressure against time. Select 18 (pressure) and 0 (to terminate input)
        """, language="bash")
        st.write("Note: Following energy minimization, the insect receptor–lipid bilayer complex system was subjected to a six-step equilibration protocol consisting of initial NVT and subsequent NPT equilibration phases using GROMACS. The first two equilibration stages were performed under the canonical ensemble (NVT) for 125 ps each using a 1 fs integration timestep, whereas the subsequent four stages were conducted under the isothermal–isobaric ensemble (NPT) with semiisotropic pressure coupling. The first NPT stage employed a 1 fs timestep for 125 ps, followed by three additional NPT stages performed with a 2 fs timestep for 500 ps each, resulting in a total equilibration time of 1.875 ns. Temperature was maintained at 300 K using the velocity-rescale thermostat with separate coupling groups defined for the solute, membrane, and solvent components. Pressure equilibration was carried out using the stochastic cell-rescale barostat with semiisotropic coupling at 1 bar and a compressibility of 4.5 × 10⁻⁵ bar⁻¹. Long-range electrostatic interactions were calculated using the Particle Mesh Ewald (PME) method under periodic boundary conditions with a real-space cutoff of 1.2 nm, while van der Waals interactions were treated using a force-switch cutoff scheme between 1.0 and 1.2 nm. Hydrogen-containing bonds were constrained using the LINCS algorithm. During equilibration, positional restraints on protein backbone atoms, protein side chains, lipid molecules, and dihedral angles were gradually reduced across the six equilibration stages to allow progressive relaxation and stabilization of the membrane protein system prior to the production molecular dynamics simulation. The six-stage NVT and NPT equilibration protocol was performed to gradually relax and stabilize the insect receptor–lipid bilayer complex system prior to production molecular dynamics simulation. Initial NVT equilibration allowed the system temperature to stabilize while maintaining a constant volume, whereas subsequent NPT equilibration enabled pressure and system density adjustment under semiisotropic conditions suitable for membrane systems. The gradual reduction of positional and dihedral restraints across the six equilibration stages minimized structural distortion, prevented membrane instability, and allowed progressive adaptation of the protein, lipid bilayer, water molecules, and ions toward thermodynamically stable conditions.")
        st.write("Note: The energy command is executed during the NVT equilibration to monitor temperature of the system whereas the energy command is executed during the NPT equilibration to monitor the pressure, density, box dimension, and potential energy of the system throughout the trajectory")

        st.write("###")

        st.write("✔️run the MD production simulation for 100 ns first, followed by extending the MD simulation to 500 ns and to 1000 ns")
        st.code("""
        nohup gmx grompp -f step7_production.mdp -c step6.6.gro -t step6.6.cpt -p topol.top -n index.ndx -o md_100ns.tpr > step7.0_tpr.log 2>&1 & # generate the portable binary run input file:
        nohup gmx mdrun -s md_100ns.tpr -deffnm md_100ns -nb gpu -pme gpu -bonded gpu -update gpu > step7.0_100ns_mdrun.log 2>&1 & # run the 100 ns production simulation with GPU acceleration
        nohup gmx convert-tpr -s md_100ns.tpr -extend 400000 -o md_500ns.tpr > extend_400000_mdrun.log 2>&1 & # extend the simulation from 100 ns → 500 ns
        nohup gmx mdrun -s md_500ns.tpr -deffnm md_500ns -cpi md_100ns.cpt -noappend -nb gpu -pme gpu -bonded gpu -update gpu > step7.0_500ns_mdrun.log 2>&1 & # continue running the simulation up to 500 ns with GPU acceleration (to preserve the atom's velocities, thermostat state, barostat state, random seeds, and trajectory continuity)
        nohup gmx convert-tpr -s md_500ns.tpr -extend 500000 -o md_1000ns.tpr > extend_500000_mdrun.log 2>&1 & # extend the simulation from 500 ns → 1000 ns
        nohup gmx mdrun -s md_1000ns.tpr -deffnm md_1000ns -cpi md_500ns.cpt -noappend -nb gpu -pme gpu -bonded gpu -update gpu > step7.0_1000ns_mdrun.log 2>&1 & # continue running the simulation up to 1000 ns with GPU acceleration
        nohup gmx trjcat -f md_100ns.xtc md_500ns.part0002.xtc md_1000ns.part0003.xtc -o md_full_1000ns.xtc > concatenate_trajectory.log 2>&1 & # concatenate trajectories together
        nohup gmx eneconv -f md_100ns.edr md_500ns.part0002.edr md_1000ns.part0003.edr -o md_full_1000ns.edr > concatenate_energy.log 2>&1 & # concatenate energy files together
        gmx trjconv -s md_100ns.tpr -f md_full_1000ns.xtc -o md_reimage.xtc -pbc mol -ur compact -center # Select group for centering: Protein & Select group for output: System. The md_reimage.xtc file is the corrected whole system trajectory. Under periodic boundary conditions, when molecules cross the box boundaries, the protein may look split, the lipid may look split, and the water may appear discontinuous, but physically this is nothing wrong as it is just visualization artifact. Hence, the goal of reimaging is to reconstruct the whole molecule, center the protein, and keep the membrane visually intact
        gmx trjconv -s md_100ns.tpr -f md_reimage.xtc -o protein_fit.xtc -fit rot+trans # Select group for least squares fit: Backbone & Select group for output: Protein. The protein_fit.xtc becomes the main analysis trajectory file used for analyzing RMSD, RMSF, Hbond, PCA, clustering, and ensemble docking since it has PBC artifacts removed, protein centered, global rotation removed, and global translation removed. Even after reimaging, the protein may still rotate and the whole system may still drift, which interferes with the downstream RMSD analysis. The goal of fitting is to remove global rotation and global translation while preserving the internal conformational motions. Fitting to backbone removes overall drifting, box movement, and rotational artifacts while preserving biologically meaningful motions
        gmx convert-tpr -s md_100ns.tpr -o protein.tpr # Select for group to be written out: Protein to extract protein-only atoms, protein-only topology information, protein-only coordinates, and protein groups. The protein.tpr file can be used for subsequent visualization and analysis
        watch -n 1 nvidia-smi # monitor the GPU usage (monitor GPU utilization, VRAM usage, temperature, and power draw)
        """, language="bash")
        st.write("Note: Following equilibration, a 100 ns production molecular dynamics simulation of the insect receptor–lipid bilayer complex system was performed using GROMACS with a 2 fs integration timestep under periodic boundary conditions. Temperature was maintained at 300 K using the velocity-rescale thermostat, while pressure was maintained at 1 bar using the semiisotropic stochastic cell-rescale barostat suitable for membrane systems. Long-range electrostatic interactions were calculated using the Particle Mesh Ewald (PME) method with a real-space cutoff of 1.2 nm, whereas van der Waals interactions were treated using a force-switch cutoff scheme between 1.0 and 1.2 nm. Hydrogen-containing bonds were constrained using the LINCS algorithm. The production simulation was subsequently extended to 500 ns and 1000 ns using checkpoint continuation to preserve system velocities and simulation states for long-timescale conformational analysis.")
        st.write("Note: The gmx convert-tpr -extend does not restart from zero, erase your previous trajectory, randomize velocities again, or create a new independent simulation, but instead it only modifies the total allowed runtime stored inside the .tpr file")

        st.write("###")

        st.write("✔️create an index group for TM helices and three ECL residues (residue 188, 190 & 194) only (since they are ligand binding residues), followed by extracting the topology information and MD simulation trajectory of TM helices and three ECL residues (residue 188, 190 & 194) only. Expect to observe lower RMSD values since loops no longer dominate the motion, termini no longer inflate RMSD, and only the GPCR core is measured")
        st.code("""
        gmx make_ndx -f protein.tpr -o tm_helices_ecl.ndx # create TM and ECL index. Type, r 31-54 | r 66-87 | r 104-125 | r 146-167 | r 188 | r 190 | r 194 | r 195-215 | r 279-298 | r 311-334, and press "Enter", followed by typing, name 10 TM_ecl_helices, to rename the newly created index group as TM_ecl_helices and position the newly created index group as 10. Type "q" to quit
        gmx convert-tpr -s protein.tpr -n tm_helices_ecl.ndx -o tm_helices_ecl.tpr # Select 10 (TM_ecl_helices) to extract TM_ecl-only topology
        gmx trjconv -s protein.tpr -f protein_fit.xtc -n tm_helices_ecl.ndx -o tm_helices_ecl_fit.xtc # Select 10 (TM_ecl_helices) to extract TM_ecl-only trajectory
        """, language="bash")

        st.write("###")

        st.write("✔️analyze the MD simulation of the CPB OctB2R-lipid bilayer complex system using XMGrace, VMD, and GNU Plot. RMSD: Is the protein stable? RMSF: Which residues fluctuate the most? Hbond: Are stabilizing interactions maintained? DSSP: Does secondary structure remain stable? PCA: What are the dominant large-scale motions?")
        st.write("✔️calculate the TM-helices only's and protein's backbone RMSD to investigate how much did the protein structure deviate from the reference structure over time to determine structural stability, equilibration quality, large conformational changes, and simulation convergence. To answer questions like does the insect receptor structure maintain its structure? does the binding pocket collapse? does the transmembrane region stay stable?")
        st.code("""
        gmx rms -s tm_helices_ecl.tpr -f tm_helices_ecl_fit.xtc -tu ns -o tm_ecl_rmsd.xvg # Select group for least squares fit: Backbone and Select group for RMSD calculation: Backbone
        gmx rms -s protein.tpr -f protein_fit.xtc -tu ns -o whole_protein_rmsd.xvg # Select group for least squares fit: Backbone and Select group for RMSD calculation: Backbone (Choose backbone because backbone represents the overall protein fold and it is less noisy than the side chains). 
        xmgrace tm_ecl_rmsd.xvg # expect to see lower RMSD values since there are no flexible loops, C-terminus, N-terminus, intracellular loops (ICLs), and extracellular loops (ECLs). Visualize the RMSD plot inside a graphical user interface (GUI) in SCREEN mode, not in REMOTE TERMINAL mode.
        xmgrace whole_protein_rmsd.xvg # expect to see large RMSD values since there are presence of flexible loops, C-terminus, N-terminus, intracellular loops (ICLs), and extracellular loops (ECLs)
        """, language="bash")
        st.write("✔️calculate the TM-helices only's and protein's C-alpha RMSF to investigate how much each residue fluctuates during the simulation to identify flexible regions, stable regions, loops, terminal flexibility, and potentially functional motions (Normally, loops have high RMSF, termini has high RMSF, and helices/core has low RMSF)")
        st.code("""
        gmx rmsf -s tm_helices_ecl.tpr -f tm_helices_ecl_fit.xtc -res -o tm_ecl_rmsf.xvg # Select group: C-alpha (using C-alpha gives residue-level flexibility, reduces side chain noise, and easier biological interpretation.)
        gmx rmsf -s protein.tpr -f protein_fit.xtc -res -o whole_protein_rmsf.xvg # Select group: C-alpha 
        xmgrace tm_ecl_rmsf.xvg
        xmgrace whole_protein_rmsf.xvg
        """, language="bash")
        st.write("✔️investigate the compactness of the protein during the simulation")
        st.code("""
        gmx gyrate -s tm_helices_ecl.tpr -f tm_helices_ecl_fit.xtc -tu ns -o tm_ecl_rg.xvg # Select Protein
        gmx gyrate -s protein.tpr -f protein_fit.xtc -tu ns -o whole_protein_rg.xvg # Select Protein
        xmgrace tm_ecl_rg.xvg
        xmgrace whole_protein_rg.xvg 
        """, language="bash")
        st.write("✔️run sasa analysis to investigate solvent accessibility of protein")
        st.code("""
        gmx sasa -s tm_helices_ecl.tpr -f tm_helices_ecl_fit.xtc -tu ns -o tm_ecl_sasa.xvg # Select Protein 
        gmx sasa -s protein.tpr -f protein_fit.xtc -tu ns -o whole_protein_sasa.xvg # Select Protein 
        xmgrace tm_sasa.xvg
        xmgrace tm_ecl_sasa.xvg
        xmgrace sasa.xvg
        """, language="bash")
        st.write("✔️count the total number of hydrogen bonds formed (since hydrogen bonds help stabilize secondary structure, ligand binding, protein binding, and membrane protein conformations. Stable hydrogen bond numbers usually indicates stable protein fold and preserved secondary structure)")
        st.code("""
        gmx hbond -s tm_helices_ecl.tpr -f tm_helices_ecl_fit.xtc -num tm_ecl_hb_bb.xvg -tu ns # Select MainChain+H for both groups
        gmx hbond -s protein.tpr -f protein_fit.xtc -num whole_protein_hb_bb.xvg -tu ns #  Select MainChain+H for both groups. You ran this command only for your master study. Count the total number of backbone hydrogen bonds (intraprotein backbone Hbonds and secondary structure stability). 
        gmx hbond -s protein.tpr -f protein_fit.xtc -num whole_protein_hb_all.xvg -tu ns # Select Protein for both groups. Count the total number of protein hydrogen bonds. 
        xmgrace tm_ecl_hb_bb.xvg
        xmgrace whole_protein_hb_bb.xvg
        xmgrace whole_protein_hb_all.xvg
        """, language="bash")
        st.write("✔️investigate the time evolution of secondary structure stability throughout the MD simulation. Calculate the average fluctuation statistics of secondary structure content. This is to investigate whether the helices remain stable, any unfolding process occurs, and any secondary structure transition happens, which is very important for membrane proteins")
        st.code("""
        gmx dssp -h # output the usage manual of the GROMACS-installed DSSP
        gmx dssp -s tm_helices_ecl.tpr -f tm_helices_ecl_fit.xtc -tu ns -hmode dssp -o tm_ecl_dssp.dat -num tm_ecl_dssp_num.xvg # Select Protein. The ss.xpm file refers to the secondary structure map across time and residues. The scount.xvg file refers to the count of helix, sheet, bend, turn, loop etc over time
        gmx analyze -f tm_ecl_dssp_num.xvg # calculate the average and error estimates of the assigned secondary structure elements throughout the trajectory
        gmx dssp -s protein.tpr -f protein_fit.xtc -tu ns -hmode dssp -o whole_protein_dssp.dat -num whole_protein_dssp_num.xvg # Select Protein. The ss.xpm file refers to the secondary structure map across time and residues. The scount.xvg file refers to the count of helix, sheet, bend, turn, loop etc over time
        gmx analyze -f whole_protein_dssp_num.xvg # calculate the average and error estimates of the assigned secondary structure elements throughout the trajectory
        """, language="bash")
        st.write("✔️compute the covariance matrix and eigenvectors with eigenvalues to investigate how atoms move during the simulation, either they move upwards together (positive covariance) or one moves upward and one moves downward (negative covariance) during the simulation to perform PCA analysis to determine the dominant collective atomic motions. PCA does not prove biological mechanism automatically or prove protein activation directly, instead it only reveals the dominant collective motions")
        st.code("""
        gmx covar -s protein.tpr -f protein_fit.xtc -o eigenval.xvg -v eigenvec.trr -av average.pdb # Select a group for covariance analysis: Backbone. In this instance, average.pdb is the average structure during simulation, eigenval.xvg is the importance of each motion, and eigenvec.trr is the direction of motion
        """, language="bash")
        st.write("✔️visualize the dominant motions by projecting the trajectories motion onto PC1 and PC2, followed by extracting the extreme conformations (the largest conformational changes/dominant motion) to investigate what does the protein actually do along these motions and how much does the protein move along a particular motion direction. Large projection indicates strong conformational change while small projection indicates little movement. PC1 (eigenvector 1) is the largest motion while PC2 (eigenvector 2) is the second largest motion")
        st.code("""
        gmx anaeig -v eigenvec.trr -s protein.tpr -f protein_fit.xtc -first 1 -last 2 -extr extreme.pdb -proj proj.xvg # Select backbone. In this instance, the proj.xvg is the projection of trajectory along PC1/PC2 whereas the extreme.pdb file is the extreme conformation/motion. The proj.xvg file refers to the trajectory file projected onto PCs. The extreme1.pdb file refers to the extreme motion along PC1 (structure at maximum/minimum along PC1) whereas the extreme2.pdb file refers to the extreme motion along PC2 (structure at maximum/minimum along PC2). The two extreme.pdb files show how the protein moves in its most dominant functional direction. You can answer questions like does PC1 affect binding site? does motion open/close pocket? does ligand stabilize or restrict PC motion (you can compare between apo system PCA vs ligand bound PCA)?
        vmd protein.tpr extreme1.pdb # visualize the protein structure via VMD
        vmd protein.tpr extreme2.pdb
        """, language="bash")
        st.write("Note: To remove the top axis and right axis of the plot rendered by XMGRACE, you can do the following steps. Go to Plot -> Axis Properties -> Click the 'Tick marks' section -> Edit: X-axis -> Draw on: Normal side (located at the 'Placement' section) -> Click 'Apply'. Do the same thing for Edit: Y-axis. After that, Go to Plot -> Graph appearance -> Frame -> Change the Frame Type to Half Open -> Accept. If you notice the four black small square boxes at the four corners of the plot rendered by XMGRACE, it actually does not matter, because after you print out the plot, the boxes will automatically disappear (you do not need to manually or intentionally remove them)")
        st.markdown("[Visit XMGRACE User's Guide v0.5](https://www.uoxray.uoregon.edu/local/manuals/grace/grace-5.0.5/doc/UsersGuide.html)")

        st.write("###")

        st.write("✔️extract the trajectory of the protein-lipid bilayer complex model and protein-only model (without the lipid bilayer) during the period where the RMSD values of TM helix bundle reach plateau (stabilize) (since loops fluctuate naturally, loops inflate RMSD artificially, TM core reflects actual receptor stability, but after determining stability, you still keep the whole receptor, including ECLs and ICLs for biologically realistic docking).")
        st.code("""
        gmx trjconv -s md_100ns.tpr -f md_reimage.xtc -b 900000 -e 1000000 -o whole_protein_lipid_bilayer_900_1000ns.xtc # Select System (the whole protein lipid bilayer complex model). The last 100ns trajectory is considered stable. In this case, the value for -b and -e are in picosecond (ps).
        gmx trjconv -s protein.tpr -f protein_fit.xtc -b 900000 -e 1000000 -o whole_protein_only_900_1000ns.xtc # Select Protein. The last 100ns trajectory is considered stable. In this case, the value for -b and -e are in picosecond (ps).
        """, language="bash")
        st.write("Note: Specifying the tag optionally, -skip 5, to skip frames reduces redundancy, speeds clustering, and usually preserves major conformational states since clustering a full 1000ns trajectory can be computationally heavy")
        st.write("Note: The early unstable frames may create clusters representing relaxation artifacts, non-equilibrated states, unrealistic pocket geometries since you dont want unstable starting conformations, equilibration artifacts, membrane adaptation structures,, and distorted early states. Clustering of stable region gives equilibrated receptor conformations, more reliable binding pockets, more physically meaningful states, and less noise in clustering. Instead of clustering the entire unstable trajectory blindly, this is often preferable for ensemble docking.")

        st.write("###")

        st.write("✔️perform RMSD-based GROMOS clustering analysis using a cutoff of 0.8 Å on Cα atoms to cluster structurally similar protein conformations with similar RMSD values together to identify the dominant receptor conformations throughout the MD trajectory. The three most populated clusters were selected as representative structures for subsequent ensemble docking analyses.")
        st.code("""
        gmx cluster -s protein.tpr -f whole_protein_only_900_1000ns.xtc -n tm_helices_ecl.ndx -method gromos -cutoff 0.12 -o whole_protein_clusters_0.12_900_1000ns.xpm -g whole_protein_cluster_0.12_900_1000ns.log -cl whole_protein_cluster_centers_0.12_900_1000ns.pdb -dist whole_protein_rmsd_dist_0.12_900_1000ns.xvg # Select group for least squares fit: TM_ecl_helices and select group for RMSD calculation: TM_ecl_helices. Select Protein for output. The cutoff value is in the unit of nm and 1 nm = 10 angstrom. In this instance, the clusters.xpm file contains the cluster assignment map; cluster.log file contains the cluster statistics; cluster_centers.pdb file contains the representative structures; and the rmsd_dist.xvg file plots the RMSD distribution graph. Found 4 clusters.
        gmx cluster -s protein.tpr -f whole_protein_only_900_1000ns.xtc -n tm_helices_ecl.ndx -method gromos -cutoff 0.175 -o whole_protein_clusters_0.175_900_1000ns.xpm -g whole_protein_cluster_0.175_900_1000ns.log -cl whole_protein_cluster_centers_0.175_900_1000ns.pdb -dist whole_protein_rmsd_dist_0.175_900_1000ns.xvg # Select group for least squares fit: TM_ecl_helices and select group for RMSD calculation: TM_ecl_helices. Select Protein for output. The cutoff value is in the unit of nm and 1 nm = 10 angstrom. In this instance, the clusters.xpm file contains the cluster assignment map; cluster.log file contains the cluster statistics; cluster_centers.pdb file contains the representative structures; and the rmsd_dist.xvg file plots the RMSD distribution graph. Found 1 cluster only.
        """, language="bash")
        st.write("Note: Using representative MD-derived protein conformations helps account for protein flexibility, which improves docking realism")
        st.write("Note: Clustering analysis was performed using the GROMOS algorithm with an RMSD cutoff of 0.08 nm. A total of 17 conformational clusters were identified from 4001 trajectory frames. The first cluster was highly dominant, containing 3030 structures (~75.7% of the total frames), indicating that the protein predominantly occupied a stable conformational state throughout the simulation. The average RMSD among clustered structures was 0.0847 nm, with an overall RMSD range of 0.040–0.158 nm, suggesting minimal structural deviation and high conformational stability during the MD simulation. Smaller clusters likely represent transient local fluctuations or minor conformational rearrangements rather than major structural transitions.")

        st.write("###")

        st.write("✔️convert both the protein topology and the whole protein-lipid bilayer topology files (.tpr) to a .gro files using gmx editconf")
        st.code("""
        gmx editconf -f protein.tpr -o protein.gro
        gmx editconf -f md_100ns.tpr -o md_100ns.gro
        """, language="bash")
        st.write("Note: Never load the protein.gro file and protein.xtc file as two separate molecules in VMD, otherwise you wont be able to visualize the trajectory simulation of the protein molecule itself, where the protein molecule will just remain static as the time goes on. This is the most common mistake done by VMD users ! A .gro file contains coordinates and atom names, while an .xtc file contains only moving coordinates with no atom identities. If you load them as two separate molecules, VMD shows a static .gro file and an invisible .xtc file. The correct way is to click the New molecule only once to first load the protein.gro file, followed by clicking the Load button to load the protein trajectory file and remember do not click New molecule again to load the protein trajectory file !")

        st.write("###")

        st.write("✔️visualize the trajectory of only the protein model throughout the 900ns-1000ns simulation using both protein.gro and whole_protein_only_900_1000ns.xtc files")

        st.write("###")

        st.write("✔️visualize the trajectory of the whole protein-lipid bilayer model throughout the 900ns-1000ns simulation using both md_100ns.gro and whole_protein_lipid_bilayer_900_1000ns.xtc files")

        st.write("###")

        st.write("✔️visualize the trajectory of the whole protein-lipid bilayer model throughout the whole 1000ns simulation using both md_100ns.gro and md_reimage.xtc files")

        st.write("###")

        st.write("✔️color and label the TM1-TM7 region of the 7 clusters of the simulated AlphaFold3-predicted CcOctB2R model using UCSF ChimeraX")
        st.code("""
        color #1.1-7:31-54 red target c # color the clustered TM1 region (cartoon representation) as red
        2dlabel create TM1 text "TM1" color red size 15
        color #1.1-7:66-87 light sea green target c # color the clustered TM2 region (cartoon representation) as light sea green
        2dlabel create TM2 text "TM2" color light sea green size 15
        color #1.1-7:104-125 orange red target c # color the clustered TM3 region (cartoon representation) as orange red
        2dlabel create TM3 text "TM3" color orange red size 15
        color #1.1-7:146-167 green target c # color the clustered TM4 region (cartoon representation) as green
        2dlabel create TM4 text "TM4" color green size 15
        color #1.1-7:195-215 blue target c # color the clustered TM5 region (cartoon representation) as blue
        2dlabel create TM5 text "TM5" color blue size 15
        color #1.1-7:279-298 violet target c # color the clustered TM6 region (cartoon representation) as violet
        2dlabel create TM6 text "TM6" color violet size 15
        color #1.1-7:311-334 indigo target c # color the clustered TM7 region (cartoon representation) as indigo
        2dlabel create TM7 text "TM7" color indigo size 15
        color #1.1-7:188,194 lime target c # color the two predicted ligand-binding residues, 188 and 194 (cartoon representation) as gray
        2dlabel create ECL text "ECL" color lime size 15
        2dlabel delete # delete all the labels at once
        """, language="bash")

        st.write("###")

        st.write("**23. set grid box using AutoDock MGL Tools, followed by installing and running AutoDock Vina on Linux to dock the same set of 10 ligands to the only one CPBOctB2R cluster**")
        st.write("✔️determine the CPU architecture of your computer, followed by downloading both the vina_1.2.7_linux_x86_64 and vina_split_1.2.7_linux_x86_64 release packages and making both the packages executable")
        st.code("""
        uname -m # output the CPU architecture of your computer. In my case, it returned x86_64, which indicates 64-bit.
        chmod +x vina_1.2.7_linux_x86_64 # make the release package of AutoDock Vina executable
        chmod +x vina_split_1.2.7_linux_x86_64
        ./vina_1.2.7_linux_x86_64 --help # output the command-line usage manual of AutoDock Vina. Remember to navigate to the working directory, /media/raid/Wee/WeeYeZhi/output/autodockvina_installation, first before you run this command
        ./vina_1.2.7_linux_x86_64 --help_advanced # output the command-line usage manual of AutoDock Vina with advanced option. 
        ./vina_1.2.7_linux_x86_64 --version
        ./vina_split_1.2.7_linux_x86_64 --help # this release package can be used to split the PDBQT ligands
        ./vina_split_1.2.7_linux_x86_64 --version
        """, language="bash")
        st.write("✔️load the refined, energy-minimized AlphaFold3-predicted CcOctB2R into the AutoDock MGL Tools interface, followed by selecting all the ligand binding residues and setting a grid box to cover all of the selected residues to output a grid dimension file containing the size and xyz coordinates of the grid box")
        st.code("""
        select site, resi 105+106+109+110+113+188+190+194+195+198+199+202+290+293+294+297+313+317+321
        print(cmd.centerofmass("site")) # you obtained the xyz coordinates of the grid box [-0.6753963792821812, 3.5131539825982245, -8.578288930509114], where these xyz coordinates of the grid box can be straight away used as an input for both AutoDock4 and AutoDock Vina
        print(cmd.get_extent("site")) # you obtained the xyz minimum and maximum sizes of the grid box [[-10.736000061035156, -7.1539998054504395, -19.104000091552734], [12.78600025177002, 12.302000045776367, 3.3350000381469727]]. Base Size Formula:Size(base) = Coordinate(max) - Coordinate(min). Vina Size Formula (with Buffer):Size(Vina) = Size(base) + Buffer(5 angstrom). Grid Points = Size(Vina)/Grid Spacing(angstrom). 
        """, language="bash")
        st.write("✔️perform molecular docking via AutoDock Vina")
        st.code("""
        mkdir -p Conformer3D_COMPOUND_CID_1480785 Conformer3D_COMPOUND_CID_18526103 Conformer3D_COMPOUND_CID_19606232 Conformer3D_COMPOUND_CID_255273 Conformer3D_COMPOUND_CID_26451 Conformer3D_COMPOUND_CID_26752 Conformer3D_COMPOUND_CID_4184 Conformer3D_COMPOUND_CID_4581 Conformer3D_COMPOUND_CID_577782 Conformer3D_COMPOUND_CID_76145148 Conformer3D_COMPOUND_CID_36324 Conformer3D_COMPOUND_CID_36326 # create a list of ligand directories first before running AutoDock Vina to organize the docking results systematicaly
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_36324.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_36324/Conformer3D_COMPOUND_CID_36324_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_36324/Conformer3D_COMPOUND_CID_36324_CcOctB2R.log 2>&1 & # run this command in the working directory containing the installed executable AutoDock Vina package, /media/raid/Wee/WeeYeZhi/output/autodockvina_installation
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_36324/Conformer3D_COMPOUND_CID_36324_CcOctB2R.pdbqt
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_36326.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_36326/Conformer3D_COMPOUND_CID_36326_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_36326/Conformer3D_COMPOUND_CID_36326_CcOctB2R.log 2>&1 & # run this command in the working directory containing the installed executable AutoDock Vina package, /media/raid/Wee/WeeYeZhi/output/autodockvina_installation
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_36326/Conformer3D_COMPOUND_CID_36326_CcOctB2R.pdbqt
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_1480785.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_1480785/Conformer3D_COMPOUND_CID_1480785_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_1480785/Conformer3D_COMPOUND_CID_1480785_CcOctB2R.log 2>&1 & # run this command in the working directory containing the installed executable AutoDock Vina package, /media/raid/Wee/WeeYeZhi/output/autodockvina_installation
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_1480785/Conformer3D_COMPOUND_CID_1480785_CcOctB2R.pdbqt
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_18526103.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_18526103/Conformer3D_COMPOUND_CID_18526103_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_18526103/Conformer3D_COMPOUND_CID_18526103_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_18526103/Conformer3D_COMPOUND_CID_18526103_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_19606232.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_19606232/Conformer3D_COMPOUND_CID_19606232_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_19606232/Conformer3D_COMPOUND_CID_19606232_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_19606232/Conformer3D_COMPOUND_CID_19606232_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_255273.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_255273/Conformer3D_COMPOUND_CID_255273_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_255273/Conformer3D_COMPOUND_CID_255273_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_255273/Conformer3D_COMPOUND_CID_255273_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_26451.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_26451/Conformer3D_COMPOUND_CID_26451_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_26451/Conformer3D_COMPOUND_CID_26451_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_26451/Conformer3D_COMPOUND_CID_26451_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_26752.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_26752/Conformer3D_COMPOUND_CID_26752_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_26752/Conformer3D_COMPOUND_CID_26752_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_26752/Conformer3D_COMPOUND_CID_26752_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_4184.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_4184/Conformer3D_COMPOUND_CID_4184_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_4184/Conformer3D_COMPOUND_CID_4184_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_4184/Conformer3D_COMPOUND_CID_4184_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_4581.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_4581/Conformer3D_COMPOUND_CID_4581_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_4581/Conformer3D_COMPOUND_CID_4581_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_4581/Conformer3D_COMPOUND_CID_4581_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_577782.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_577782/Conformer3D_COMPOUND_CID_577782_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_577782/Conformer3D_COMPOUND_CID_577782_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_577782/Conformer3D_COMPOUND_CID_577782_CcOctB2R.pdbqt 
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs/whole_protein_cluster_centers_0.175_900_1000ns.pdbqt --ligand /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results/Conformer3D_COMPOUND_CID_76145148.pdbqt --scoring vinardo --center_x -0.675 --center_y 3.513 --center_z -8.578 --size_x 29.25 --size_y 25.5 --size_z 29.25 --spacing 0.375 --seed 12345 --out /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_76145148/Conformer3D_COMPOUND_CID_76145148_CcOctB2R.pdbqt --cpu 16 --exhaustiveness 64 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_76145148/Conformer3D_COMPOUND_CID_76145148_CcOctB2R.log 2>&1 &
        ./vina_split_1.2.7_linux_x86_64 --input /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/simulation_900_1000ns/Conformer3D_COMPOUND_CID_76145148/Conformer3D_COMPOUND_CID_76145148_CcOctB2R.pdbqt 
        """, language="bash")
        st.write("✔️Alternatively, you can perform molecular docking via AutoDock Vina in batch mode by running only one command. However, the generated log file does not show the binding affinity values for each ligand clearly, thats why you chose to run one by one instead in the end")
        st.code("""
        nohup ./vina_1.2.7_linux_x86_64 --receptor /media/raid/Wee/WeeYeZhi/output/yasara_energy_minimization_results/alphafold3_prediction/energy_minimized_alphafold3_predicted_CPB_octopamine_beta2_receptor.pdbqt --batch /media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results --scoring vinardo --center_x 6.263 --center_y -13.193 --center_z 3.302 --size_x 50 --size_y 60 --size_z 50 --spacing 0.375 --dir /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results --cpu 16 --seed 12345 --exhaustiveness 32 --energy_range 5 --num_modes 15 > /media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results/autodock_vina_results.log 2>&1 &
        """, language="bash")
        st.markdown("[Visit the latest GitHub release of AutoDock Vina for installation](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7)")
        st.markdown("[Visit the AutoDock Vina Tutorial Webpage](https://autodock-vina.readthedocs.io/en/latest/docking_requirements.html)")
        st.write("Note: While setting the grid box, make sure the grid box points size being set for AutoDock4 is even integers. You can add a 5 angstrom buffer to the grid box size after running the Pymol commands to enable the rotation and conformational changes/search of the ligands")
        st.write("Molecular docking simulations were performed using AutoDock Vina (version 1.2.7) employing the Vinardo scoring function to evaluate ligand-receptor binding affinities. The search space was centered on the geometric center of mass of the predefined binding pocket residues, establishing grid center coordinates at X = -0.675, Y = 3.513, and Z = -8.578. To define the search volume, the Cartesian extrema of the selected residues were calculated, and an approximate 5 Å buffer was added to ensure unrestricted conformational sampling of the ligand. These dimensions were optimized to align with a standard 0.375 Å grid spacing, resulting in a final grid box size of 29.25 × 25.50 × 29.25 Å. To achieve thorough conformational sampling, the search algorithm was executed with an elevated exhaustiveness parameter of 64. The simulation was configured to generate a maximum of 15 binding poses (num_modes = 15) with a maximum energy penalty of 5 kcal/mol relative to the top-ranked pose (energy_range = 5). All calculations were parallelized across 16 CPU cores, and a constant random seed (12345) was applied to ensure full reproducibility of the docking results.")

        st.write("###")

        st.write("**24. Predict the 3D structure of the CPBOctB2R_octopamine_complex via AlphaFold3 locally**")
        st.code("""
        docker run --user 0:0 -d --name alphafold3_prediction_of_CPBOctB2R_octopamine_complex -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/input:/root/af_input -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/output:/root/af_output -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/alphafold3_model_parameters:/root/models -v /media/raid/Wee/WeeYeZhi/output/alphafold3_prediction_results/downloaded_alphafold3_public_databases:/root/public_databases --gpus all alphafold3 python run_alphafold.py --json_path=/root/af_input/CPBOctB2R_octopamine_complex/CPBOctB2R_octopamine_complex.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output/CPBOctB2R_octopamine_complex
        """, language="bash")
        # ----LOAD ALPHAFOLD3 CPBOctB2R_octopamine_complex JSON FILE----
        # Check if the file exists before reading
        if CPBOctB2R_octopamine_complex.exists():
            with open(CPBOctB2R_octopamine_complex, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download CPBOctB2R_octopamine_complex JSON File",
                data=script_byte,
                file_name=CPBOctB2R_octopamine_complex.name,  # Extract just the file name
                mime="application/json",  # MIME type for shell scripts
            )
        else:
            st.error(f"{CPBOctB2R_octopamine_complex.name} does not exist.")
        st.write("✔️perform structural superimposition between the AlphaFold3-predicted CPBOctB2R and AlphaFold3-predicted CPBOctB2R-octopamine complex")
        st.code("""
        open CPB_OctB2R_Octopamine_seed-1_sample-1_model.cif # This will be #1 model
        open CPB_OctB2R_seed-1_sample-2_model.cif # This will be your #2 model
        matchmaker #1 to #2 # The refined RMSD value is 0.424 angstroms
        """, language="bash")

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

# Submit Transcriptome Assembly to ENA

if selected == "Submit Transcriptome Assembly to ENA":
    with st.container():
        st.write("---")
        st.header("Time to submit CPB transcriptome assembly to the European Nucleotide Archive (ENA) using either Webin-CLI Docker or Webin-CLI Java 🤖")

        st.write("###")

        st.write("**1. Download the manifest text file and the gzipped CPB transcriptome assembly fasta file to be used for validation and submission to the ENA database.**")
        # ----LOAD MANIFEST FILE OF CPB TRANSCRIPTOME ASSEMBLY----
        # Check if the file exists before reading
        if ENA_manifest.exists():
            with open(ENA_manifest, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download Manifest Text File",
                data=script_byte,
                file_name=ENA_manifest.name,  # Extract just the file name
                mime="application/x-sh",  # MIME type for shell scripts
            )
        else:
            st.error(f"{ENA_manifest.name} does not exist.")
        # ----LOAD CPB TRANSCRIPTOME ASSEMBLY FASTA FILE----
        # Check if the file exists before reading
        if CPB_transcriptome_assembly.exists():
            with open(CPB_transcriptome_assembly, "rb") as script_file:
                script_byte = script_file.read()

            # Add download button
            st.download_button(
                label="Download CPB Transcriptome Assembly Fasta File",
                data=script_byte,
                file_name=CPB_transcriptome_assembly.name,  # Extract just the file name
                mime="application/gzip",  # MIME type for shell scripts
            )
        else:
            st.error(f"{CPB_transcriptome_assembly.name} does not exist.")

        st.write("###")

        st.write("**2. Run the Webin-CLI program using Java**")
        st.write("✔️download the latest release of the Webin-CLI jar file (webin-cli-9.0.3.jar) manually from the official Webin-CLI GitHub page and save the file inside your desired working Webin-CLI installation directory")
        st.code("""
        java -jar /media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/webin-cli-9.0.3.jar -help
        java -jar /media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/webin-cli-9.0.3.jar -version
        """, language="bash")
        st.write("✔️validate both the manifest file and the CPB transcriptome assembly fasta file to ENA using the Webin-CLI submission service (remember to use your own mobile hotspot network to run both the validation and submission process of the transcriptome assembly file)")
        st.code("""
        java -jar /media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/webin-cli-9.0.3.jar -validate -centerName='Universiti Kebangsaan Malaysia' -context=transcriptome -inputDir=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/input -manifest=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/input/manifest.txt -outputDir=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/output -passwordFile=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/input/password.txt -userName=Webin-53180
        """, language="bash")
        st.write("✔️submit both the manifest file and the CPB transcriptome assembly fasta file to ENA using the Webin-CLI submission service")
        st.code("""
        java -jar /media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/webin-cli-9.0.3.jar -submit -centerName='Universiti Kebangsaan Malaysia' -context=transcriptome -inputDir=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/input -manifest=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/input/manifest.txt -outputDir=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/output -passwordFile=/media/raid/Wee/WeeYeZhi/output/ENA_transcriptome_submission/input/password.txt -userName=Webin-53180
        """, language="bash")

        st.write("###")

        st.write("**3. Alternatively, run the Webin-CLI program using Docker**")
        st.write("✔️pull the Webin-CLI Docker image from DockerHub")
        st.code("""
        docker pull enasequence/webin-cli # you are pulling the latest enasequence/webin-cli docker image from DockerHub
        """, language="bash")
        st.write("✔️obtain the digest ID of the Webin-CLI docker image")
        st.code("""
        docker image inspect --format '{{.RepoDigests}}' enasequence/webin-cli:latest # you will obtain enasequence/webin-cli@sha256:62fbfba3aaa3172d51beb99a102427d6d31ed7598deb167e2aaa18521e8924d3 where other scientists can run the command, docker pull trinityrnaseq/transdecoder@sha256:80fa372fe94bc0ac4440fea7905c47e05187d13b1e837f8dca72c0ca49ab1bb4 to use the same docker transdecoder image with the same software environment
        """, language="bash")
        st.write("✔️check whether the Webin-CLI docker image has been successfully pulled")
        st.code("""
        docker images
        """, language="bash")
        st.write("✔️output the command-line usage manual of the Webin-CLI Docker Image")
        st.code("""
        docker run --rm enasequence/webin-cli -help
        docker run --rm enasequence/webin-cli -version
        """, language="bash")
        st.write("✔️run the docker container to validate both the manifest file and the CPB transcriptome assembly fasta file to ENA using the Webin-CLI submission service")
        st.code("""
        docker run --rm --user 1000:1000 -v /media/raid/Wee/WeeYeZhi/output:/data enasequence/webin-cli -validate -centerName='Universiti Kebangsaan Malaysia' -context=transcriptome -inputDir=/data/ENA_transcriptome_submission/input -manifest=/data/ENA_transcriptome_submission/input/manifest.txt -outputDir=/data/ENA_transcriptome_submission/output -passwordFile=/data/ENA_transcriptome_submission/input/password.txt -userName=Webin-53180
        """, language="bash")
        st.write("✔️run the docker container to submit both the manifest file and the CPB transcriptome assembly fasta file to ENA using the Webin-CLI submission service")
        st.code("""
        docker run --rm --user 1000:1000 -v /media/raid/Wee/WeeYeZhi/output:/data enasequence/webin-cli -submit -centerName='Universiti Kebangsaan Malaysia' -context=transcriptome -inputDir=/data/ENA_transcriptome_submission/input -manifest=/data/ENA_transcriptome_submission/input/manifest.txt -outputDir=/data/ENA_transcriptome_submission/output -passwordFile=/data/ENA_transcriptome_submission/input/password.txt -userName=Webin-53180
        """, language="bash")
        st.markdown("[Read how to submit transcriptome assembly fasta file to ENA database](https://ena-docs.readthedocs.io/en/latest/submit/assembly/transcriptome.html)")
        st.markdown("[Read how to run Webin-CLI programme via Java and Docker](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html)")
        st.markdown("[Obtain the latest version of Webin-CLI Java from GitHub](https://github.com/enasequence/webin-cli/releases/tag/)")
        st.markdown("[Visit the Webin-CLI Docker Image in the DockerHub](https://hub.docker.com/r/enasequence/webin-cli)")
        st.markdown("[Visit the latest version of Aspera-CLI Official GitHub Page](https://github.com/IBM/aspera-cli/releases/tag/v4.26.1)")
        st.write("Note: Avoid using your institution's computer to validate and submit assembly files as it has firewall issues which interrupts the assembly file upload process to the ENA database via the FTP server later (due to possible internet reconnection issue by peers), where the complete file upload process will always fail later. Instead, try using your own home laptop, home wifi or your own mobile hotspot network to make the assembly submission process to bypass the firewall issue. Bear in mind that genome and transcriptome assemblies can only be submitted using the Webin-CLI submission interface. Remember to gzip your transcriptome assembly fasta file first before making submission. You have to run the validate command first to validate the content of the files, followed by running the submit command to submit the Webin-CLI-validated files to the ENA database.")

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

        st.write("30. Use the keyboard shortcut, Ctrl + Y, to remove a specific line of code within the PyCharm IDE and use the keyboard shortcut, Ctrl + C and Ctrl + V, to copy and paste the line of code in PyCharm IDE")

        st.write("###")

        st.write("###")