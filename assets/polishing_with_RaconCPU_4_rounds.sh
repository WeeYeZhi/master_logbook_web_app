#!/bin/bash

# Align PacBio reads to assembly for the first time
echo "Starting round 1: Align PacBio reads to initial assembly with minimap2..."
minimap2 -x map-pb -t 24 /media/raid/Wee/WeeYeZhi/output/MaSuRCA_results/MaSuRCA_raw_assembly_latestmodifiedPE_gapclosing_trimmedadapter_results/CA.mr.67.17.15.0.02/primary.genome.scf.fasta /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_1/long_read_alignment_round1.paf
echo "Minimap2 alignment round 1 completed."

# Polish assembly with racon for the first time
echo "Starting round 1: Polishing assembly with racon..."
racon -t 24 /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_1/long_read_alignment_round1.paf /media/raid/Wee/WeeYeZhi/output/MaSuRCA_results/MaSuRCA_raw_assembly_latestmodifiedPE_gapclosing_trimmedadapter_results/CA.mr.67.17.15.0.02/primary.genome.scf.fasta > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_1/polished1.fasta
echo "Racon polishing round 1 completed."

# Align PacBio reads to assembly for the second time
echo "Starting round 2: Align PacBio reads to polished assembly with minimap2..."
minimap2 -x map-pb -t 24 /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_1/polished1.fasta /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_2/long_read_alignment_round2.paf
echo "Minimap2 alignment round 2 completed."

# Polish assembly with racon for the second time
echo "Starting round 2: Polishing assembly with racon..."
racon -t 24 /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_2/long_read_alignment_round2.paf /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_1/polished1.fasta > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_2/polished2.fasta
echo "Racon polishing round 2 completed."

# Align PacBio reads to assembly for the third time
echo "Starting round 3: Align PacBio reads to polished assembly with minimap2..."
minimap2 -x map-pb -t 24 /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_2/polished2.fasta /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_3/long_read_alignment_round3.paf
echo "Minimap2 alignment round 3 completed."

# Polish assembly with racon for the third time
echo "Starting round 3: Polishing assembly with racon..."
racon -t 24 /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_3/long_read_alignment_round3.paf /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_2/polished2.fasta > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_3/polished3.fasta
echo "Racon polishing round 3 completed."

# Align PacBio reads to assembly for the fourth time
echo "Starting round 4: Align PacBio reads to polished assembly with minimap2..."
minimap2 -x map-pb -t 24 /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_3/polished3.fasta /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_4/long_read_alignment_round4.paf
echo "Minimap2 alignment round 4 completed."

# Polish assembly with racon for the fourth time
echo "Starting round 4: Polishing assembly with racon..."
racon -t 24 /media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_4/long_read_alignment_round4.paf /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_3/polished3.fasta > /media/raid/Wee/WeeYeZhi/output/racon-CPU_results/racon_polishing_4/polished4.fasta
echo "Racon polishing round 4 completed."

echo "All 4 polishing rounds finished successfully!"

