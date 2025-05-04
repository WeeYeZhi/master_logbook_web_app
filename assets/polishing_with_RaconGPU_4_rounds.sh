#!/bin/bash

READS="/media/raid/Wee/WeeYeZhi/resources_from_LKM/raw_pacbio_read/PacBio.fq" # Your long reads file
INITIAL_ASSEMBLY="/media/raid/Wee/WeeYeZhi/output/MaSuRCA_results/MaSuRCA_raw_assembly_latestmodifiedPE_gap_closing_trimmedadapter_results/CA.mr.67.17.15.0.02/primary.genome.scf.fasta" # Your initial assembly file
ROUNDS=4    # Number of polishing rounds
THREADS=32  # Number of CPU threads to use

# Set up working directory
WORKDIR="/media/raid/Wee/WeeYeZhi/output/racon-GPU_results/racon_polishing"
mkdir -p $WORKDIR
cp $INITIAL_ASSEMBLY $WORKDIR/assembly_round0.fasta

cd $WORKDIR

for (( i=1; i<=$ROUNDS; i++ ))
do
    PREV_ROUND=$((i-1))
    echo "=== Polishing round $i ==="
    
    # Create a round-specific .paf file name
    PAF_FILE="alignment_round${i}.paf"
    
    # Align reads to the current assembly
    minimap2 -t $THREADS -x map-pb assembly_round${PREV_ROUND}.fasta ../$READS > $PAF_FILE
    
    # Run racon-gpu
    racon -t $THREADS ../$READS $PAF_FILE assembly_round${PREV_ROUND}.fasta > assembly_round${i}.fasta

done

echo "All $ROUNDS rounds of polishing complete."
echo "Final polished assembly: $WORKDIR/assembly_round${ROUNDS}.fasta"
