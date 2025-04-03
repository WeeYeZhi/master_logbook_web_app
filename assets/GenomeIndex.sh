#!/bin/bash

# Define the directory where the genome files of CPB are located
directory="/home/cbr16/Documents/WeeYeZhi/output/GenomeIndexresults/GenomeIndex"

# Create STAR index for the Conopomorpha cramerella reference genome
STAR \
   --runThread 8 \
   --runMode genomeGenerate \
   --genomeDir "$directory/" \
   --genomeFastaFiles "$directory/CPB_insect_draft_assembly.v4.fna" \
   --sjdbGTFfile "$directory/CPB_insect_draft_assembly.v4.gff3" \
   --sjdbOverhang 99