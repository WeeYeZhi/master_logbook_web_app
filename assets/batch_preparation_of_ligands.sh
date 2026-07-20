#!/bin/bash

MGLTOOLS=/media/raid/Wee/WeeYeZhi/output/autodock_mgl_tools_installation/mgltools_x86_64Linux2_1.5.7
INPUT_DIR=/media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/openbabel_results
OUTPUT_DIR=/media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.mol2; do
    
    # Extract filename only (no path, no extension)
    base=$(basename "${file%.mol2}")
    
    # Define output path
    output="$OUTPUT_DIR/${base}.pdbqt"

    echo "Processing: $file"

    $MGLTOOLS/bin/pythonsh \
    $MGLTOOLS/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
    -l "$file" -o "$output" -v -A hydrogens -U nphs

    if [ $? -ne 0 ]; then
        echo "❌ Failed: $file"
    else
        echo "✅ Done: $output"
    fi

done