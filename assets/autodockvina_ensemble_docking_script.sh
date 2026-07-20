#!/bin/bash

# AutoDock Vina executable
VINA="./vina_1.2.7_linux_x86_64"

# Directories
RECEPTOR_DIR="/media/raid/Wee/WeeYeZhi/output/charmmgui_results/charmm-gui-7882792398/gromacs"
LIGAND_DIR="/media/raid/Wee/WeeYeZhi/output/latest_ligand_preparation_results/autodock_mgltools_results"
OUTPUT_DIR="/media/raid/Wee/WeeYeZhi/output/molecular_docking_results/autodock_vina_results"

# Docking parameters
CENTER_X=6.263
CENTER_Y=-13.193
CENTER_Z=3.302

SIZE_X=50
SIZE_Y=60
SIZE_Z=50

CPU=16
EXHAUSTIVENESS=32
NUM_MODES=15
ENERGY_RANGE=5
SPACING=0.375

# Loop through every ligand
for ligand in "$LIGAND_DIR"/*.pdbqt
do
    ligand_name=$(basename "$ligand" .pdbqt)

    # Dock against apo_cluster1, apo_cluster2 and apo_cluster3
    for cluster in 1 2 3
    do
        receptor="$RECEPTOR_DIR/apo_cluster${cluster}.pdbqt"

        output="$OUTPUT_DIR/${ligand_name}_apocluster${cluster}.pdbqt"

        logfile="$OUTPUT_DIR/${ligand_name}_apocluster${cluster}.log"

        echo "Docking ${ligand_name} against apo_cluster${cluster}..."

        nohup "$VINA" \
            --receptor "$receptor" \
            --ligand "$ligand" \
            --scoring vinardo \
            --center_x $CENTER_X \
            --center_y $CENTER_Y \
            --center_z $CENTER_Z \
            --size_x $SIZE_X \
            --size_y $SIZE_Y \
            --size_z $SIZE_Z \
            --spacing $SPACING \
            --seed 12345
            --cpu $CPU \
            --exhaustiveness $EXHAUSTIVENESS \
            --energy_range $ENERGY_RANGE \
            --num_modes $NUM_MODES \
            --out "$output" \
            > "$logfile" 2>&1 &
    done
done

wait

echo "All docking jobs submitted."