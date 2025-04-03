#!/bin/bash

# Ensure conda is initialized
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed or not in PATH. Please install Conda first."
    exit 1
fi

# Create and activate a new Conda environment (optional, change 'braker_env' as needed)
env_name="braker_env"
echo "Creating Conda environment: $env_name"
conda create -y -n $env_name python=3.9
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $env_name

# Install dependencies
modules=(
    "anaconda perl"
    "anaconda biopython"
    "bioconda perl-app-cpanminus"
    "bioconda perl-file-spec"
    "bioconda perl-hash-merge"
    "bioconda perl-list-util"
    "bioconda perl-module-load-conditional"
    "bioconda perl-posix"
    "bioconda perl-file-homedir"
    "bioconda perl-parallel-forkmanager"
    "bioconda perl-scalar-util-numeric"
    "bioconda perl-yaml"
    "bioconda perl-class-data-inheritable"
    "bioconda perl-exception-class"
    "bioconda perl-test-pod"
    "bioconda perl-file-which"
    "bioconda perl-mce"
    "bioconda perl-threaded"
    "bioconda perl-list-util"
    "bioconda perl-math-utils"
    "bioconda cdbtools"
    "eumetsat perl-yaml-xs"
    "bioconda perl-data-dumper"
)

echo "Installing Perl dependencies..."
for module in "${modules[@]}"; do
    echo "Installing: $module"
    conda install -y -c $module
    if [ $? -ne 0 ]; then
        echo "Failed to install: $module"
        exit 1
    fi
done

echo "All dependencies installed successfully."
