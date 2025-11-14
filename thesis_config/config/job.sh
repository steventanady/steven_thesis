#!/bin/bash
#SBATCH --job-name=cosmosis_run
#SBATCH -n 5
#SBATCH --partition=cpu
#SBATCH --mem-per-cpu=3000
#SBATCH --time=03:00:00
#SBATCH --output=cosmosis_%j.out
#SBATCH --error=cosmosis_%j.err
#SBATCH --mail-type=ALL

echo "Starting SLURM job on $(hostname)"
echo "Current PATH: $PATH"
echo "Loading Conda..."

# Activate Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ~/env || { echo "Conda env failed to activate"; exit 1; }

# Debug: Check if cosmosis exists
which cosmosis || { echo "CosmoSIS not found in env"; exit 127; }

# Source CosmoSIS setup
echo "Sourcing CosmoSIS configure..."
source cosmosis-configure || { echo "Failed to source cosmosis-configure"; exit 1; }

# Change to config directory
cd ~/thesis/thesis_config/config/tester || { echo "Directory not found"; exit 1; }

# Run CosmoSIS
echo "Running CosmoSIS..."
mpirun -n 5 cosmosis --mpi tester_oms8.ini
