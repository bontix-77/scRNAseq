#!/bin/bash
#SBATCH --job-name=Seurat
#SBATCH --output=convert_%j.out
#SBATCH --error=convert_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1       # The old tool is single-threaded
#SBATCH --mem=64G                # Uses less RAM
#SBATCH --time=24:00:00         # Needs MORE time (it is slower)

module load nextflow
module load singularity

# Define paths
WORK_DIR="/home/abontemp/nf-scRNA"


cd "$WORK_DIR"

echo "Starting Seurat"
nextflow run main_HPC.nf -offline

echo "Done."