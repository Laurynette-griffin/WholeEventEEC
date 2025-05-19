#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --job-name=ThermalPhiStudy
#SBATCH --array=1 # Submit jobs in the array
#SBATCH --output=logs/ThermalPhiStudy_%A_%a.log
#SBATCH --account=rke_group

# Load environment
source /home/griffl6/setup_pythia.sh

# Navigate to working directory
cd /home/griffl6/work/wEECstudies

# Get current date 
DATE=$(date +%Y-%m-%d)

# unique output file name
OUTPUT_FILE="/home/griffl6/work/wEECstudies/outfiles/thermalbackgroundphistudy_10k_${DATE}_${SLURM_ARRAY_TASK_ID}.root"

# Run with output file as argument
./ThermalPhiStudy $OUTPUT_FILE

echo "Job $SLURM_ARRAY_TASK_ID finished"