#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=pt2variable
#SBATCH --array=1-50 # Submit jobs in the array
#SBATCH --output=logs/pt2variable_%A_%a.log
#SBATCH --account=rke_group

# Load environment
source /home/griffl6/setup_pythia.sh

# Navigate to working directory
cd /home/griffl6/work/pythia/wEECstudies

# Get current date 
DATE=$(date +%Y-%m-%d)

#output file name
OUTPUT_FILE="/home/griffl6/work/pythia/outfiles/pythia_pp_rhic_fullevent_15pthat70_pt2variable_1mil_${DATE}_${SLURM_ARRAY_TASK_ID}.root"

# Run with output file as argument
./CentralityStudy2530 CentralityStudypt1.cmnd $OUTPUT_FILE

echo "Job $SLURM_ARRAY_TASK_ID finished"
