#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --job-name=ctztheta
#SBATCH --array=1-20               # Submit 100 jobs in the array
#SBATCH --output=logs/clzt_%A_%a.log
#SBATCH --account=rke_group

# Load environment
source /home/griffl6/setup_pythia.sh

# Navigate to working directory
cd /home/griffl6/work/pythia/wEECstudies

# Get current date 
DATE=$(date +%Y-%m-%d)

# unique output file name
OUTPUT_FILE="/home/griffl6/work/pythia/outfiles/pythia_pp_rhic_fullevent_5pthat60_closuretrialztheta_100k__${DATE}_${SLURM_ARRAY_TASK_ID}.root"

# Run with output file as argument
./closuretrialztheta centralitymbs02.cmnd $OUTPUT_FILE

echo "Job $SLURM_ARRAY_TASK_ID finished"
