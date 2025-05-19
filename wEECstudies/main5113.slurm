#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=pt1gwholeevent
#SBATCH --array=1-50                # Submit 10 jobs in the array
#SBATCH --output=logs/main5113_%A_%a.log
#SBATCH --account=rke_group

# Load environment
source /home/griffl6/setup_pythia8312.sh

# Navigate to working directory
cd /home/griffl6/work/pythia/examples

# Get current date 
DATE=$(date +%Y-%m-%d)

# unique output file name
OUTPUT_FILE="/home/griffl6/work/pythia/outfiles/pythia_pp_rhic_fullevent_pt1_on_${DATE}_${SLURM_ARRAY_TASK_ID}.root"

# Run with output file as argument
./main5113 main5113.cmnd $OUTPUT_FILE

echo "Job $SLURM_ARRAY_TASK_ID finished"