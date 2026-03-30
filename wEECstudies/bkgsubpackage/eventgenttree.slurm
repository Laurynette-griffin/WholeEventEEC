#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --job-name=eventgen
#SBATCH --array=1-100 # Submit jobs in the array
#SBATCH --output=logs/eventgen_%A_%a.log
#SBATCH --account=rke_group

# Load environment
source /home/griffl6/setup_pythia.sh

# Navigate to working directory
cd /home/griffl6/work/pythia/wEECstudies/bkgsubpackage

# Get current date 
DATE=$(date +%Y-%m-%d)

# unique output file name
OUTPUT_FILE="/home/griffl6/work/pythia/outfiles/ttrees/pythia_pt1_31_40_pt2_20_31_pthat15_60_b2b_3pi_4_chargedevent_event_sphenixcuts_singlejet_20_30_escheme_wta_hardcone_${DATE}_${SLURM_ARRAY_TASK_ID}.root"

# Run with output file as argument
./eventgenttree eventgenttree.cmnd $OUTPUT_FILE

echo "Job $SLURM_ARRAY_TASK_ID finished"
