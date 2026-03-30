#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=2

#SBATCH --mem=32G

#SBATCH --time=06:00:00

#SBATCH --job-name=wholeevent_EEC_ee
#SBATCH --array=1-50                # Submit 10 jobs in the array
#SBATCH --output=logs/ee_eec_%A_%a.log
#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/wEECstudies/


# Get current date 
DATE=$(date +%Y-%m-%d)

#output file name
OUTPUT_FILE="/home/griffl6/work/pythia/outfiles/pythia_ee_aleph_fullevent_1mil_${DATE}_${SLURM_ARRAY_TASK_ID}.root"

./AlephwEEC AlephwEEC.cmnd $OUTPUT_FILE

echo "Job Finished"