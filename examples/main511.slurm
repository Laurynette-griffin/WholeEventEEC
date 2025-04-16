#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=4

#SBATCH --mem=64G

#SBATCH --time=04:00:00

#SBATCH --job-name=8wholeevent

#SBATCH --output=../offline/pythia_pp_rhic_8GeV_fullevent_5pthat60_Apr9_dr_thirdjet.root

#SBATCH --output=logs/main511.log

#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/examples

./main511 main511.cmnd 

echo "Job Finished"