#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=2

#SBATCH --mem=64G

#SBATCH --time=05:00:00

#SBATCH --job-name=bothoffwholeevent_EEC_45

#SBATCH --output=../offline/pythia_pp_rhic_45GeV_fullevent_10pthat60_Apr9_bothoff.root

#SBATCH --output=logs/main5111bothoff.log

#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/examples

./main5111bothoff main5111bothoff.cmnd 

echo "Job Finished"