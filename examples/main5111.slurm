#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=4

#SBATCH --mem=64G

#SBATCH --time=05:00:00

#SBATCH --job-name=45wholeevent_EEC

#SBATCH --output=../offline/pythia_wholeeventEEC_15GeV_starEEC_10pthat60_Feb21.root

#SBATCH --output=logs/main5111.log

#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/examples

./main5111 main5111.cmnd 

echo "Job Finished"