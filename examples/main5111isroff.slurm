#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=2

#SBATCH --mem=32G

#SBATCH --time=05:00:00

#SBATCH --job-name=isroffwholeevent_EEC_45

#SBATCH --output=../offline/pythia_wholeeventEEC_15GeV_starEEC_10pthat60_Feb21.root

#SBATCH --output=logs/main5111isroff.log

#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/examples

./main5111isroff main5111isroff.cmnd 

echo "Job Finished"