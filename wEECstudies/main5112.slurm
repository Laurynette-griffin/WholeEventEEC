#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=2

#SBATCH --mem=32G

#SBATCH --time=04:00:00

#SBATCH --job-name=wholeevent_EEC_ee

#SBATCH --output=../offline/pythia_wholeeventEEC_ee_starEEC_Feb21.root

#SBATCH --output=logs/main5112.log

#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/examples

./main5112 main5112.cmnd 

echo "Job Finished"