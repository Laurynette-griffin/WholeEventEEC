#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --nodes=1

#SBATCH --cpus-per-task=2

#SBATCH --mem=32G

#SBATCH --time=05:00:00

#SBATCH --job-name=mpioffwholeevent_EEC_45

#SBATCH --output=../offline/pythia_pp_rhic_45GeV_fullevent_10pthat60_Apr9_mpioff.root

#SBATCH --output=logs/main5111mpioff.log

#SBATCH --account=rke_group

source /home/griffl6/setup_pythia8312.sh

cd /home/griffl6/work/pythia/examples

./main5111mpioff main5111mpioff.cmnd 

echo "Job Finished"