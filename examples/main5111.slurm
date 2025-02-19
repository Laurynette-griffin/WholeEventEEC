#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu

#SBATCH --mail-type=ALL

#SBATCH --cpus-per-task=4

#SBATCH --time=04:00:00

#SBATCH --job-name=wholeevent_EEC

#SBATCH --mem=32G

#SBATCH --output=pythia_wholeeventEEC_starEEC_Feb17.out

#SBATCH --account=rke_group

export LD_LIBRARY_PATH=/data/rke_group/Software/LHAPDF6/local/lib:$LD_LIBRARY_PATH

export PATH=/data/rke_group/Software/PYTHIA/pythia8312_local

main5111 main511.cmnd 

echo "Job Finished"