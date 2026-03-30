#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --job-name=weec_hydrotest
#SBATCH --array=1-10               # Submit 10 jobs in the array
#SBATCH --output=logs/weec_hydrotest%A_%a.log
#SBATCH --account=rke_group

# Load environment
source /home/griffl6/setup_pythia.sh

# Navigate to working directory
cd /home/griffl6/work/pythia/wEECstudies

# Get current date 
DATE=$(date +%Y-%m-%d)

EXECUTABLE_NAME="./wEEC_hydro"

OUTPUT_FILE="/home/griffl6/work/pythia/paper_outfiles/pythia_pp_rhic_fullevent_5pthat60_CentralityStudy010_psi2match_trackptcut1gev_copycopy_${DATE}_${SLURM_ARRAY_TASK_ID}.root"

# Input Pythia TTree File 5M events ran 5k events passed cut
PYTHIA_FILE="/home/griffl6/work/pythia/wEECstudies/outfiles/ttreepythia_pt1_31_40_pt2_20_31_pthat15_60_b2b_3pi_4_chargedevent_event_sphenixcuts_singlejet_20_30_escheme_wta_hardcone_Feb23_5k.root"

CONFIG_FILE="wEEC_hydro.cmnd"
# List of Thermal ASCII Files (10 files) !!!
THERMAL_BASE_DIR="/data/rke_group/AuAuOnlyHydroFiles/Simulations/Centrality0_10"
THERMAL_FILES=()
for i in {000..099}; do
    # This constructs the full path for each of your i thermal files
    THERMAL_FILES+=("${THERMAL_BASE_DIR}/job_AuAu_only_soft_sector_hydro_10_events_centr0_10_test_${i}/event_0_10_final_state_hadrons.dat")
done

TOTAL_EVENTS=4900 
NUM_JOBS=10

# Calculate chunk size (integer division)
CHUNK_SIZE=$(( TOTAL_EVENTS / NUM_JOBS ))

# Calculate start and end indices for the current job
START_EVENT=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

# The end event is calculated for the next job's start, or the total events for the last job
if [ "${SLURM_ARRAY_TASK_ID}" -eq $(( NUM_JOBS - 1 )) ]; then
    END_EVENT=${TOTAL_EVENTS}
else
    END_EVENT=$(( (SLURM_ARRAY_TASK_ID + 1) * CHUNK_SIZE ))
fi

# Join the thermal files into a space-separated string for the command line
THERMAL_FILE_LIST="${THERMAL_FILES[*]}"
    
    "${EXECUTABLE_NAME}" \
    "${CONFIG_FILE}" \
    "${OUTPUT_FILE}" \
    "${PYTHIA_FILE}" \
    "${START_EVENT}" \
    "${END_EVENT}" \
    ${THERMAL_FILE_LIST}

echo "Job finished."