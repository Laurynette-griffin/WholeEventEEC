#!/bin/bash

#SBATCH --mail-user=laurynette.griffin@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=hadtester        # Job name
#SBATCH --output=logs/array_%A_%a.out  # Standard output file
#SBATCH --error=logs/array_%A_%a.err   # Standard error file
#SBATCH --array=0                   # Create 10 jobs (indices 0 through 9)
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1                     
#SBATCH --mem=32G                      
#SBATCH --time=10:30:00                
#SBATCH --account=rke_group

source /home/griffl6/setup_pythia.sh

# Navigate to working directory where the compiled executable is located
cd /home/griffl6/work/pythia/wEECstudies/bkgsubpackage

# Define the name of your compiled executable (assuming your C++ file is compiled)
EXECUTABLE_NAME="./HadronstoTree"

# Input Pythia TTree File 5M events
PYTHIA_FILE="outfiles/ttreepythia_pt1_31_40_pt2_20_31_pthat15_60_b2b_3pi_4_sphenixcuts_Dec3_5M.root"

# List of Thermal ASCII Files (10 files) !!!
THERMAL_BASE_DIR="/data/rke_group/AuAuOnlyHydroFiles/Simulations/Centrality0_10"
THERMAL_FILES=()
for i in {0..1}; do
    # This constructs the full path for each of your i thermal files
    THERMAL_FILES+=("${THERMAL_BASE_DIR}/job_AuAu_only_soft_sector_hydro_10_events_centr0_10_test_00${i}/event_0_10_final_state_hadrons.dat")
done

# Join the thermal files into a space-separated string for the command line
THERMAL_FILE_LIST="${THERMAL_FILES[*]}"

# Analysis parameters
N_POINT=2  # N-point correlator (adjust if running N=2, etc.)
OUTPUT_BASE_NAME="outfiles/raw_nothadded/eec_bkgsub_Test1_${SLURM_ARRAY_TASK_ID}" # Base name for final ROOT file

# =========================================================================
# EVENT INDEX CALCULATION (for parallelization)
# =========================================================================

# # of events in the Pythia TTree 

TOTAL_EVENTS=1000 
NUM_JOBS=1

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

# =========================================================================
# EXECUTION
# =========================================================================

echo "Running EEC Analysis for N=${N_POINT}"
echo "Pythia File: ${PYTHIA_FILE}"
echo "Total Thermal Files to be loaded: ${#THERMAL_FILES[@]}"
echo "Processing Pythia Events: ${START_EVENT} to ${END_EVENT}"
echo "Output Base Name: ${OUTPUT_BASE_NAME}"
echo "Job Index: ${SLURM_ARRAY_TASK_ID}"

# The C++ executable expects arguments in this order:
# 1. <PythiaROOTFile>
# 2. <N-point>
# 3. <OutputName>
# 4. <StartEvent>
# 5. <EndEvent>
# 6. <ThermalFile1> [ThermalFile2] ...

"${EXECUTABLE_NAME}" \
    "${PYTHIA_FILE}" \
    "${N_POINT}" \
    "${OUTPUT_BASE_NAME}" \
    "${START_EVENT}" \
    "${END_EVENT}" \
    ${THERMAL_FILE_LIST}

echo "Job finished."