#!/bin/bash
#SBATCH -J mesa_grid
#SBATCH -p cca
#SBATCH --array=0-17            # 18 models: indices 0–17 for masses 15–100
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

# Load MESA + OpenBLAS
source ~/.bashrc
mesa-12778
module load openblas 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Base directories
BASE="/mnt/home/lvanson/ASBE-course/practicum/session3/solutions_session3/data"
MESA_BASE="$BASE/HW3_base"
MESA_INLIST="$MESA_BASE/inlist"
OUTPUT_DIR="$BASE/HW3_output"

# Define your mass grid (must match create_grid_run.py)
MASSES=(15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100)
DUTCH=1

# Get this job’s mass from array index
MASS=${MASSES[$SLURM_ARRAY_TASK_ID]}
RUN_DIR="$OUTPUT_DIR/M${MASS}_DW${DUTCH}"

# Sanity check
echo "Starting MESA run for M=${MASS} Msun in $RUN_DIR"
echo "Running on $(hostname) at $(date)"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

cd "$RUN_DIR" || { echo "Error: directory $RUN_DIR not found"; exit 1; }

# Run MESA
/usr/bin/time -p "$MESA_BASE/star" &> MESArun.out

# Log completion time
echo "Finished M=${MASS} at $(date)" >> "$OUTPUT_DIR/run_times.txt"
