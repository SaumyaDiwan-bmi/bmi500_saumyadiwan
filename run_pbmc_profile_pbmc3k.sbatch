#!/bin/bash

#SBATCH --job-name=pbmc_profile
#SBATCH --account=general
#SBATCH --partition=c64-m512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-user=saumya.diwan@emory.edu

set -euo pipefail

# ---------------- USER SETTINGS ----------------
PROJECT_DIR=/users/sdiwan2/code/bmi500_project/bmi500
DATA_DIR=/users/sdiwan2/code/bmi500_project/bmi500/data
DATASET=pbmc3k
OUT_DIR=${PROJECT_DIR}/outputs/${DATASET}
PYTHON_SCRIPT=${PROJECT_DIR}/scanpy_pbmc_profiled.py
# ------------------------------------------------

cd "$PROJECT_DIR"

mkdir -p logs "$OUT_DIR"

source /opt/anaconda/etc/profile.d/conda.sh
conda activate bmi500

THREADS=${SLURM_CPUS_PER_TASK}

export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"

echo "============================================================"
echo "Dataset:      $DATASET"
echo "Data dir:     $DATA_DIR"
echo "Output dir:   $OUT_DIR"
echo "Python file:  $PYTHON_SCRIPT"
echo "CPUs:         $THREADS"
echo "Node:         $(hostname)"
echo "Start time:   $(date)"
echo "============================================================"

COARSE_FILE="${OUT_DIR}/${DATASET}.time_v.txt"
CONSOLE_FILE="${OUT_DIR}/${DATASET}.console.txt"

/usr/bin/time -v -o "$COARSE_FILE" \
    python "$PYTHON_SCRIPT" \
        --data-dir "$DATA_DIR" \
        --data-set "$DATASET" \
        --out-dir "$OUT_DIR" \
        --num-threads "$THREADS" \
    2>&1 | tee "$CONSOLE_FILE"

echo
echo "================ COARSE-GRAIN SUMMARY ====================="
grep -E \
"User time \(seconds\)|System time \(seconds\)|Percent of CPU this job got|Elapsed \(wall clock\) time|Maximum resident set size \(kbytes\)" \
"$COARSE_FILE" || true

echo
echo "================ INSTRUMENTED TIMING ======================"
grep "^\[PROFILE\]" "$CONSOLE_FILE" || true

echo
echo "Finished: $(date)"
echo "============================================================"
