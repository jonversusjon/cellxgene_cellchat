#!/bin/bash
#SBATCH -J cellchat_disease_subsets
#SBATCH -p cpu_high_mem
#SBATCH --mem=600G
#SBATCH --cpus-per-task=13
#SBATCH --time=24:00:00

set -euo pipefail

# -------------------- config validation --------------------
if [ -z "${1-}" ]; then
  echo "[ERROR] Config file required. Usage: sbatch $0 <config.yaml>" >&2
  exit 1
fi
CFG="$1"
if [ ! -f "$CFG" ]; then
  echo "[ERROR] Config file not found: $CFG" >&2
  exit 1
fi

# -------------------- derive run script path --------------------
# Prefer scripts/run_cellchat_batch.R, else run_cellchat_batch.R
RUN_R="scripts/run_cellchat_batch.R"
[ -f "$RUN_R" ] || RUN_R="run_cellchat_batch.R"
if [ ! -f "$RUN_R" ]; then
  echo "[ERROR] Cannot locate run_cellchat_batch.R (looked in scripts/ and CWD)" >&2
  exit 1
fi

# -------------------- extract outdir from YAML --------------------
# Looks for any 'outdir:' key (e.g., under run:), keeps last occurrence
OUTDIR_BASE=$(grep -E '^[[:space:]]*outdir:' "$CFG" | tail -n1 | sed -E 's/.*outdir:[[:space:]]*"?([^"]+)"?.*/\1/' | tr -d '"')
if [ -z "${OUTDIR_BASE}" ]; then
  echo "[ERROR] Could not find 'outdir:' in config file" >&2
  exit 1
fi

# Append SLURM job ID
OUTDIR="${OUTDIR_BASE}_${SLURM_JOB_ID:-manual}"

# Create directory structure
mkdir -p "${OUTDIR}/slurm_logs" "${OUTDIR}/logs"

# Env overrides used by the R code
export CELLCHAT_OUTDIR="$OUTDIR"
export HDF5_USE_FILE_LOCKING=FALSE  # safer on NFS/HPC

# -------------------- logging --------------------
TS=$(date +%Y%m%d_%H%M%S)
exec > >(tee -a "${OUTDIR}/slurm_logs/stdout_${SLURM_JOB_ID:-noid}_${TS}.log") \
     2> >(tee -a "${OUTDIR}/slurm_logs/stderr_${SLURM_JOB_ID:-noid}_${TS}.log" >&2)

echo "[INFO] Host: $(hostname)"
echo "[INFO] PWD : $(pwd)"
echo "[INFO] Job : ${SLURM_JOB_ID:-interactive}"
echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] SLURM logs: ${OUTDIR}/slurm_logs/"
echo "[INFO] R logs: ${OUTDIR}/logs/"
echo "[INFO] SLURM allocated CPUs: ${SLURM_CPUS_PER_TASK:-1}"
echo "[INFO] SLURM allocated Mem  : ${SLURM_MEM_PER_NODE:-unknown} MB"

# -------------------- avoid over-threading in math libs --------------------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# -------------------- conda activation inside Slurm --------------------
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
  source "/opt/conda/etc/profile.d/conda.sh"
else
  echo "[ERROR] conda.sh not found. Load your conda module or fix the path." >&2
  exit 1
fi

conda activate r-4.4-cellchatv2
echo "[INFO] Conda env: $(conda info --envs | awk '/\*/{print $1, $NF}')"
R --version | head -n1 || true

# -------------------- resolve cores safely --------------------
# Prefer SLURM_CPUS_PER_TASK; fallback to SLURM_TRES_PER_TASK (parse cpu=K); then scontrol; then 1
CORES="${SLURM_CPUS_PER_TASK:-}"
if [ -z "${CORES}" ] && [ -n "${SLURM_TRES_PER_TASK:-}" ]; then
  CORES="$(echo "${SLURM_TRES_PER_TASK}" | sed -n 's/.*cpu=\([0-9]\+\).*/\1/p')"
fi
if [ -z "${CORES}" ] && [ -n "${SLURM_JOB_ID:-}" ]; then
  CORES="$(scontrol show job "${SLURM_JOB_ID}" | tr ' ' '\n' | sed -n 's/^CpusPerTask=\([0-9]\+\)$/\1/p' | head -n1)"
fi
: "${CORES:=1}"
echo "[INFO] Resolved cores: ${CORES}"

# -------------------- avoid srun env-var conflict --------------------
# (Your cluster complains when both SLURM_CPUS_PER_TASK and SLURM_TRES_PER_TASK are set.)
unset SLURM_CPUS_PER_TASK
unset SLURM_TRES_PER_TASK

# -------------------- config echo (reflects new pipeline) --------------------
echo "------------------------------------------------------------"
echo "[INFO] Using config: $CFG"
echo "[INFO] Mode: SingleCellExperiment + zellkonverter (HDF5-backed, no Seurat)"
echo "[INFO] Key YAML fields (for quick sanity):"
# run + input
grep -nE '^[[:space:]]*(run:|input:|metadata:|database:|cell_communication:|comm_prob:)' -n "$CFG" || true
# details we care about now
grep -nE '^[[:space:]]*(outdir:|species:|files:|dir:|pattern:|cell_type_col:|donor_col:|disease_col:|control_value:|disease_level:)' "$CFG" || true
grep -nE '^[[:space:]]*disease_comparison:' -n "$CFG" || true
grep -nE '^[[:space:]]*(enable:|min_cells_per_condition:|min_cells_per_analysis:)' "$CFG" || true
grep -nE '^[[:space:]]*analyses:' -n "$CFG" || true
echo "------------------------------------------------------------"

# Friendly note about subsets actually used
if grep -qE '^[[:space:]]*disease_comparison:[[:space:]]*$' "$CFG" && \
   awk 'f&&/enable:/ {print; exit} /^[[:space:]]*disease_comparison:/ {f=1}' "$CFG" | grep -q 'enable:[[:space:]]*true'; then
  echo "[INFO] This run will produce subset analyses if possible:"
  echo "       - full dataset"
  echo "       - normal-only (by donor_context)"
  echo "       - disease-only (by donor_context)"
else
  echo "[INFO] disease_comparison.enable is not TRUE; running full dataset only."
fi

# -------------------- run --------------------
# Pass allocated CPUs as CLI override so the R runner prioritizes it (CLI > config > SLURM > default)
echo "[INFO] Launching: Rscript ${RUN_R} ${CFG} ${CORES}"
Rscript "${RUN_R}" "${CFG}" "${CORES}"
EXIT=$?

echo "[INFO] Finished with exit code: $EXIT"
exit $EXIT