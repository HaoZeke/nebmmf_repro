#!/usr/bin/env bash
# Outer Slurm allocation + Flux + Snakemake flux executor.
# Copied shape from TheochemUI/gpr_optim/bench/elja/scripts/sbatch_flux_in_slurm.sh.
# Morse OptBench: CI-NEB vs OCINEB (config_roneb.ini). Account chem-ui.
#
#   cd case_studies/optbench
#   export EONCLIENT=$HOME/micromamba/envs/eon320/bin/eonclient
#   FLUX_SBATCH_PARTITION=any_cpu FLUX_SBATCH_CPUS=16 FLUX_SNAKE_JOBS=16 \
#     sbatch elja_flux_in_slurm.sh
set -euo pipefail
export PATH="${HOME}/.pixi/bin:${HOME}/.local/bin:${PATH}"

PARTITION=${FLUX_SBATCH_PARTITION:-any_cpu}
ACCOUNT=${FLUX_SBATCH_ACCOUNT:-chem-ui}
CPUS=${FLUX_SBATCH_CPUS:-16}
WALL=${FLUX_SBATCH_TIME:-04:00:00}
MEM=${FLUX_SBATCH_MEM:-32G}
SNAKE_JOBS=${FLUX_SNAKE_JOBS:-16}

bench_root() {
  if [[ -n "${FLUX_BENCH_ROOT:-}" ]]; then
    cd "$FLUX_BENCH_ROOT" && pwd
  elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    cd "$SLURM_SUBMIT_DIR" && pwd
  else
    cd "$(dirname "$0")" && pwd
  fi
}

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  HERE=$(bench_root)
  mkdir -p "${HERE}/logs"
  export FLUX_BENCH_ROOT="$HERE"
  SELF=$(readlink -f "$0")
  exec sbatch \
    --job-name=ocineb-flux \
    --account="$ACCOUNT" \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$WALL" \
    --hint=nomultithread \
    --chdir="$HERE" \
    --export=ALL,FLUX_BENCH_ROOT="$HERE",EONCLIENT="${EONCLIENT:-}" \
    --output="${HERE}/logs/flux-%j.out" \
    --error="${HERE}/logs/flux-%j.err" \
    "$SELF"
fi

HERE=$(bench_root)
cd "$HERE"
if command -v pixi >/dev/null 2>&1; then
  set +u
  eval "$(pixi shell-hook -e flux 2>/dev/null || pixi shell-hook 2>/dev/null || true)"
  set -u
fi
if ! command -v flux >/dev/null 2>&1; then
  for cand in \
    "${HERE}/../../.pixi/envs/flux/bin" \
    "${HOME}/Git/Github/TheochemUI/gpr_optim/bench/elja/.pixi/envs/flux/bin"
  do
    if [[ -x "${cand}/flux" ]]; then
      export PATH="${cand}:${PATH}"
      break
    fi
  done
fi
command -v flux >/dev/null
command -v snakemake >/dev/null
if [[ -n "${EONCLIENT:-}" ]]; then
  export PATH="$(dirname "$EONCLIENT"):${PATH}"
fi

snakemake --unlock --snakefile workflow/Snakefile --configfile config/general_config.yml >/dev/null 2>&1 || true

SLOTS=${FLUX_SBATCH_CPUS:-${SLURM_CPUS_PER_TASK:-16}}
export FLUX_F58_FORCE_ASCII=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "=== flux-in-slurm start $(date -Iseconds) job=${SLURM_JOB_ID} slots=${SLOTS} ==="
flux start --test-size="$SLOTS" \
  snakemake all \
    --snakefile workflow/Snakefile \
    --configfile config/general_config.yml \
    --profile profile/flux \
    --jobs "$SNAKE_JOBS" \
    --rerun-incomplete \
    --keep-going
echo "=== flux-in-slurm end $(date -Iseconds) rc=$? ==="
