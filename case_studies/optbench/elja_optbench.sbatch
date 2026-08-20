#!/bin/bash
# Elja CPU array: Morse Pt-island cineb vs mmf (OCINEB; was RONEB).
# Shipped INI files and endpoints. Account/partition from existing chem-ui jobs.
#SBATCH --job-name=ocineb-optbench
#SBATCH --account=chem-ui
#SBATCH --partition=any_cpu
#SBATCH --time=02:00:00
#SBATCH --hint=nomultithread
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=0-58%20
#SBATCH --output=logs/optbench-%A_%a.out
#SBATCH --error=logs/optbench-%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"

EONCLIENT="${EONCLIENT:-}"
if [[ -z "$EONCLIENT" ]]; then
  if command -v eonclient >/dev/null 2>&1; then
    EONCLIENT=$(command -v eonclient)
  elif [[ -x "$HOME/.pixi/envs/eon/bin/eonclient" ]]; then
    EONCLIENT="$HOME/.pixi/envs/eon/bin/eonclient"
  elif [[ -x "$ROOT/../../.pixi/envs/eon/bin/eonclient" ]]; then
    EONCLIENT="$ROOT/../../.pixi/envs/eon/bin/eonclient"
  else
    echo "eonclient not on PATH; set EONCLIENT" >&2
    exit 1
  fi
fi

# 59 products product_0.con .. product_58.con; both methods in one task
# (chem-ui MaxSubmit=100; a 118-element array is rejected).
PROD=${SLURM_ARRAY_TASK_ID}
SYS="pt_prod_${PROD}"
for spec in cineb:resources/config_nebci.ini mmf:resources/config_roneb.ini; do
  METHOD="${spec%%:*}"
  INI="${spec##*:}"
  OPATH="results/03_neb/${SYS}/${METHOD}"
  mkdir -p "$OPATH"
  cp "$INI" "$OPATH/config.ini"
  cp "resources/pt-island-con/reactant.con" "$OPATH/reactant.con"
  cp "resources/pt-island-con/product_${PROD}.con" "$OPATH/product.con"
  ( cd "$OPATH" && "$EONCLIENT" )
  echo "DONE ${SYS} ${METHOD} host=$(hostname) eon=${EONCLIENT}"
done
