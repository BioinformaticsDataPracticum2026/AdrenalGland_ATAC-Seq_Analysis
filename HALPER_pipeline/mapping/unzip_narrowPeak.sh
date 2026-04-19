#!/bin/bash
#SBATCH -A bio230007p        # your project/account
#SBATCH -p RM-shared          # or general / RM-shared partition
#SBATCH --time=06:00:00       # max runtime
#SBATCH -n 1               # 1 CPU
#SBATCH --mem=2000         # 2000 MB total (~2 GB)
#SBATCH -o unzip_narrowpeak_%j.out      # stdout log
#SBATCH -e unzip_narrowpeak_%j.err      # stderr log
#SBATCH -J unzip_narrowpeak # job name

set -euo pipefail

BASE=""
# parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base) BASE="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 --base DIR"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: $0 --base DIR" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${BASE}" ]]; then
  echo "Error: --base DIR is required." >&2
  echo "Usage: $0 --base DIR" >&2
  exit 1
fi
BASE="${BASE%/}"
# load the conda environment
module load anaconda3/2024.10-1
conda activate "${BASE}/hal"
export PATH="${BASE}/repos/hal/bin:${PATH}"
export PYTHONPATH="${BASE}/repos/halLiftover-postprocessing:${PYTHONPATH:-}"
# unzip the narrowPeak file
cd "${BASE}/MouseAtac/AdrenalGland/peak/idr_reproducibility"
zcat idr.conservative_peak.narrowPeak.gz > idr.conservative_peak.narrowPeak
# unzip the narrowPeak file
cd "${BASE}/HumanAtac/peak/idr_reproducibility"
zcat idr.conservative_peak.narrowPeak.gz > idr.conservative_peak.narrowPeak
