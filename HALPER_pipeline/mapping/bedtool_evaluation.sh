#!/bin/bash
#SBATCH -A bio230007p        #  project/account
#SBATCH -p RM-shared          # RM-shared partition
#SBATCH --time=01:00:00       # max runtime
#SBATCH -n 1               # 1 CPU
#SBATCH --mem=2000         # 2000 MB total 
#SBATCH -o generate_narrowPeak_%j.out      # stdout log
#SBATCH -e generate_narrowPeak_%j.err      # stderr log
#SBATCH -J generate_narrowPeak # job name

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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HALPER_OUT="${SCRIPT_DIR}/results/conservative/narrowPeak"
HUMAN_CONSERVATIVE_NARROWPEAK="${BASE}/HumanAtac/peak/idr_reproducibility/idr.conservative_peak.narrowPeak"

mkdir -p "${HALPER_OUT}"
cd "${HALPER_OUT}"

# psc module load bedtools
module load bedtools
# Unzip the target file with HALPER filter
gunzip -c idr.conservative_peak.MouseToHuman.HALPER.narrowPeak.gz > mouse_to_human_conservative.narrowPeak

# can sort for quick use, optional
sort -k1,1 -k2,2n mouse_to_human_conservative.narrowPeak > mouse_to_human_conservative.sorted.narrowPeak

sort -k1,1 -k2,2n "${HUMAN_CONSERVATIVE_NARROWPEAK}" > human_conservative.sorted.narrowPeak

# shared mapped
bedtools intersect -a mouse_to_human_conservative.sorted.narrowPeak -b human_conservative.sorted.narrowPeak -u > shared_peaks_conservative.narrowPeak

# Mouse specific (actually mouse_nonreciprocal_peaks)
bedtools intersect -a mouse_to_human_conservative.sorted.narrowPeak -b human_conservative.sorted.narrowPeak -v > mouse_specific_peaks_conservative.narrowPeak

# Human specific
bedtools intersect -a human_conservative.sorted.narrowPeak -b mouse_to_human_conservative.sorted.narrowPeak -v > human_specific_peaks_conservative.narrowPeak
