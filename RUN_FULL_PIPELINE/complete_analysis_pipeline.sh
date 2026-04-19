#!/bin/bash

#SBATCH --job-name=full_ATAC_peak_analysis
#SBATCH --output=ATAC_%j.log
#SBATCH --error=ATAC_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        
#SBATCH --mem=8000M                
#SBATCH --time=15:00:00          
#SBATCH --account=bio230007p

set -euo pipefail

SCRIPT_DIR="$SLURM_SUBMIT_DIR"
WD="${3:-}"

# NEED TO UPDATE: make sure we have all path convensions are standardized, finalize SLURM job requirements
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <narrowpeak_dir> <bin_path> [wd]"
    echo "Example: $0 ./narrowPeak /ocean/projects/bio230007p/mccreary/bin ./"
    exit 1
fi


# NEED TO COPY ALL RELEVANT SCRIPTS INTO DIRECTORY
echo "Starting mapping with integrate_halper.sh..."
bash "$SCRIPT_DIR/integrate_halper.sh" # inputs...
echo "Starting motif enrichment and peak annotation with run_homer_on_dir.sh..."
bash "$SCRIPT_DIR/run_full_annotatePeaks_findMotifs.sh" # inputs...
echo "Starting GO BP enrichment with rGREAT run_GO_pipeline.sh..."
bash "$SCRIPT_DIR/run_GO_pipeline.sh" # inputs...
echo "Done."