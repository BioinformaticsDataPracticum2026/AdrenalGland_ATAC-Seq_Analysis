#!/bin/bash
#SBATCH -A bio230007p        # your project/account
#SBATCH -p RM-shared          # or general / RM-shared partition
#SBATCH --time=06:00:00       # max runtime
#SBATCH -n 2               # 2 CPUs (halLiftover is mostly single-threaded; -n 1 is often enough)
#SBATCH --mem=2000         # 2000 MB total (~2 GB); raise if the step OOMs
#SBATCH -o halper_%j.out      # stdout log
#SBATCH -e halper_%j.err      # stderr log
#SBATCH -J halper_mouse2human_adrenal # job name

module load anaconda3/2024.10-1
conda activate /ocean/projects/bio230007p/wli27/hal
export PATH=/ocean/projects/bio230007p/wli27/repos/hal/bin:${PATH}
export PYTHONPATH=/ocean/projects/bio230007p/wli27/repos/halLiftover-postprocessing:${PYTHONPATH}

bash /ocean/projects/bio230007p/wli27/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/wli27/MouseAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak \
  -o /ocean/projects/bio230007p/wli27/output/Mouse/mapping/new/ \
  -s Mouse \
  -t Human \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
