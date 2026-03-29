#!/bin/bash
#SBATCH -A bio230007p        # your project/account
#SBATCH -p RM-shared          # or general / RM-shared partition
#SBATCH --time=01:00:00       # max runtime
#SBATCH -n 1               # 1 CPU
#SBATCH --mem=2000         # 2000 MB total (~2 GB); raise if the step OOMs
#SBATCH -o generate_narrowPeak_%j.out      # stdout log
#SBATCH -e generate_narrowPeak_%j.err      # stderr log
#SBATCH -J generate_narrowPeak # job name



# Run from the directory that contains the liftover outputs and human.bed
cd /ocean/projects/bio230007p/wli27/output/Mouse/mapping/new

# psc module load bedtools
module load bedtools
# Unzip the target file
gunzip -c idr.optimal_peak.MouseToHuman.HALPER.narrowPeak.gz > mouse_to_human.narrowPeak

# can sort for quick use, optional
sort -k1,1 -k2,2n mouse_to_human.narrowPeak > mouse_to_human.sorted.narrowPeak

sort -k1,1 -k2,2n /ocean/projects/bio230007p/wli27/HumanAtac/peak/idr_reproducibility/idr.optimal_peak.narrowPeak > human.sorted.narrowPeak

# shared mapped
bedtools intersect -a mouse_to_human.sorted.narrowPeak -b human.sorted.narrowPeak -u > shared_peaks.narrowPeak

# Mouse mapped not open in human peaks (not overlapping human)
bedtools intersect -a mouse_to_human.sorted.narrowPeak -b human.sorted.narrowPeak -v > mouse_specific_peaks.narrowPeak

# Human mapped not open in mouse peaks (not overlapping mouse)
bedtools intersect -a human.sorted.narrowPeak -b mouse_to_human.sorted.narrowPeak -v > human_specific_peaks.narrowPeak
