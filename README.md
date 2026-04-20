# Bioinformatics Data Practicum Project
The purpose of this pipeline is to provide complete analysis of IDR-conservative ATAC-Seq peak data for mouse and human samples. It takes heirarchial alignment files (.hal) and conducts peak mapping, annotation of promoters and enhancers, motif enrichment, and GO analysis.

This analysis used human and mouse ATAC-Seq data from healthy adrenal gland tissue in female subjects. Human data was from the ENCODE database (ENCSR241OBO and ENCSR864ADD) and mouse data was from Liu et al., Scientific Data, 2019. Reference genomes used throughout analysis were hg38 and mm10 respectively. 

## Installation
### Dependencies:
This tool was designed for a Linux SLURM cluster. To ensure smooth execution of complete_analysis_pipeline.sh, install the following to your cluster environment bin before running:

#### HALPER
* Python version 3.6 or 3.7 (https://www.python.org/downloads/release/python-371/)
* Python libraries `matplotlib` and `numpy`
	* numpy (http://www.numpy.org/)
		* HALPER has been tested using numpy versions 1.14.3, 1.16.0, 1.16.4, 1.16.6, and 1.18.2
	* matplotlib (https://matplotlib.org/downloads.html)
		* HALPER has been tested using matplotlib versions 1.5.1, 2.2.3, 3.2.1
* HALPER has been tested on Linux (CentOS 6, CentOS 7, and Ubuntu 18.04.4), Windows (Windows 10), and Mac

#### HOMER
- R (including ggplot2, tidyverse)
- HOMER (including the genomes hg38 and mm10)

## Each step of the pipeline can also be used individually:

### HALPER
The detailed steps for installaion can be seen in [install_hal](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)

Since the whole process is used in PSC, the part for Anaconda installation can be replaced with 
```
module load anaconda3/2024.10-1
```

### rGREAT
The detailed steps for installaion can be seen in [rGREAT](https://github.com/jokergoo/rgreat).
It is hard to install R in psc, so for isolated analysis we just use the local script in the local computer.

### HOMER
The individual scripts for running and analyzing HOMER output can be found in [HOMER](https://github.com/BioinformaticsDataPracticum2026/looking_spleentacular/blob/main/HOMER/README)

## Usage

### input
To run the full pipeline, users can copy the "RUN_FULL_PIPELINE" directory into a Linux SLURM cluster. This directory contains the script "COMPLETE_ANALYSIS_PIPELINE.sh" which can be run as such:

```
sbatch COMPLETE_ANALYSIS_PIPELINE.sh <.hal filepath> <halper_map_peak_orthologs.sh path> <bin_path>"
```

### output
COMPLETE_ANALYSIS_PIPELINE.sh will conduct peak mapping, annotation, motif enrichment, and GO analysis. The output will be organized as follows:

```
RUN_FULL_PIPELINE/
├── COMPLETE_ANALYSIS_PIPELINE.sh
├── scripts called by COMPLETE_ANALYSIS_PIPELINE.sh...
│
├── results/    									# output from mapping                              
│   ├── conservative/                            
│       ├── narrowPeak/
│       	├── shared_peaks_conservative.narrowPeak
│       	└── …
│
├── narrowPeak/										# select narrowPeak files, renamed for downstream analysis
│   ├── shared_peaks.narrowPeak
│   └── ...
│
├── homer_results/									# output from HOMER annotatePeaks.pl and findMotifsGenome.pl
│   ├── human_specific/
│   	├── annotated_peaks.txt
│   	├── annotatedPeaks.log
│   	├── findMotifs.log
│   	├── motif_instances.bed
│   	├── peaks_homer.bed
│   	└── motifs/
│   		├── nonRedundant.motifs
│   		└── ...
│   ├── mouse_specific/
│   	└── ...
│   └── shared_peaks/
│   	└── ...
│
├── filtered_annotations/										# downstream analysis from HOMER
│   ├── human_specific_motif_table.tsv
│   ├── human_specific.png
│   ├── human_specific.txt
│   └── ...
│
├── motif_annotation_split/										# additional motif analysis from HOMER
│   ├── human_enhancers_motifs.txt
│   ├── human_promoters_motifs.txt
│   └── ...
│
└── rGREAT_results/												# output from GO enrichment analysis
    ├── dot plots...
    ├── pngs...
    └── tsvs...
```

## Data Structure
### BASE/ (tools and raw peaks only)

```
BASE/
├── hal/                                    # conda env (--conda-env), e.g. conda activate …/hal
│   └── …
│
├── repos/                                  # HAL tools + HALPER checkout
│   ├── hal/bin/                            # HAL binaries (--hal-bin); on PATH (halLiftover, …)
│   └── halLiftover-postprocessing/         # HALPER (--halper-pp for PYTHONPATH)
│       ├── halper_map_peak_orthologs.sh    # default --halper-map (full path to this file)
│       └── …
│
├── MouseAtac/AdrenalGland/peak/idr_reproducibility/
│   ├── idr.conservative_peak.narrowPeak.gz
│   └── idr.conservative_peak.narrowPeak
│
└── HumanAtac/peak/idr_reproducibility/
    ├── idr.conservative_peak.narrowPeak.gz
    └── idr.conservative_peak.narrowPeak
```

#### Mapping Results
mouse_to_human_conservative.narrowPeak: Mapped Peaks from mouse into human genome
```
chr1	960807	960976	chr4:156234671-156235043:194	-1	.	-1	-1	-1	15
```
human_specific_peaks_conservative.narrowPeak: Peaks only in human data
```
chr1	180738	181589	.	590	.	2.31421	15.43961	13.78570	398

```
mouse_specific_peaks_conservative.mouse_coords.narrowPeak: Peak only in mouse data
```
chr4	156031565	156032576	.	663	.	5.01463	24.19597	21.65876	591
```
shared_peaks_conservative.narrowPeak: Peaks in both mouse and human data
```
chr1	960807	960976	chr4:156234671-156235043:194	-1	.	-1	-1	-1	15
```

### Full Repository Structure `looking_spleentacular/` (mapped peaks + analyses)

```
looking_spleentacular/
├── HALPER_pipeline/mapping/
│   ├── *.sh                            # jobs
│   └── results/conservative/           # HALPER gzip + bedtools products used here
│       ├── idr.conservative_peak.MouseToHuman.HALPER.narrowPeak.gz
│       └── narrowPeak/
│           ├── mouse_to_human_conservative.narrowPeak
│           ├── mouse_to_human_conservative.sorted.narrowPeak
│           ├── human_conservative.sorted.narrowPeak
│           ├── shared_peaks_conservative.narrowPeak
│           ├── mouse_specific_peaks_conservative.narrowPeak
│           ├── human_specific_peaks_conservative.narrowPeak
│           └── mouse_specific_peaks_conservative.mouse_coords.narrowPeak
├── HOMER/
│   ├── narrowPeak/                   
│   ├── homer_results/<category>/     # findMotifs + annotatePeaks (motifs/, logs, BEDs)
│   └── filtered_annotations/           # promoter/enhancer labels + motif tables + plots
├── HOMER_evaluation/              
├── GO_pipeline/
│   ├── my_peaks/                    
│   └── rGREAT_results/               
├── GO_plot/
├── liftOver_pipeline/
│   ├── scripts/                    
│   ├── peak_sets/
│   └── go_analysis/
├── rGREAT/
├── rGREAT_0415/
└── *.ipynb
```


