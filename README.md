# septool
Integrated automatic tool for sepsis analysis.

This repository contains an automatic integrated software Septool - a snakemake pipeline for analyzing Nanopore FASTQ files for presence of microbial sources of sepsis.

## Prerequisites

1. Ubuntu 24.04 or higher (other Linux distros should work, although they were not tested)
2. Conda and Mamba installed
3. Git

## Installation
1. Clone the repository:
```
git clone https://github.com/wernerkrampl/septool.git
cd septool
```
2. Create conda environment:
```
mamba create septool-env snakemake
mamba activate septool-env
```

## Usage
1. Prepare input files:
   - Place your FASTQ files in the repository root directory
   - Ensure you have a reference genome named reference.fasta, which contains reference sequences for microbiome of interest (note that Septool will check several databases, this reference is used predominantly for validation of Septool findings)
2. Run the Septool pipeline
```
snakemake --cores <number_of_cores> --use-conda
```

## Output
- Reports and processed files will be generated in the following directories:
  - `reports/fastqc`: Quality control reports
  - `processed`: Trimmed and mapped files
  - `reports/kraken2`: Taxonomic classification reports
  - `reports/resfinder`: Resistance gene detection reports
  - `reports/amrfinder`: Antimicrobial resistance analysis reports
  - `reports/qualimap`: Mapping quality reports

