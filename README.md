# Fertility_screen

This repository contains the analysis of bar-seq data aimed at studying fertility phenotypes in a strain of Plasmodium berghei, a rodent malaria parasite.

## Overview
We used barcoded PlasmoGEM vectors to mutagenize P. berghei strains that produce exclusively fertile male or female gametocytes. Our comprehensive screening covered over 1200 targetable genes, allowing us to probe sex-specific phenotypes. The outcome of our study revealed the identification of numerous genes with specific functions in sexual reproduction.

Within this repository, you will find scripts and datasets detailing our analysis of raw bar-seq data and the computation of fertility phenotype metrics. Additionally, we have included the scripts used to generate figures featured in our research paper.


## Prerequisites

Users need to install before using the Snakemake workflow.

- Python (>=3.7)
- Snakemake (7.32.4)

## Installation

Install Snakemake using pip.
  ~~~
  pip install snakemake
  ~~~

## Usage

Provide instructions for running the Snakemake workflow. Include sample commands and arguments, if applicable.


## Descriptions
Barseq experiment was done in 7 pools from pool1 to pool7. We have pool5 and pool7 was done two times. For combining all pools from pool1 to pool7 one could use script form codes ‘combine_all_pool.py’.

Scripts from the location ‘analysisCodes/scripts’ are used to plot figures in papers.


1)	For plotting ordered relative fertility error one could use script ‘plot_s_curve.py’ which is in folder: analysisCodes/scripts.
2)	MPMP pathway enrichment analysis was done by using the script: mpmp_enrichment.py form the folder: ‘analysisCodes/scripts’  
3)	 Comparison of fertility screen data with different screens we used (violin plots are generated) scripts: ‘violin_and_proteomics.py’ and ‘phospho_proteome.py’
4)	 combined error analysis to see noise in data we used script ‘error_combine_analysis.py’
