# Fertility Screen

This repository contains the analysis of BARseq data aimed at studying fertility phenotypes in a strain of Plasmodium berghei, a rodent malaria parasite.

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
## Concepts and Descriptions
In this section, we will explain the terminology and concepts that were employed in the calculation of male and female fertility.

### Relative abundance
To calculate the abundance of mutants in each sample, we averaged the forward and reverse reads. Subsequently, we computed the abundance for each pool. For example, for a pool of 3 mutants/genes, the count matrix can be generated as shown in the following table.

| Genes | Sample  1 | Sample  2 |
|----------|----------|----------|
| Gene1 | 10.5 | 11.5 |
| Gene2 | 30 | 10.5 |
| Gene3|  20.5| 5.5 |



We determine the relative abundance by dividing the abundance of a specific feature by the total abundance in the sample. The resulting relative abundance table is as
follows:

| Genes | Sample  1 | Sample  2 |
|----------|----------|----------|
| Gene1 | 10.5/60 | 11.5/31 |
| Gene2 | 30/60 | 10/31 |
| Gene3|  20.5/60| 5.5/31 |

Throughout the screens, we assessed biological variance by transmitting each pool at least twice from different infected mice. We further monitored sampling errors by subsampling each transmission experiment in two batches of mosquitoes and by amplifying and counting barcodes twice from each DNA sample. Data from all types of replicates were entered into the final analysis and used to determine independently for each sex the confidence with which a mutant’s transmission rate was determined.
![Sample description](https://github.com/vpandey-om/Fertility_screen/blob/master/output/sample.png)
## Usage
### Convert Fastq to count matrix
To convert paired forward and reverse reads from BARseq experiments into a count matrix, we employ the following command.
~~~
snakemake --use-conda --conda-frontend conda raw_barseq  --cores 1 -F --config input_file=barseq-raw/testdata/sequence barcode_file=barcode_to_gene_210920.csv output_dir=barseq-raw/testRes -p  
~~~
The process requires two key inputs: a directory (e.g., `barseq-raw/testdata/sequence`) where all the fastq files for the samples are stored and a CSV file (e.g., `barcode_to_gene_210920.csv`) containing barcode information for each gene or mutant. The resulting output is directed to another directory (e.g., `barseq-raw/testRes`), where both the mutant count matrix file and a log file are generated.

To identify and remove mutants with zero counts in all samples, use the following command. This will generate the file `removed_zeros_barcode_counts_table.csv`
in the `barseq-raw/testRes` directory.
~~~
snakemake --use-conda --conda-frontend conda remove_zeros  --cores 1 -F --config output_dir=barseq-raw/testRes -p
~~~




### Combine fertility screen data

In the Barseq experiment, we analyzed over 1200 mutants distributed across seven pools labeled as Pool1 to Pool7. Notably, Pool5 and Pool7 were studied twice. To consolidate the data from all these pools (Pool1 to Pool7), you can use the following command.
~~~
snakemake --use-conda --conda-frontend conda combine_pools  --cores 1 -F
~~~
This command/script will calculate male and female fertility, their variances, and male and female phenotypes (e.g., Reduced, Not reduced).

### Motility screen data
We conducted a BARseq experiment, focusing on the pool with the most strongest male fertility phenotypes. Subsequently, we collected barcodes from purified microgametes as part of a motility screen.

~~~
snakemake --use-conda --conda-frontend conda motility_screen  --cores 1 -F
~~~
This command/script will calculate motility rate, their variances, and motility phenotype (e.g., Reduced, Not reduced).

### Generate Visualizations
To visualize ranked male fertility and female fertility, as well as create a scatter plot of male vs. female fertility, you can utilize the following command.
~~~
snakemake --use-conda --conda-frontend conda plot_fertility  --cores 1 -F
~~~

We performed an enrichment analysis using the Malaria Parasite Metabolic Pathways (MPMP) for male and female-only fertility phenotypes. For this use following command.

~~~
snakemake --use-conda --conda-frontend conda mpmp_enrichment  --cores 1 -F
~~~

We visualize the enrichment of top-ranked pathways (male/female) using violin plots with the following commands.

~~~
snakemake --use-conda --conda-frontend conda mpmp_violin  --cores 1 -F
~~~

To investigate the relationship between fertility phenotypes and gene expression, we visualized genes associated with male and female-only fertility phenotypes within the Malaria Cell Atlas (mca) gene expression-based clusters.

~~~
snakemake --use-conda --conda-frontend conda mca_gene_plot  --cores 1 -F
~~~

To filter out noisy data, we generated a plot of relative abundance against relative error. Our analysis revealed that data with a relative abundance below the cutoff value of log2(-12) exhibited unacceptably high errors, making it unsuitable for further analysis.

~~~
snakemake --use-conda --conda-frontend conda error_plot  --cores 1 -F
~~~
