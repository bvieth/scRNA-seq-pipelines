
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scRNA-seq-pipelines

This repository is a compendium to “A Systematic Evaluation of Single
Cell RNA-Seq Analysis Pipelines” by Vieth et al. 2019, currently on
[bioRxiv](https://www.biorxiv.org/content/10.1101/583013v1).

All code in the repository is distributed under the GPL-3 license.

For examples on how to use the power analyis tool powsimR, please see
[powsimR](https://github.com/bvieth/powsimR).

For any questions or issues with the code in this repository, please use
the “Issues” tab.

# Getting started

Below you will find a brief outline of the analysis. Please refer to the
corresponding folders for detailed information.

## Dependencies

You will need the following software:

  - samtools version 1.9
  - BWA version 0.7.12
  - kallisto 0.43.1
  - zUMIs 2.4.5
  - powsimR (see prepackaged version in simulation folder)

## Data aquisition

### scRNA-seq data sets

The main input for the simulations were the scRNA-seq experiments from
[Ziegenhain et
al., 2017](https://www.sciencedirect.com/science/article/pii/S1097276517300497?via%3Dihub).
The fastq files can be downloaded from
[GSE75790](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75790).
The cDNA reads were cut to 45 bp length and reads were randomly
downsampled to 1 million reads per cell to be comparable. In addition,
the fastq files of 1k 1:1 mixture of HEK293T and NIH3T cells from [10X
Genomics
support](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/1k_hgmm_v3)
were downloaded and reads from the QC-passed NIH3T cells were extracted.
For the comparison of pipelines, we processed the expression profiles of
~ 1000 human
[PBMCs](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3)
from 10X Genomics.

## Annotation

The annotation folder contains detailed information on how to get the
necessary files from Gencode, Vega and RefSeq and processing them in
order to use them for alignments. Additionally, we provide ERCC spike-in
sequences.

## Alignment

The alignment folder contains the commands per aligner for indexing,
mapping, counting expression per gene and subsequently summarising the
expression into one gene count matrix per protocol, aligner and
annotation combination. In particular, it contains R scripts for
filtering out multi-mapped reads for bwa and kallisto alignments.
Additionally, we provide an RNA-seq experiment output of one cell (fastq
file) as exemplary input for bwa and kallisto.

## Simulations

We use powsimR for our simulations. The folder simulation contains a
walkthrough R script to estimate the distributional parameters necessary
to simulate gene expression, to set up differential expression and to
finally run the simulations. Additionally, we provide exemplary input of
SCRB-seq (UMI) and Smart-seq2 data to run simulations.
