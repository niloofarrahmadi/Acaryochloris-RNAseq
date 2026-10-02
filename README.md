# Acaryochloris marina RNA-seq Analysis

This repository contains the computational workflow and results for a focused RNA-seq analysis of the short-term response of *Acaryochloris marina* MBIC11017 to orange and far-red light.

## Research question

Does a one-hour shift to orange light produce a coordinated transcriptional response in phycobilisome-associated genes relative to a one-hour shift to far-red light?

## Biological background

*Acaryochloris marina* MBIC11017 is a photosynthetic cyanobacterium containing chlorophyll d and phycobilisomes. This project focuses on whether a short change in light quality is associated with increased transcript levels of genes encoding structural components of the phycobilisome.

## Data source

The RNA-seq data were obtained from the publicly available dataset of Kashimoto et al. (2020).

The experiment included six paired-end RNA-seq samples:

| Run | Condition | Replicate |
|---|---|---:|
| DRR194261 | FarRed | 1 |
| DRR194262 | Orange | 1 |
| DRR194263 | FarRed | 2 |
| DRR194264 | Orange | 2 |
| DRR194265 | FarRed | 3 |
| DRR194266 | Orange | 3 |

Cells were shifted from dim white light to either orange or far-red light for one hour before RNA collection. The two conditions were tested at the same light intensity.

Paper:
Kashimoto et al. (2020), *Journal of General and Applied Microbiology*

DOI: 10.2323/jgam.2019.11.008

The raw sequencing data are available through the NCBI Sequence Read Archive (SRA) under the six run accessions listed above.

## Analysis workflow

The analysis was performed in two environments:

- Linux/WSL for sequence quality control, read alignment, and gene-level counting
- RStudio with R 4.3.3 for statistical analysis and visualization

The main workflow was:

```text
Raw paired-end FASTQ
        ↓
FastQC quality control
        ↓
Bowtie2 genome alignment
        ↓
Samtools sorting and indexing
        ↓
featureCounts gene-level counting
        ↓
DESeq2 differential expression analysis
        ↓
Focused analysis of 15 structural phycobilisome-associated genes
        ↓
PCA and heatmap visualization
        ↓
Sensitivity analysis

Based on the quality-control results, the raw reads were retained without adapter trimming.

Reference genome and annotation

Reference assembly:

GCF_000018105.1 (ASM1810v1)

Organism:

Acaryochloris marina MBIC11017

Gene-level counting was performed using the corresponding NCBI RefSeq GFF annotation.

Differential expression analysis

Gene-level raw counts were analyzed using DESeq2.

Design:

~ condition

Main contrast:

Orange vs FarRed

FarRed was used as the reference condition.

Low-count filtering in the main analysis retained genes with at least 10 counts in at least 3 samples.

For each gene, DESeq2 provided a log2 fold change and an adjusted p-value.

The main thresholds were:

adjusted p-value < 0.05
|log2 fold change| > 1
Focused phycobilisome analysis

Rather than focusing only on genes with the largest differential-expression values, the analysis examined a targeted set of 15 structural phycobilisome-associated genes identified from the genome annotation.

The set included annotated phycocyanin subunits and linker proteins.

Main results

For the 15 selected structural phycobilisome-associated genes:

15/15 had positive log2 fold changes under orange light
8/15 had log2FC > 1
5/15 had adjusted p-value < 0.05
5/15 met both thresholds
Mean log2 fold change = 2.02
Median log2 fold change = 1.21

The strongest responses were observed for several phycocyanin subunit genes, with log2 fold changes of approximately 4.

The PCA showed separation of the orange- and far-red-light samples mainly along PC1. PC1 and PC2 explained 57.3% and 23.4% of the variance, respectively.

Sensitivity analysis

To test whether the main result depended on the low-count filtering choice, the DESeq2 analysis was repeated using a stricter filter requiring at least 10 counts in at least 4 samples.

The overall pattern of the 15 selected genes did not change:

15/15 remained positive
8/15 still had log2FC > 1
5/15 still had adjusted p-value < 0.05
5/15 still met both thresholds

This indicates that the main conclusion was not strongly affected by this filtering choice.

Repository structure
Acaryochloris-RNAseq/
├── .gitignore
├── README.md
├── metadata/
│   └── samples.csv
├── counts/
│   ├── gene_counts.txt
│   └── gene_counts.txt.summary
├── scripts/
│   ├── 01_fastqc.sh
│   ├── 02_alignment.sh
│   ├── 03_featureCounts.sh
│   ├── 04_deseq2.R
│   └── 05_phycobilisome_analysis.R
└── results/
    ├── PCA.png
    ├── phycobilisome_heatmap.png
    ├── deseq2_results.csv
    ├── normalized_counts.csv
    ├── pca_data.csv
    ├── phycobilisome_structural_genes.csv
    ├── phycobilisome_summary.csv
    ├── sensitivity_15_genes.csv
    ├── sensitivity_summary.csv
    ├── size_factors.csv
    └── software_versions.txt

Large raw FASTQ files, BAM files, reference sequences, annotation files, and FastQC output files are not stored in this repository. The public data accessions and reference information are documented above.

Software
FastQC 0.12.1
Bowtie2 2.5.4
Samtools 1.21
featureCounts 2.0.8
R 4.3.3
RStudio
DESeq2
pheatmap
Reproducibility

The analysis scripts are provided in the scripts/ directory.

The workflow can be followed in order:

01_fastqc.sh — quality control
02_alignment.sh — Bowtie2 alignment and BAM processing
03_featureCounts.sh — gene-level read counting
04_deseq2.R — DESeq2 analysis and PCA data
05_phycobilisome_analysis.R — focused gene-set analysis, heatmap, and sensitivity analysis

The main input count matrix and sample metadata are provided in counts/ and metadata/, while analysis outputs are provided in results/.
