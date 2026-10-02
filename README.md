# Acaryochloris marina RNA-seq Analysis

RNA-seq analysis of the short-term transcriptional response of *Acaryochloris marina* MBIC11017 to orange and far-red light.

## Research question

Does a one-hour shift to orange light produce a coordinated transcriptional response in phycobilisome-associated genes relative to a one-hour shift to far-red light?

## Data source

RNA-seq data from Kashimoto et al. (2020).

SRA runs:

- DRR194261
- DRR194262
- DRR194263
- DRR194264
- DRR194265
- DRR194266

The experiment included three biological replicates under orange light and three under far-red light.

Paper:
Kashimoto et al. (2020), Journal of General and Applied Microbiology.
DOI: 10.2323/jgam.2019.11.008

## Analysis workflow

Raw paired-end FASTQ files were analyzed using the following workflow:

FASTQ quality control  
→ genome alignment  
→ gene-level read counting  
→ DESeq2 differential expression analysis  
→ focused analysis of 15 structural phycobilisome-associated genes

Reads were retained without adapter trimming based on the quality-control results.

## Reference genome

Organism: *Acaryochloris marina* MBIC11017

NCBI assembly:
GCF_000018105.1 (ASM1810v1)

Gene annotation:
NCBI RefSeq GFF annotation associated with GCF_000018105.1

## Main methods

- FastQC 0.12.1
- Bowtie2 2.5.4
- Samtools 1.21
- featureCounts 2.0.8
- R 4.3.3 in RStudio
- DESeq2

DESeq2 design:

`~ condition`

Main contrast:

`Orange vs FarRed`

Low-count filtering in the main analysis:
at least 10 counts in at least 3 samples.

A sensitivity analysis used a stricter filter requiring at least 10 counts in at least 4 samples.

## Main focused result

The focused analysis included 15 structural phycobilisome-associated genes identified from the genome annotation.

In the main analysis:

- 15/15 genes had positive log2 fold changes
- 8/15 had log2FC > 1
- 5/15 had adjusted p-value < 0.05
- 5/15 met both criteria

The same overall pattern was observed after applying the stricter low-count filter.

## Repository contents

`counts/`  
Gene-level count matrix and featureCounts summary.

`metadata/`  
Sample condition information.

`scripts/`  
Analysis scripts used for the RNA-seq workflow.

`results/`  
Main figures and focused analysis results.

Large raw FASTQ, BAM, reference genome, annotation, and quality-control files are not stored in this repository. Their public sources and accession information are documented above.
