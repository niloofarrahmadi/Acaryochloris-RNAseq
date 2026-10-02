#!/usr/bin/env bash
set -euo pipefail

mkdir -p counts

featureCounts \
  -T 4 \
  -p --countReadPairs \
  -B -C \
  -s 0 \
  -t gene \
  -g ID \
  -a annotation/GCF_000018105.1_ASM1810v1_genomic.gff.gz \
  -o counts/gene_counts.txt \
  alignment/DRR194261.sorted.bam \
  alignment/DRR194262.sorted.bam \
  alignment/DRR194263.sorted.bam \
  alignment/DRR194264.sorted.bam \
  alignment/DRR194265.sorted.bam \
  alignment/DRR194266.sorted.bam

echo "featureCounts completed."
