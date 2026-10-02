#!/usr/bin/env bash
set -euo pipefail

mkdir -p alignment logs

REFERENCE="reference/Acaryochloris_MBIC11017"

for SAMPLE in \
  DRR194261 \
  DRR194262 \
  DRR194263 \
  DRR194264 \
  DRR194265 \
  DRR194266
do
  echo "Aligning ${SAMPLE}..."

  bowtie2 \
    -x "${REFERENCE}" \
    -1 "raw_fastq/${SAMPLE}_1.fastq.gz" \
    -2 "raw_fastq/${SAMPLE}_2.fastq.gz" \
    2> "logs/${SAMPLE}_bowtie2.log" \
    | samtools sort \
        -@ 4 \
        -o "alignment/${SAMPLE}.sorted.bam" -

  samtools index "alignment/${SAMPLE}.sorted.bam"
  samtools quickcheck "alignment/${SAMPLE}.sorted.bam"

  echo "${SAMPLE} completed."
done

echo "Alignment and BAM processing completed."
