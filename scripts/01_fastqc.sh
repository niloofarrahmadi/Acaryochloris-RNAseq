#!/usr/bin/env bash
set -euo pipefail

mkdir -p qc_before

fastqc \
  -t 4 \
  -o qc_before \
  raw_fastq/*.fastq.gz

echo "FastQC completed."
