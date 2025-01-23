#!/bin/bash
#$1 - FASTQ1
#$2 - FASTQ2

for file in ./*.fastq; do
# Remove 10 bp from the start and the end of fastq file 1
	cutadapt -j 128 -u 10 -u -10 -o "trimmed_10_bp.fastq${file#??}" $file
done
# Remove 10 bp from the start and the end of fastq file 2
#cutadapt  -u 10 -u -10 -o #minus_10_bp_${2}" $2
