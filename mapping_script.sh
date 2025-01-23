#!/bin/bash

BOWTIE2_INDEX=~/human_reference_genome_index_for_bowtie2_name-hg38/hg38

echo "index file for hg38 genome in ~/human_reference_genome_index_for_bowtie2_name-hg38/hg38"

echo "write the path to your fastq file 1"

read FASTQ_FILE_1

echo "write the path to your fastq file 2"

read FASTQ_FILE_2

echo "write number of threads for bowtie2"

read THREADS

echo "write output_dir for bam_sorted_file"

read OUTPUT_DIR

echo "write the name for output SORTED BAM file (use .bam at the end)"

read OUTPUT_FILE_NAME

bowtie2 --threads $THREADS -x $BOWTIE2_INDEX --very-sensitive-local -1 $FASTQ_FILE_1 \
-2 $FASTQ_FILE_2 | \
samtools view -@ $THREADS -S -b | \
samtools sort -@ $THREADS > $OUTPUT_DIR/$OUTPUT_FILE_NAME
