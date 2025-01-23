#!/bin/bash

echo 'write the path to bam file'

read BAM_FILE_PATH

echo 'write output file FULL path (add .vcf at the end)'

read OUTPUT_FILE_PATH

samtools mpileup -Ou -f ~/human_reference_genome_index_for_bowtie2_name-hg38/Homo_sapiens.GRCh38.dna.toplevel.fa $BAM_FILE_PATH | \
 bcftools call --threads 64 -mv -Ov -o $OUTPUT_FILE_PATH
