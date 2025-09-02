#!/bin/bash

# Ensure output directories exist
mkdir -p ./results/Pavlova_mRNA_2020_cos1
mkdir -p ./results/Pavlova_mRNA_2020_cos2

# Check and set ulimit for open files
ulimit -n 4096

# Map first single-end read file
STAR --runThreadN 24 \
--readFilesIn ./data/RNA_reads_controls/200207_HSGA.Pavlova_mRNA_2020.cos_contr1.fastq \
--genomeDir ./data/STAR_index \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./results/Pavlova_mRNA_2020_cos1/

# Map second single-end read file
STAR --runThreadN 24 \
--readFilesIn ./data/RNA_reads_controls/200207_HSGA.Pavlova_mRNA_2020.cos_contr2.fastq \
--genomeDir ./data/STAR_index \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./results/Pavlova_mRNA_2020_cos2/
