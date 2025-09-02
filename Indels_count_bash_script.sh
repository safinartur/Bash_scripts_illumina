#!/bin/bash

# INPUTS:
# - THREADS: Number of threads to use for Bowtie2 and other commands (int)
# - BOWTIE2_INDEX: Path to the Bowtie2 index to be used (string)
# - FASTQ_1: Path to the first FASTQ file (string)
# - FASTQ_2: Path to the second FASTQ file (string)
# - OUTPUT_DIR: Directory to save output files (string)
# - OUTPUT_FILE_NAME: Base name for output files (string)
# - REGION: The region of interest in the format "chrX:start-end" (string)
mkdir "${5}"

bowtie2 --threads $1 -x $2 --sensitive -1 $3 -2 $4 -S "${5}/${6}.sam"


samtools view -@ $1 -S -b -F 4 "${5}/${6}.sam" | samtools sort -@ $1 > "${5}/${6}_sorted.bam"


samtools index -@ $1 "${5}/${6}_sorted.bam"

#MQ > 20, remove header lines, remove reads without indels
samtools view -h -q 20 "${5}/${6}_sorted.bam" | awk '$1 ~ "^@" || $6 ~ "I|D"' > "${5}/temp.sam"


samtools view -b "${5}/temp.sam" > "${5}/output_indels.bam"

samtools view -c "${5}/output_indels.bam" > "${5}/number_of_reads_with_indels"


samtools view -c "${5}/${6}_sorted.bam" > "${5}/total_number_of_reads"

rm "${5}/temp.sam"
rm "${5}/${6}.sam"
