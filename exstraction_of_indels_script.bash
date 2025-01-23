echo "input_bam_file"
read $BAM_FILE
echo "region (ex: chr7:55,049,219-55,049,440 or chr9:107,422,261-107,422,454 or chr14:102,085,631-102,085,812 or chr10:97,641,214-97,641,375)"
read $REGION

samtools view -h -q 20 $BAM_FILE $REGION \
| awk '$1 ~ "^@" || $6 ~ "I|D"' \
| samtools view -b - > output_indels.bam
