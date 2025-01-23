#!bin/bash

# Read the numbers from the two files
num_with_indels=$(cat number_of_reads_with_indels)
total_num_reads=$(cat total_number_of_reads)

# Calculate the division using the bc command
result=$(echo "scale=5; $num_with_indels / $total_num_reads" | bc)

# Print the result to the terminal
echo "NHEJ percentage: $result"

# Save the result to a file
echo "$result" > NHEJ_percentage
