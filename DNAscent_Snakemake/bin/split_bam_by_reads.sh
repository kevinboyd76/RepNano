#!/bin/bash

#
# BAM File Splitter and Indexer
#
# Description:
#   This script splits a BAM file into a specified number of chunks and indexes the resulting files.
#   Each chunk retains the BAM header to ensure downstream compatibility.
#
# Usage:
#   ./split_bam_by_chunks.sh <input.bam> <output_dir> <num_chunks> <threads>
#
# Arguments:
#   <input.bam>   - Path to the input BAM file.
#   <output_dir>  - Output directory.
#   <num_chunks>  - Number of chunks to split the BAM file into.
#   <threads>     - Number of threads to use for samtools operations.
#
# Requirements:
#   - samtools must be installed and accessible in the environment.
#   - Ensure sufficient disk space for output files.
#
# Author:
#   Chris Sansam
#   2025-02-06
#

set -e  # Exit on error

# Check for required arguments
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <input.bam> <output_dir> <num_chunks> <threads>"
  exit 1
fi

# Assign arguments to variables
input_bam="$1"
output_dir="$2"
num_chunks="$3"
threads="$4"

base_name=$(basename -s .bam "$input_bam")

# Create output directory
mkdir -p "$output_dir"

# Extract and store the header
header_file="${output_dir}/${base_name}_header.sam"
samtools view -H "$input_bam" > "$header_file"

# Determine the total number of reads in the BAM file
total_reads=$(samtools view -c "$input_bam")
reads_per_chunk=$(( (total_reads + num_chunks - 1) / num_chunks ))  # Round up division

echo "Total reads in BAM: $total_reads"
echo "Splitting into $num_chunks chunks, each containing approximately $reads_per_chunk reads."

# Split BAM file into chunks
for ((i=0; i<num_chunks; i++)); do
  start=$(( i * reads_per_chunk ))
  end=$(( start + reads_per_chunk ))

  output_chunk="${output_dir}/${base_name}_chunk_${i}.bam"

  # Extract reads in the specified range and combine with the header to create a valid BAM file
  samtools view -@ "$threads" "$input_bam" | awk -v start="$start" -v end="$end" 'NR > start && NR <= end' | \
    cat "$header_file" - | samtools view -b -@ "$threads" > "$output_chunk"

  echo "Created: $output_chunk"
done

echo "BAM splitting complete: Chunks stored in $output_dir"

# **Index all BAM files using the specified number of threads**
echo "Indexing BAM files with $threads threads..."
for bam_file in "$output_dir"/*.bam; do
  echo "Indexing $bam_file..."
  samtools index -@ "$threads" "$bam_file"
done

echo "All BAM files have been indexed successfully!"
