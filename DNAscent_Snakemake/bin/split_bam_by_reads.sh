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
#   2025-02-17
#

#!/bin/bash

#!/bin/bash

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

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Get the total number of reads in the BAM file
total_reads=$(samtools view -c "$input_bam")
reads_per_chunk=$(( (total_reads + num_chunks - 1) / num_chunks ))  # Round up division

echo "Total reads in BAM: $total_reads"
echo "Splitting into $num_chunks chunks, each containing approximately $reads_per_chunk reads."

# Extract BAM header once
samtools view -H "$input_bam" > "${output_dir}/${base_name}_header.sam"

# Split BAM file into chunks
for ((i=0; i<num_chunks; i++)); do
  start=$(( i * reads_per_chunk ))
  end=$(( start + reads_per_chunk ))

  output_chunk="${output_dir}/${base_name}_chunk_${i}.bam"

  echo "Processing chunk $i ($start to $end)..."

  # Extract reads and save to chunk
  samtools view -@ "$threads" "$input_bam" | awk -v s="$start" -v e="$end" 'NR > s && NR <= e' | \
    cat "${output_dir}/${base_name}_header.sam" - | samtools view -b -@ "$threads" > "$output_chunk"

  echo "Created: $output_chunk"
done

echo "BAM splitting complete: Chunks stored in $output_dir"

# Index all BAM files
echo "Indexing BAM files with $threads threads..."
for bam_file in "$output_dir"/*.bam; do
  if [ -s "$bam_file" ]; then
    echo "Indexing $bam_file..."
    samtools index -@ "$threads" "$bam_file"
  else
    echo "WARNING: Skipping empty BAM file: $bam_file"
  fi
done

echo "All BAM files have been indexed successfully!"
