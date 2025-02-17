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

#!/bin/bash

set -e  # Exit on error

# Debugging: Print the user running the script
echo "Running as user: $(whoami)"
echo "Current working directory: $(pwd)"

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

# Debugging: Check if input BAM file exists and is readable
if [ ! -f "$input_bam" ]; then
  echo "ERROR: Input BAM file does not exist: $input_bam"
  exit 1
elif [ ! -r "$input_bam" ]; then
  echo "ERROR: Input BAM file is not readable: $input_bam"
  ls -l "$input_bam"
  exit 1
fi

# Debugging: Check output directory permissions
if [ ! -d "$output_dir" ]; then
  echo "ERROR: Output directory does not exist: $output_dir"
  exit 1
elif [ ! -w "$output_dir" ]; then
  echo "ERROR: Output directory is not writable: $output_dir"
  ls -ld "$output_dir"
  exit 1
fi

# Debugging: Print file ownership and permissions
echo "File permissions for input BAM:"
ls -l "$input_bam"
echo "Directory permissions for output:"
ls -ld "$output_dir"

# Create output directory if it does not exist
mkdir -p "$output_dir"

# Determine the total number of reads in the BAM file
total_reads=$(samtools view -c "$input_bam")
if [ $? -ne 0 ]; then
  echo "ERROR: Failed to count reads in BAM file. Check samtools installation and permissions."
  exit 1
fi

reads_per_chunk=$(( (total_reads + num_chunks - 1) / num_chunks ))  # Round up division

echo "Total reads in BAM: $total_reads"
echo "Splitting into $num_chunks chunks, each containing approximately $reads_per_chunk reads."

# Split BAM file into chunks
for ((i=0; i<num_chunks; i++)); do
  start=$(( i * reads_per_chunk ))
  end=$(( start + reads_per_chunk ))

  output_chunk="${output_dir}/${base_name}_chunk_${i}.bam"

  # Debugging: Check if the output file is writable before writing
  touch "$output_chunk" 2>/dev/null || { echo "ERROR: Cannot write to $output_chunk"; exit 1; }

  # Extract reads and save to chunk
  samtools view -b -@ "$threads" "$input_bam" | \
    samtools split -u -f "${output_dir}/${base_name}_chunk_%d.bam" - "$num_chunks"

  echo "Created: $output_chunk"
done

echo "BAM splitting complete: Chunks stored in $output_dir"

# **Index all BAM files using the specified number of threads**
echo "Indexing BAM files with $threads threads..."
for bam_file in "$output_dir"/*.bam; do
  if [ -s "$bam_file" ]; then  # Check if file exists and is non-empty
    echo "Indexing $bam_file..."
    samtools index -@ "$threads" "$bam_file"
  else
    echo "WARNING: Skipping empty BAM file: $bam_file"
  fi
done

echo "All BAM files have been indexed successfully!"

