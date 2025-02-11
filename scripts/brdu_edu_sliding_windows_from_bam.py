#!/usr/bin/env python3
"""
Script Name: brdu_edu_sliding_windows_from_bam.py

Description:
This script processes a BAM file to compute probability-weighted BrdU and EdU 
modification frequencies using a sliding window approach. It outputs BEDGRAPH 
files for visualization. The number of reads to process (from start of file) 
or a read query name can be provided to produce a bedgraph for a subset of 
the data. If num_reads or query_name are omitted all reads will be processed.

Usage:
    python brdu_edu_sliding_windows_from_bam.py -b <input.bam> -o1 <brdu_output.bedgraph> -o2 <edu_output.bedgraph>

Arguments:
    -b, --bam          Path to input BAM file (required)
    -o1, --brdu_output Path to output BrdU BEDGRAPH file
    -o2, --edu_output  Path to output EdU BEDGRAPH file
    -p, --processes    Number of CPU cores to use (default: auto-detect)
    -n, --num_reads    Number of reads to process (default: all)
    -w, --window_bp    Window size in base pairs (default: 100)
    -s, --step_size_bp Step size in base pairs (default: 10)
    -r, --query_name   Process only the specified read name

Author:
    Chris Sansam

Date:
    February 2025
"""

import pysam
import argparse
import sys
import os
import itertools
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

### Auto-Detect Project Root ###
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# Ensure `lib/` is in `sys.path`
LIB_PATH = os.path.join(SCRIPT_DIR, "lib")
sys.path.insert(0, LIB_PATH)

# Import helper functions
import slidingWindowsBedgraphs_utils as sw


# Function: Process a Single Read
def process_read(read_data):
    """
    Processes extracted read data.

    Parameters:
    - read_data: Tuple containing (reference_name, reference_start, reference_end,
      is_reverse, MM_tag, ML_tag, window_size, step_size)

    Returns:
    - List of tuples (chrom, start, end, BrdU_freq, EdU_freq)
    """
    if read_data is None:
        return []

    (
        reference_name,
        reference_start,
        reference_end,
        is_reverse,
        MM,
        ML,
        window_size,
        step_size,
    ) = read_data

    # Parse modifications and compute window frequencies
    mod_positions = sw.parse_ml_mm(MM, ML)
    window_results = sw.compute_sliding_windows(mod_positions, window_size, step_size)

    return sw.convert_relative_to_abs_positions(
        reference_name, reference_start, reference_end, is_reverse, window_results
    )


# Function: Process BAM & Write BEDGRAPH
def process_bam_and_write_bedgraphs(
    bam_path,
    brdu_output,
    edu_output,
    window_size,
    step_size,
    query_name=None,
    num_processes=None,
    num_reads=None,
):
    """
    Uses multiprocessing to process BAM reads in parallel and writes
    probability-weighted modification frequencies to BrdU and EdU BEDGRAPH files.

    Parameters:
    - bam_path: Path to the input BAM file.
    - brdu_output: Path to the BrdU BEDGRAPH file.
    - edu_output: Path to the EdU BEDGRAPH file.
    - window_size: Size of window in bp for calculating weighted incorporation rate.
    - step_size: Step size in bp for sliding window.
    - query_name: Read query name to process only a specific read.
    - num_processes: Number of CPU processes to use (default: auto-detect).
    - num_reads: Number of reads to process (default: all).
    """
    if num_processes is None:
        num_processes = max(1, cpu_count() - 1)  # Use all but one CPU core

    print(f"Using {num_processes} CPU cores for multiprocessing...")

    # Ensure output directories exist
    os.makedirs(os.path.dirname(brdu_output), exist_ok=True)
    os.makedirs(os.path.dirname(edu_output), exist_ok=True)

    # Open BAM file and fetch reads
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        if query_name:
            # Fetch only the read with matching query_name
            read_chunks = [read for read in bamfile if read.query_name == query_name]
            total_reads = len(read_chunks)
            if total_reads == 0:
                raise ValueError(
                    f"Error: Read with query_name '{query_name}' not found in BAM file."
                )
        else:
            read_chunks = bamfile.fetch()
            total_reads = bamfile.mapped if num_reads is None else num_reads

        # Extract only necessary data to avoid pickling issues
        read_data_list = (
            sw.extract_read_data(read, window_size, step_size) for read in read_chunks
        )
        if num_reads:
            read_data_list = itertools.islice(
                read_data_list, num_reads
            )  # Limit to first `num_reads`

        # Process reads in parallel
        with Pool(processes=num_processes) as pool:
            results = list(
                tqdm(
                    pool.imap(process_read, read_data_list),
                    total=total_reads,
                    desc="Processing BAM Reads",
                    unit="read",
                )
            )

    # Flatten results
    all_window_results = [item for sublist in results for item in sublist]

    # Write final BEDGRAPH files
    sw.write_bedgraphs(brdu_output, edu_output, all_window_results)


# Main Execution (CLI)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process BAM file and generate BrdU & EdU BEDGRAPH files."
    )

    # Required Arguments
    parser.add_argument(
        "-b", "--bam", required=True, help="Path to input BAM file"
    )

    # Output Files
    parser.add_argument(
        "-o1",
        "--brdu_output",
        default=os.path.join(PROJECT_ROOT, "results", "BrdU.bedgraph"),
        help="Output file for BrdU BEDGRAPH",
    )
    parser.add_argument(
        "-o2",
        "--edu_output",
        default=os.path.join(PROJECT_ROOT, "results", "EdU.bedgraph"),
        help="Output file for EdU BEDGRAPH",
    )

    # Processing Parameters
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=None,
        help="Number of CPU cores to use (default: auto-detect)",
    )
    parser.add_argument(
        "-n",
        "--num_reads",
        type=int,
        default=None,
        help="Number of reads to process (default: all)",
    )

    # Sliding Window Parameters
    parser.add_argument(
        "-w",
        "--window_bp",
        type=int,
        default=100,
        help="Width of windows in base pairs",
    )
    parser.add_argument(
        "-s",
        "--step_size_bp",
        type=int,
        default=10,
        help="Step size of windows in base pairs",
    )

    # Specific Read Processing
    parser.add_argument(
        "-r",
        "--query_name",
        default=None,
        help="Read query name to return BEDGRAPH for that read only.",
    )

    args = parser.parse_args()

    # Run processing function
    process_bam_and_write_bedgraphs(
        args.bam,
        args.brdu_output,
        args.edu_output,
        args.window_bp,
        args.step_size_bp,
        args.query_name,
        args.processes,
        args.num_reads,
    )