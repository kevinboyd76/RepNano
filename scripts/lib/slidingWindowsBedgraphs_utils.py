"""
Utilities for processing BAM files with DNA modifications.

This module provides helper functions for parsing DNA modification tags, computing 
modification frequencies using a sliding window approach, converting relative 
positions to genomic coordinates, and exporting results in BEDGRAPH format.

Functions:
----------
- parse_ml_mm(MM, ML): Parses MM and ML tags from BAM files to extract BrdU and EdU 
  modification probabilities at each base.
- extract_read_data(read, window_size, step_size): Extracts relevant information 
  (such as modification tags) from a single BAM read.
- compute_sliding_windows(mod_positions, window_size=100, step_size=10): Computes 
  probability-weighted BrdU and EdU modification frequencies in a sliding window 
  manner across a sequence.
- convert_relative_to_abs_positions(reference_name, reference_start, reference_end, 
  is_reverse, window_results): Converts relative modification positions from 
  sliding windows into absolute genomic positions.
- write_bedgraphs(brdu_output_file, edu_output_file, new_windows): Writes the 
  modification frequencies into BEDGRAPH files for visualization.
  
Usage:
------
This module is intended to be imported and used within a BAM processing pipeline. 
Example usage within a script:

    from slidingWindowsBedgraphs_utils import parse_ml_mm, compute_sliding_windows

    mod_positions = parse_ml_mm(MM_tag, ML_tag)
    window_results = compute_sliding_windows(mod_positions)

Author:
-------
Chris Sansam

Date:
-----
February 2025
"""


def parse_ml_mm(MM, ML):
	"""
	Parses MM and ML tags, correctly handling relative positioning and matching ML values.

	Returns:
	- mod_positions: Dictionary {position: {'BrdU': prob, 'EdU': prob, 'None': prob}}
	"""
	mod_positions = {}
	ml_probs = iter(ML)  # Convert ML array to iterator

	# Parse MM tag
	mm_entries = [entry for entry in MM.split(";") if entry]

	for entry in mm_entries:
		parts = entry.split(",")
		mod_type = parts[0]  # "N+b?" or "N+e?"
		rel_positions = list(map(int, parts[1:]))  # Extract relative positions

		abs_position = 0  # Start at position 0, accumulate to get absolute positions

		for rel_pos in rel_positions:
			abs_position += rel_pos + 1  # **Fix: Add +1 to account for correct base spacing**

			# Get probability from ML, convert 0-256 → float 0-1
			prob = next(ml_probs) / 256.0

			if abs_position not in mod_positions:
				mod_positions[abs_position] = {"BrdU": 0, "EdU": 0, "None": 1}

			if mod_type == "N+b?":
				mod_positions[abs_position]["BrdU"] = prob
			elif mod_type == "N+e?":
				mod_positions[abs_position]["EdU"] = prob

	# Compute probability of *no modification* after processing all positions
	for pos in mod_positions:
		total_mod_prob = mod_positions[pos]["BrdU"] + mod_positions[pos]["EdU"]
		mod_positions[pos]["None"] = max(0, 1 - total_mod_prob)  # Ensure no negative values

	return mod_positions

def extract_read_data(read,window_size,step_size):
	"""
	Extracts necessary information from a pysam AlignedSegment object.

	Returns:
	- Tuple containing (reference_name, reference_start, is_reverse, MM_tag, ML_tag)
	"""
	if not read.has_tag("MM") or not read.has_tag("ML"):
		return None  # Skip reads without modification data

	return (
		read.reference_name,
		read.reference_start,
		read.reference_end,
		read.is_reverse,
		read.get_tag("MM"),
		read.get_tag("ML"),
		window_size,
		step_size
	)

def compute_sliding_windows(mod_positions, window_size=100, step_size=10):
	"""
	Computes probability-weighted modification frequencies for BrdU and EdU
	using a sliding window approach.

	Parameters:
	- mod_positions: Dictionary {position: {'BrdU': prob, 'EdU': prob, 'None': prob}}
	- window_size: Size of the sliding window (default=100 bp)
	- step_size: Step size for the sliding window (default=10 bp)

	Returns:
	- window_results: List of tuples (window_start, window_end, BrdU_freq, EdU_freq)
	"""

	# Get all modification positions and sort them
	sorted_positions = sorted(mod_positions.keys())

	if not sorted_positions:
		return []  # Return empty if no modifications exist

	min_pos = sorted_positions[0]
	max_pos = sorted_positions[-1]

	window_results = []

	# Slide windows along the entire sequence range
	for window_start in range(min_pos, max_pos - window_size + 1, step_size):
		window_end = window_start + window_size
		brdu_sum = 0.0
		edu_sum = 0.0
		t_count = 0

		for pos in range(window_start, window_end):
			if pos in mod_positions:
				brdu_sum += mod_positions[pos].get("BrdU", 0)  # BrdU probability sum
				edu_sum += mod_positions[pos].get("EdU", 0)    # EdU probability sum
				t_count += 1  # Assuming all modifications occur at T bases

		if t_count > 0:
			brdu_freq = brdu_sum / t_count
			edu_freq = edu_sum / t_count
		else:
			brdu_freq = 0.0
			edu_freq = 0.0

		window_results.append((window_start, window_end, brdu_freq, edu_freq))

	return window_results



def convert_relative_to_abs_positions(reference_name, reference_start, reference_end, is_reverse, window_results):
	"""
	Converts relative positions from sliding window output into absolute genomic positions.

	Parameters:
	- reference_name: Chromosome or read name
	- reference_start: Start position of the read
	- is_reverse: Boolean indicating strand direction
	- window_results: List of tuples (relative_start, relative_end, BrdU_freq, EdU_freq)

	Returns:
	- List of tuples (chrom, start, end, BrdU_freq, EdU_freq)
	"""

	pos_list = tuple(int(tup[0]) for tup in window_results)
	window_size = (window_results[0])[1] - (window_results[0])[0]

	if not is_reverse:
		starts = tuple(reference_start + 1 + x for x in pos_list)
	else:
		starts = tuple(reference_end - 1 - x for x in pos_list)

	new_windows = tuple(
		(reference_name, x, x + window_size, window_results[i][2], window_results[i][3])
		for i, x in enumerate(starts)
	)

	return new_windows



def write_bedgraphs(brdu_output_file, edu_output_file, new_windows):
		"""
		Writes probability-weighted modification frequencies to a BEDGRAPH file.

		Parameters:
		- output_file: Path to the BEDGRAPH file.
		- new_windows: List of tuples (chrom_name, start, end, BrdU_freq, EdU_freq).
		- value_index: Index of the value to write (3 for BrdU, 4 for EdU).
		"""
		with open(brdu_output_file, "w") as f:
				for window in new_windows:
						chrom, start, end, brdu_freq, edu_freq = window
						value = brdu_freq
						f.write(f"{chrom}\t{start}\t{end}\t{value:.6f}\n")
		print(f"BEDGRAPH file saved: {brdu_output_file}")

		with open(edu_output_file, "w") as f:
				for window in new_windows:
						chrom, start, end, brdu_freq, edu_freq = window
						value = edu_freq
						f.write(f"{chrom}\t{start}\t{end}\t{value:.6f}\n")
		print(f"BEDGRAPH file saved: {edu_output_file}")

def foo():
	print("foo!")
		