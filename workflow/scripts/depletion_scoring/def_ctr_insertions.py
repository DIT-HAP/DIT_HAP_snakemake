#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "loguru",
# ]
# ///

"""
Control Insertion Selection for Depletion Analysis
==================================================

Select control insertions from a transposon insertion dataset for use as a
neutral baseline in depletion analysis. Control insertions are drawn from
intergenic regions that sit far enough from any gene boundary to be considered
unaffected by selection pressure, providing a stable reference against which
depleted insertions can be scored.

The selection algorithm queries the annotation table for insertions annotated
as ``Intergenic region`` whose distance to both the upstream and downstream
region boundaries exceeds a fixed threshold (500 bp). The resulting set is then
intersected with the count matrix index so that only insertions actually present
in the counts are retained, and duplicate index entries are collapsed (keeping
the first occurrence). Retention statistics are logged for reproducibility.

Input
-----
- A tab-separated insertion count matrix with a 4-level row MultiIndex
  (columns 0-3) and a 2-level column header (rows 0 and 1).
- A tab-separated genomic annotation table with a matching 4-level row
  MultiIndex, containing at least the ``Type``, ``Distance_to_region_start``,
  and ``Distance_to_region_end`` columns.

Output
------
- A tab-separated file of selected control insertions, written with the row
  MultiIndex preserved (``sep="\\t"``, ``index=True``).

Usage
-----
    python def_ctr_insertions.py -i insertion_counts.tsv -a genomic_annotations.tsv -o control_insertions.tsv
    python def_ctr_insertions.py --input counts.tsv --annotation annotations.tsv --output controls.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
CONTROL_DISTANCE_THRESHOLD = 500  # Minimum distance (bp) from gene boundaries for control insertions

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input and output file paths for the selection pipeline."""
    input_file: Path
    annotation_file: Path
    output_file: Path

    def __post_init__(self) -> None:
        """Validate that input files exist, then ensure the output directory exists."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        if not self.annotation_file.exists():
            raise ValueError(f"Input file does not exist: {self.annotation_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class ControlSelectionResult:
    """Summary statistics of the control insertion selection run."""
    total_insertions_processed: int
    control_insertions_selected: int
    success_rate: float

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

setup_logger()

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_and_preprocess_data(counts_file: Path, annotations_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and preprocess insertion count and genomic annotation tables."""
    counts_df = pd.read_csv(
        counts_file, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t"
    )

    # Remove rows with any NA value
    counts_df = counts_df.dropna(axis=0, how="any").copy()

    # Load and process annotations
    insertion_annotations = pd.read_csv(
        annotations_file, index_col=[0, 1, 2, 3], sep="\t"
    )

    return counts_df, insertion_annotations


@logger.catch
def get_control_insertions(counts_df: pd.DataFrame, insertion_annotations: pd.DataFrame) -> pd.DataFrame:
    """Select control insertions based on stringent genomic criteria."""
    ctr_insertions = insertion_annotations.query(
        f"Type == 'Intergenic region' and Distance_to_region_start > {CONTROL_DISTANCE_THRESHOLD} and Distance_to_region_end > {CONTROL_DISTANCE_THRESHOLD}"
    )

    ctr_insertions = ctr_insertions[ctr_insertions.index.isin(counts_df.index)].drop_duplicates(keep="first")

    return ctr_insertions


@logger.catch
def save_results(control_insertions: pd.DataFrame, output_file: Path) -> None:
    """Save selected control insertions to a tab-separated file."""
    control_insertions.to_csv(output_file, sep="\t", index=True)
    logger.info(f"Saved {len(control_insertions)} control insertions to {output_file}")


@logger.catch
def process_control_insertions(config: InputOutputConfig) -> ControlSelectionResult:
    """Execute the complete control insertion selection pipeline."""
    logger.info(f"Starting processing of {config.input_file}")

    # Load data
    counts_df, insertion_annotations = load_and_preprocess_data(
        config.input_file, config.annotation_file
    )

    total_insertions = len(counts_df)
    logger.info(f"Loaded {total_insertions} insertions for processing")

    # Select control insertions
    control_insertions = get_control_insertions(counts_df, insertion_annotations)
    control_count = len(control_insertions)

    # Save results
    save_results(control_insertions, config.output_file)

    # Calculate success metrics
    success_rate = (control_count / total_insertions * 100) if total_insertions > 0 else 0.0

    return ControlSelectionResult(
        total_insertions_processed=total_insertions,
        control_insertions_selected=control_count,
        success_rate=success_rate,
    )

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Select control insertions for transposon depletion analysis"
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to tab-separated insertion count table")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Path to tab-separated genomic annotation table")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output path for selected control insertions")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate paths, select control insertions, and report metrics."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        # Validate input and output paths
        config = InputOutputConfig(
            input_file=args.input,
            annotation_file=args.annotation,
            output_file=args.output,
        )

        # Run the core analysis
        results = process_control_insertions(config)

        logger.success(f"Analysis complete. Results: {json.dumps(asdict(results))}")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
