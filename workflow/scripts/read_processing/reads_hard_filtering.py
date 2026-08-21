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
Hard Filtering of Insertion Reads
=================================

Filter insertion reads by a minimum read-count threshold measured at a chosen
initial timepoint. Each sample is processed independently: an insertion is
retained for a sample only when its read count at the initial timepoint meets
or exceeds the cutoff, so a sparsely-covered insertion in one sample does not
suppress the same locus in a well-covered sample.

The core algorithm groups the wide, multi-indexed count matrix by ``Sample``
(the first column level), builds a boolean mask from the ``>= cutoff`` test on
the initial-timepoint column, keeps the masked rows per sample, and concatenates
the surviving per-sample blocks back into a single matrix. Retention statistics
(total, retained, removed, retention rate, samples processed) are logged for
reproducibility.

Input
-----
- A tab-separated counts matrix with a 4-level row MultiIndex (columns 0-3) and
  a 2-level column MultiIndex ``(Sample, Timepoint)`` on header rows 0 and 1.

Output
------
- A tab-separated matrix of retained insertions, written with the row index and
  the 2-level column header preserved (``sep="\\t"``, ``header=True``,
  ``index=True``). Same structure as the input, filtered by row.

Usage
-----
    python reads_hard_filtering.py -i raw_reads.tsv -o filtered_reads.tsv -itp 0h -c 5
    python reads_hard_filtering.py -i counts.tsv -o filtered.tsv --init-timepoint YES0 --cutoff 10 --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402
from read_processing.hard_filtering import (  # noqa: E402
    AnalysisResult,
    load_insertion_data,
    validate_timepoint_exists,
    apply_hard_filtering,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Configuration for hard filtering parameters and paths."""
    input_file: Path
    output_file: Path
    initial_timepoint: str
    cutoff_threshold: int

    def __post_init__(self) -> None:
        """Validate paths and timepoint, then ensure the output directory exists."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        if not self.initial_timepoint.strip():
            raise ValueError("Initial timepoint cannot be empty")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


setup_logger()

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Filter insertion reads by hard filtering based on read count thresholds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reads_hard_filtering.py -i raw_reads.tsv -o filtered_reads.tsv -itp 0h -c 5
  python reads_hard_filtering.py -i counts.tsv -o filtered.tsv --init-timepoint YES0 --cutoff 10
        """
    )

    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input TSV file with insertion reads"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output TSV file for filtered reads"
    )
    parser.add_argument(
        "-itp", "--init-timepoint",
        type=str,
        required=True,
        help="Initial timepoint column name for filtering"
    )
    parser.add_argument(
        "-c", "--cutoff",
        type=int,
        required=True,
        help="Minimum read count threshold at initial timepoint"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    return parser.parse_args()


def main() -> int:
    """Main orchestrator: load, filter, and write hard-filtered insertion reads."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output,
            initial_timepoint=args.init_timepoint.strip(),
            cutoff_threshold=args.cutoff,
        )

        logger.info(f"Starting processing of {config.input_file}")

        df = load_insertion_data(config.input_file)
        filtered_df, results = apply_hard_filtering(
            df, config.initial_timepoint, config.cutoff_threshold
        )

        logger.info("Saving filtered results...")
        filtered_df.to_csv(config.output_file, sep="\t", header=True, index=True)

        logger.info("=" * 70)
        logger.info("FILTERING SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Total insertions: {results.total_insertions:,}")
        logger.info(f"Retained insertions: {results.retained_insertions:,}")
        logger.info(f"Removed insertions: {results.removed_insertions:,}")
        logger.success(f"Retention rate: {results.retention_rate:.2f}%")
        logger.info(f"Samples processed: {results.samples_processed}")
        logger.success(f"Analysis complete. Results saved to {config.output_file}")

    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
