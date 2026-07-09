#!/usr/bin/env python3

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
Merge Similar Timepoints by Summing Read Counts
===============================================

Merge a set of near-identical timepoint columns in a wide-format insertion
count matrix into a single consolidated column. The specified similar
timepoints are summed row-wise into one merged column, the redundant source
columns are dropped, and the remaining columns are re-sorted alphabetically.

This step collapses technical or biological replicate timepoints that should
be treated as one condition downstream, keeping the count matrix compact
while preserving the (Chr, Coordinate, Strand, Target) row MultiIndex.

Input
-----
- A tab-separated count matrix with a 4-level row MultiIndex
  (Chr, Coordinate, Strand, Target) and one column per timepoint.

Output
------
- A tab-separated matrix in the same format, with the similar timepoint
  columns replaced by a single merged column and columns sorted by name.

Usage
-----
    python merge_similar_timepoints.py -i counts.tsv -o merged.tsv \
        -s 0h YES0 -m YES0 -d 0h
    python merge_similar_timepoints.py -i counts.tsv -o merged.tsv \
        -s 0h YES0 -m YES0 -d 0h --verbose

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

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class Config:
    """Configuration for merging similar timepoints."""
    input_file: Path
    output_file: Path
    similar_timepoints: list[str]
    merged_timepoint: str
    drop_columns: list[str]

    def __post_init__(self) -> None:
        """Validate the input exists and ensure the output directory is present."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level=log_level,
        colorize=False,
    )

setup_logger()

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def merge_timepoints(config: Config) -> None:
    """Merge similar timepoints by summing read counts into one column."""
    df = pd.read_csv(config.input_file, sep="\t", index_col=[0, 1, 2, 3])
    logger.info(f"Loaded {df.shape[0]:,} rows, columns: {list(df.columns)}")

    logger.info(f"Summing {config.similar_timepoints} -> '{config.merged_timepoint}'")
    df[config.merged_timepoint] = df[config.similar_timepoints].sum(axis=1)

    logger.info(f"Dropping columns: {config.drop_columns}")
    df.drop(columns=config.drop_columns, inplace=True)
    df.sort_index(axis=1, inplace=True)

    df.to_csv(config.output_file, sep="\t", index=True)
    logger.success(f"Merged -> {df.shape[0]:,} rows, columns: {list(df.columns)}")
    logger.success(f"Saved to {config.output_file}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Merge similar timepoints by summing read counts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python merge_similar_timepoints.py \\
      -i counts.tsv -o merged.tsv \\
      -s 0h YES0 -m YES0 -d 0h
        """,
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Input counts TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output merged counts TSV file")
    parser.add_argument("-s", "--similar-timepoints", nargs="+", required=True,
                        help="Timepoint column names to sum together")
    parser.add_argument("-m", "--merged-timepoint", type=str, required=True,
                        help="Name for the merged timepoint column")
    parser.add_argument("-d", "--drop-columns", nargs="+", required=True,
                        help="Column names to drop after merging")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()

def main() -> int:
    """Merge similar timepoints for a single count matrix and write the result."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    try:
        config = Config(
            input_file=args.input,
            output_file=args.output,
            similar_timepoints=args.similar_timepoints,
            merged_timepoint=args.merged_timepoint,
            drop_columns=args.drop_columns,
        )
        logger.info(f"Merging timepoints in {config.input_file}")
        merge_timepoints(config)
        logger.success("Script completed successfully!")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
