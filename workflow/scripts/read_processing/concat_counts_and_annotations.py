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
Concatenate Per-Sample Counts and Annotations
==============================================

Reads multiple per-sample insertion count TSVs (one per sample/condition
combination) and multiple annotation TSVs, then concatenates them into a single
wide count matrix and a single deduplicated annotation table.

The count matrix carries a (Sample, Timepoint) MultiIndex on its columns, keyed
by the leading dot-delimited token of each input file name. The annotation table
is the union of all per-sample annotations with duplicate rows removed, indexed
by the (Chr, Coordinate, Strand, Target) insertion coordinate.

Input
-----
- Per-sample count TSVs with a (Chr, Coordinate, Strand, Target) row MultiIndex.
- Per-sample annotation TSVs with the same row MultiIndex.

Output
------
- Combined counts TSV (wide, (Sample, Timepoint) column MultiIndex).
- Deduplicated annotations TSV.

Usage
-----
    python concat_counts_and_annotations.py \
        -i s1_PBL.tsv s1_EMM.tsv \
        -a s1_PBL.annotated.tsv s1_EMM.annotated.tsv \
        -oc raw_reads.tsv -oa annotations.tsv

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
    """Configuration for concatenating per-sample counts and annotations."""
    counts_files: list[str]
    annotation_files: list[str]
    output_counts: Path
    output_annotations: Path

    def __post_init__(self) -> None:
        """Validate all inputs exist and create output directories."""
        for f in [Path(x) for x in self.counts_files + self.annotation_files]:
            if not f.exists():
                raise ValueError(f"Input file does not exist: {f}")
        self.output_counts.parent.mkdir(parents=True, exist_ok=True)
        self.output_annotations.parent.mkdir(parents=True, exist_ok=True)

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

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def concatenate(config: Config) -> None:
    """Concatenate per-sample count and annotation TSV files."""
    counts_df = {}
    for f in config.counts_files:
        key = Path(f).name.split(".")[0]
        counts_df[key] = pd.read_csv(f, sep="\t", index_col=[0, 1, 2, 3])
        logger.info(f"Loaded counts {key}: {counts_df[key].shape[0]:,} rows")
    counts = pd.concat(counts_df, axis=1).rename_axis(["Sample", "Timepoint"], axis=1)
    logger.info(f"Combined counts shape: {counts.shape[0]:,} rows × {counts.shape[1]} columns")

    annotations_df = {}
    for f in config.annotation_files:
        key = Path(f).name.split(".")[0]
        annotations_df[key] = pd.read_csv(f, index_col=[0, 1, 2, 3], sep="\t")
        logger.info(f"Loaded annotations {key}: {annotations_df[key].shape[0]:,} rows")
    annotations = (
        pd.concat(list(annotations_df.values()), axis=0)
        .reset_index()
        .drop_duplicates()
        .set_index(["Chr", "Coordinate", "Strand", "Target"])
    )

    counts.to_csv(config.output_counts, sep="\t", index=True)
    annotations.to_csv(config.output_annotations, sep="\t", index=True)
    logger.success(f"Counts: {counts.shape[0]:,} rows → {config.output_counts}")
    logger.success(f"Annotations: {annotations.shape[0]:,} rows → {config.output_annotations}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Concatenate per-sample insertion counts and annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python concat_counts_and_annotations.py \\
      -i s1_PBL.tsv s1_EMM.tsv -a s1_PBL.annotated.tsv s1_EMM.annotated.tsv \\
      -oc raw_reads.tsv -oa annotations.tsv
        """,
    )
    parser.add_argument("-i", "--counts", nargs="+", type=Path, required=True,
                        help="Per-sample count TSV files")
    parser.add_argument("-a", "--annotations", nargs="+", type=Path, required=True,
                        help="Per-sample annotation TSV files")
    parser.add_argument("-oc", "--output-counts", type=Path, required=True,
                        help="Output combined counts TSV")
    parser.add_argument("-oa", "--output-annotations", type=Path, required=True,
                        help="Output combined annotations TSV")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()

def main() -> int:
    """Main orchestrator function for the script execution."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = Config(
            counts_files=[str(f) for f in args.counts],
            annotation_files=[str(f) for f in args.annotations],
            output_counts=args.output_counts,
            output_annotations=args.output_annotations,
        )

        logger.info(f"Concatenating {len(config.counts_files):,} count files and "
                    f"{len(config.annotation_files):,} annotation files")
        concatenate(config)
        logger.success("Script completed successfully!")

    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
