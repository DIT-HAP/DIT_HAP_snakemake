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
Extract Read-Pair Filtering Statistics from Log Files
=====================================================

Parses log files produced during read-pair filtering to extract per-sample
filtering statistics for PBL (left) and PBR (right) read pairs. For each log
it reads the "FILTERING SUMMARY" blocks, pulls out chunk counts, original /
filtered / removed read-pair counts and retention rates, then aggregates the
PBL and PBR figures into combined totals and an overall retention rate.

The parser matches summary blocks of the following shape (one per read type)::

    ============================================================
    FILTERING SUMMARY
    ============================================================
    Total chunks processed: 12345
    Original read pairs: 1,234,567
    Filtered read pairs: 1,123,456
    Removed read pairs: 111,111
    Overall retention rate: 91.05%
    Output written to: sample_name.PBL.filtered.parquet

Input
-----
- One or more filtering log files, each containing "FILTERING SUMMARY"
  sections whose "Output written to" line ends in ``.PBL.filtered.parquet`` or
  ``.PBR.filtered.parquet``. The sample name is taken from each log file's stem.

Output
------
- A tab-separated statistics table indexed by sample (index label ``Sample``)
  with per-read-type counts and retention rates plus the aggregated columns
  ``total_original_pairs``, ``total_filtered_pairs`` and
  ``overall_retention_rate``. Floats are written with ``%.2f`` precision.

Usage
-----
    python extract_mapping_filtering_statistics.py -i sample1.log sample2.log -o filtering_statistics.tsv
    python extract_mapping_filtering_statistics.py -i sample1.log -o stats.tsv -v

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
from qc.mapping_stats import (  # noqa: E402
    AnalysisResult,
    extract_summary_data,
    create_dataframe,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input log paths and the output statistics path."""
    input_files: list[Path]
    output_file: Path

    def __post_init__(self) -> None:
        """Validate that every input file exists and create the output parent directory."""
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)




setup_logger()

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Extract filtering statistics from read pair filtering log files"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        nargs="+",
        required=True,
        help="Path to the log files",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Path to save the output statistics TSV",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate inputs, extract statistics, and write the TSV."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        # Validate inputs
        config = InputOutputConfig(
            input_files=args.input,
            output_file=args.output,
        )

        logger.info("Starting filtering statistics extraction")

        # Extract statistics
        statistics = extract_summary_data(config.input_files)

        if not statistics:
            logger.error("No statistics extracted from any log files")
            return 1

        # Create DataFrame
        df = create_dataframe(statistics)

        # Save results (dissolved from the former save_results() helper)
        if df.empty:
            logger.error("Cannot save empty DataFrame")
        else:
            # Sort by sample name
            df = df.rename_axis("Sample", axis=0).sort_index()

            # Save to file
            selected_columns = [
                "total_original_pairs",
                "total_filtered_pairs",
                "overall_retention_rate",
                "pbl_pbr_ratio",
                "retention_rate_pbl",
                "retention_rate_pbr"
            ]
            df[selected_columns].to_csv(config.output_file, sep="\t", index=True, float_format="%.2f")
            logger.success(f"Statistics saved to: {config.output_file}")

            # Display summary
            summary_cols = ["total_original_pairs", "total_filtered_pairs", "overall_retention_rate"]
            available_cols = [col for col in summary_cols if col in df.columns]

            if available_cols:
                logger.info("Summary statistics:")
                logger.info(f"\n{df[available_cols].describe()}")

        # Create analysis result
        result = AnalysisResult(
            total_samples_processed=len(statistics),
            total_log_files=len(config.input_files),
            output_path=config.output_file,
        )

        logger.success("Extraction completed successfully!")
        logger.info(f"Processed {result.total_samples_processed} samples from {result.total_log_files} log files")
        logger.info(f"Output saved to: {result.output_path}")

    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
