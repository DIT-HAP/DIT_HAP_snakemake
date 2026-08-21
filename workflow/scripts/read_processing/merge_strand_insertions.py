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
Merge PBL and PBR Strand Insertions
===================================

Merges insertion count data from PBL (left primer) and PBR (right primer)
transposon reads into a single strand-resolved table. The two inputs are
joined on chromosome and coordinate; the strand-specific columns are then
recombined so that each coordinate yields a "+" strand row (the PBL "-"
column paired with the PBR "+" column) and a "-" strand row (the PBL "+"
column paired with the PBR "-" column), together with a summed total read
count.

Input
-----
- PBL insertion TSV (-i): tab-separated with a header row and columns
  Chr, Coordinate, "+", "-".
- PBR insertion TSV (-j): tab-separated with a header row and columns
  Chr, Coordinate, "+", "-".

Output
------
- Merged TSV (-o) indexed by (Chr, Coordinate, Strand) with integer columns
  PBL, PBR, and Reads (= PBL + PBR).

Usage
-----
    python merge_strand_insertions.py -i pbl_insertions.tsv -j pbr_insertions.tsv -o merged.tsv
    python merge_strand_insertions.py -i pbl_insertions.tsv -j pbr_insertions.tsv -o merged.tsv --verbose

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
from read_processing.merge import MergeResult, merge_insertion_data  # noqa: E402

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths for the merge step."""
    input_pbl: Path
    input_pbr: Path
    output_file: Path

    def __post_init__(self) -> None:
        for path in (self.input_pbl, self.input_pbr):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


setup_logger()

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Merge PBL and PBR insertion counts by coordinate and strand")
    parser.add_argument("-i", "--inputPBL", type=Path, required=True, help="Path to the PBL input file")
    parser.add_argument("-j", "--inputPBR", type=Path, required=True, help="Path to the PBR input file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Merge PBL/PBR strand insertions and report summary statistics."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = InputOutputConfig(
            input_pbl=args.inputPBL,
            input_pbr=args.inputPBR,
            output_file=args.output,
        )

        logger.info(f"Starting processing of PBL: {config.input_pbl} and PBR: {config.input_pbr}")

        # Load data
        pbl_data = pd.read_csv(config.input_pbl, sep="\t", header=0)
        pbr_data = pd.read_csv(config.input_pbr, sep="\t", header=0)

        # Merge data
        merged_df = merge_insertion_data(pbl_data, pbr_data)

        # Save results
        logger.info(f"Writing merged data to: {config.output_file}")
        merged_df.to_csv(config.output_file, sep="\t", index=True, header=True)

        total_reads = merged_df["Reads"].sum()
        pbl_reads = merged_df["PBL"].sum()
        pbr_reads = merged_df["PBR"].sum()

        logger.info(f"Total reads: {total_reads:,}")
        logger.info(f"PBL reads: {pbl_reads:,} ({pbl_reads/total_reads*100:.1f}%)")
        logger.info(f"PBR reads: {pbr_reads:,} ({pbr_reads/total_reads*100:.1f}%)")
        logger.success(f"Output saved to: {config.output_file}")

        # Calculate results summary
        results = MergeResult(
            total_sites_processed=len(merged_df),
            total_reads_merged=total_reads,
            coordinate_strand_pairs=len(merged_df),
        )

        logger.info("Processing completed:")
        logger.info(f"  - Total sites processed: {results.total_sites_processed:,}")
        logger.info(f"  - Total reads merged: {results.total_reads_merged:,}")
        logger.info(f"  - Coordinate-strand pairs: {results.coordinate_strand_pairs:,}")

        logger.success(f"Analysis complete. Results saved to {config.output_file}")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
