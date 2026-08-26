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
PBL-PBR Pairs Extraction
========================

Extract PBL-PBR pairs from DIT-HAP pipeline TSV files for downstream rendering.
Each input file is read with a three-level row index, filtered to strictly
positive PBL/PBR pairs, and written to a long-format TSV with sample/timepoint/
condition metadata derived from the filename.

Input
-----
- One or more TSV files (``-i``/``--input``), tab-separated, with a
  three-column MultiIndex (``index_col=[0, 1, 2]``) and ``PBL`` and ``PBR``
  data columns. Filenames must follow the pattern ``{sample}_{timepoint}_{condition}.tsv``.

Output
------
- A single long-format TSV file (``-o``/``--output``) with columns:
  ``sample``, ``timepoint``, ``condition``, ``chr``, ``coordinate``, ``strand``,
  ``pbl``, ``pbr``.

Usage
-----
    python PBL_PBR_correlation_analysis.py -i file1.tsv file2.tsv -o pairs.tsv
    python PBL_PBR_correlation_analysis.py -i file1.tsv -o pairs.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-13
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from qc.correlation import parse_filename, read_tsv_file  # noqa: E402

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PBLPBRCorrelationConfig:
    """Immutable config holding validated input TSV paths and the output TSV path."""
    input_files: list[Path]
    output_path: Path

    def __post_init__(self) -> None:
        """Validate that every input file exists and ensure the output directory is present."""
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Extract PBL-PBR pairs from multiple TSV files")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="Input TSV files (space-separated)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Read PBL/PBR TSVs, extract positive pairs with metadata, and write to long-format TSV."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Validate input and output paths via the config dataclass (raises ValueError on bad input)
    try:
        config = PBLPBRCorrelationConfig(
            input_files=args.input,
            output_path=args.output,
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    logger.info("=== PBL-PBR Pairs Extraction ===")
    logger.info(f"Processing {len(config.input_files)} input files...")

    # Sort files by name
    sorted_files = sorted(config.input_files, key=lambda x: x.name)
    logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")

    # Read and process files using vectorized pandas
    all_dataframes: list[pd.DataFrame] = []

    for file_path in sorted_files:
        filename = file_path.name
        logger.info(f"Reading {filename}...")

        # Parse filename to extract metadata
        metadata = parse_filename(file_path)
        if metadata is None:
            continue

        sample, timepoint, condition = metadata

        # Read and filter data
        df = read_tsv_file(file_path)
        if df is None:
            continue

        # Create long-format dataframe with vectorized operations
        pairs_df = df.rename(columns={
            'Chr': 'chr', 'Coordinate': 'coordinate', 'Strand': 'strand',
            'PBL': 'pbl', 'PBR': 'pbr',
        })
        pairs_df = pairs_df.assign(sample=sample, timepoint=timepoint, condition=condition)

        # Reorder columns to match spec: sample, timepoint, condition, chr, coordinate, strand, pbl, pbr
        pairs_df = pairs_df[['sample', 'timepoint', 'condition', 'chr', 'coordinate', 'strand', 'pbl', 'pbr']]

        all_dataframes.append(pairs_df)
        logger.info(f"  - Extracted {len(pairs_df)} pairs from {filename}")

    if not all_dataframes:
        logger.error("Error: No valid data found in any input file!")
        return 1

    # Concatenate all dataframes
    logger.info("Concatenating all pairs...")
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    # Write to TSV
    logger.info(f"Writing {len(combined_df)} pairs to {config.output_path}...")
    combined_df.to_csv(config.output_path, sep='\t', index=False, float_format='%.4f')

    logger.success(f"Extraction complete! Output saved to: {config.output_path}")
    logger.info(f"Total pairs extracted: {len(combined_df)}")
    logger.info(f"Files processed: {len(all_dataframes)}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
