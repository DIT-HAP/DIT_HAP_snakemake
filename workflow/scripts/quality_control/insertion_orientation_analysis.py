#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "loguru",
# ]
# ///

"""
Insertion Orientation Pairs Extraction
======================================

Extract paired ``+``/``-`` strand insertion counts per Sample x Timepoint from
one or more TSV files with multi-level row and column indexing, writing a
long-format TSV for downstream rendering. This is the computation half of the
insertion orientation figure; it performs no plotting.

Rows are stacked so each ``(Chr, Coordinate, Target, Timepoint, Sample)`` key
carries its ``+`` and ``-`` counts side by side, dropped where either strand is
missing, then filtered to rows whose minimum strand count is strictly positive.
That filtering matches the previous behaviour exactly and is what determines the
downstream correlation values.

Input
-----
- One or more TSV/TXT files (``-i``/``--input``) with a 4-level row index
  (``index_col=[0, 1, 2, 3]``) and a 2-level column header (``header=[0, 1]``)
  in which the ``Strand`` level holds the ``+`` / ``-`` orientation.

Output
------
- A single long-format TSV (``-o``/``--output``) with columns ``sample``,
  ``timepoint``, ``plus_count``, ``minus_count``.

Usage
-----
    python insertion_orientation_analysis.py -i file1.tsv file2.tsv -o strand_pairs.tsv
    python insertion_orientation_analysis.py -i file1.tsv -o strand_pairs.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  2.0.0
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

# 4. Local application imports (require the project ``src`` dir on sys.path)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_table  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
OUTPUT_COLUMNS = ["sample", "timepoint", "plus_count", "minus_count"]
READER_KWARGS = {"index_col": [0, 1, 2, 3], "header": [0, 1]}


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InsertionOrientationAnalysisConfig:
    """Validated input/output paths for the insertion orientation extraction."""

    input_files: list[Path]
    output_path: Path

    def __post_init__(self) -> None:
        """Validate that inputs exist with a TSV suffix and prepare the output directory."""
        if not self.input_files:
            raise ValueError("At least one input file must be provided")
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
            if file_path.suffix.lower() not in [".tsv", ".txt"]:
                raise ValueError(f"Input file must be a TSV file: {file_path}")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)




# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def extract_strand_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Pair +/- strand counts per Sample x Timepoint and keep rows with both strands positive."""
    required_levels = {"Strand", "Sample", "Timepoint"}
    available_levels = set(df.index.names) | set(df.columns.names)
    missing_levels = required_levels - available_levels
    if missing_levels:
        raise ValueError(f"Missing required index/column levels: {sorted(missing_levels)}")

    # Stack both column levels then pivot Strand back out, so each row holds the
    # +/- pair. dropna(axis=0) removes keys lacking one strand, exactly as before.
    plus_minus_pair = (
        df.stack(future_stack=True).stack(future_stack=True).unstack("Strand").dropna(axis=0)
    )
    logger.info(f"Paired {len(plus_minus_pair)} strand rows before positivity filtering")

    missing_strands = [strand for strand in ("+", "-") if strand not in plus_minus_pair.columns]
    if missing_strands:
        raise ValueError(f"Missing strand columns after unstacking: {missing_strands}")

    # Preserve the original min(axis=1) > 0 semantics: both strands strictly
    # positive, which is what makes the log-log correlation well defined.
    filtered = plus_minus_pair[plus_minus_pair.min(axis=1) > 0]
    logger.info(f"Retained {len(filtered)} rows with both strands strictly positive")

    if filtered.empty:
        logger.warning("No valid strand pairs after positivity filtering!")
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    pairs = (
        filtered[["+", "-"]]
        .rename(columns={"+": "plus_count", "-": "minus_count"})
        .reset_index()
        .rename(columns={"Sample": "sample", "Timepoint": "timepoint"})
    )

    return pairs[OUTPUT_COLUMNS]


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Extract insertion orientation (+/-) strand pairs from multiple TSV files.")
    parser.add_argument("-i", "--input", nargs="+", type=Path, required=True, help="One or more input TSV files with multi-level indexing.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV file path for the strand pairs.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Read strand-indexed TSVs, extract positive +/- pairs, and write a long-format TSV."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    try:
        config = InsertionOrientationAnalysisConfig(
            input_files=args.input,
            output_path=args.output,
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    logger.info("=== Insertion Orientation Pairs Extraction ===")
    logger.info(f"Processing {len(config.input_files)} input files...")

    sorted_files = sorted(config.input_files, key=lambda p: p.name)
    logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")

    all_pairs: list[pd.DataFrame] = []

    for file_path in sorted_files:
        logger.info(f"--- Processing file: {file_path.name} ---")

        # Per-file control flow: skip a file that fails to process, continue with the rest.
        try:
            df = read_table(file_path, **READER_KWARGS)
            pairs = extract_strand_pairs(df)
        except (ValueError, KeyError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error(f"Failed to process {file_path.name}: {e}")
            continue

        if pairs.empty:
            logger.warning(f"No valid strand pairs in {file_path.name}, skipping")
            continue

        all_pairs.append(pairs)
        logger.info(f"  - Extracted {len(pairs)} pairs from {file_path.name}")

    if not all_pairs:
        logger.error("Error: No valid data found in any input file!")
        return 1

    combined = pd.concat(all_pairs, ignore_index=True)

    logger.info(f"Writing {len(combined)} strand pairs to {config.output_path}...")
    combined.to_csv(config.output_path, sep="\t", index=False, float_format="%.4f")

    logger.success(f"Extraction complete! Output saved to: {config.output_path}")
    logger.info(f"Groups extracted: {combined.groupby(['sample', 'timepoint']).ngroups}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
