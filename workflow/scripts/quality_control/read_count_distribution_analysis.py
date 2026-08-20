#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy",
#     "pandas",
#     "loguru",
# ]
# ///

"""
Read Count Distribution Computation
===================================

Compute log10 read-count histograms and hard-filtering retention statistics from
one or more TSV files, writing both as intermediate TSVs for downstream
rendering. This is the computation half of the read count distribution figure;
it performs no plotting.

Histograms are pre-binned here rather than emitting raw long-format values: a
long format would multiply the row count by the number of time points, turning a
10 MB input into hundreds of MB. The renderer must therefore consume these bins
as-is and never re-bin.

Bins are computed over ``log10(value)`` using only strictly positive values, and
the cutoff keeps rows whose initial-time-point value is greater than or equal to
the cutoff, matching the previous behaviour exactly.

Input
-----
- One or more TSV files (``-i``/``--input``) with a 4-level row index
  (``index_col=[0, 1, 2, 3]``) and one numeric read-count column per time point.

Output
------
- ``-o``/``--output``: binned distribution TSV with columns ``sample``,
  ``timepoint``, ``bin_left``, ``bin_right``, ``count``. A group with no
  strictly positive values contributes a single marker row with empty
  ``bin_left``/``bin_right`` and ``count`` of 0, so the renderer can still lay
  out a panel for it.
- ``-s``/``--stats``: cutoff retention TSV with columns ``sample``,
  ``original_rows``, ``rows_kept``, ``pct_rows_kept``, ``original_counts``,
  ``counts_kept``, ``pct_counts_kept``.

Usage
-----
    python read_count_distribution_analysis.py -i sample1.tsv sample2.tsv -t YES0 -c 8 -o dist.tsv -s stats.tsv
    python read_count_distribution_analysis.py -i *.tsv -t YES0 -c 8 -o dist.tsv -s stats.tsv --bins 80 --verbose

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
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
DISTRIBUTION_COLUMNS = ["sample", "timepoint", "bin_left", "bin_right", "count"]
STATS_COLUMNS = [
    "sample",
    "original_rows",
    "rows_kept",
    "pct_rows_kept",
    "original_counts",
    "counts_kept",
    "pct_counts_kept",
]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ReadCountDistributionAnalysisConfig:
    """Immutable, validated configuration for read count distribution computation."""

    input_files: list[Path]
    output_path: Path
    stats_path: Path
    initial_time_point: str
    cutoff: float
    bins: int = 50

    def __post_init__(self) -> None:
        """Validate configuration values and create the output directories."""
        if not self.input_files:
            raise ValueError("At least one input file must be provided")
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
            if file_path.suffix.lower() not in [".tsv", ".txt"]:
                raise ValueError(f"Input file must be a TSV file: {file_path}")
        if self.cutoff <= 0:
            raise ValueError("Cutoff value must be positive")
        if self.bins < 5 or self.bins > 200:
            raise ValueError("Number of bins must be between 5 and 200")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.stats_path.parent.mkdir(parents=True, exist_ok=True)




# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_and_validate_data(file_path: Path) -> pd.DataFrame:
    """Load a TSV file with its 4-level row index and reject an empty table."""
    df = pd.read_csv(file_path, sep="\t", engine="python", index_col=[0, 1, 2, 3])

    if df.empty:
        raise ValueError("Empty DataFrame after reading")

    return df


@logger.catch
def parse_sample_name(file_path: Path) -> str:
    """Derive the sample label from a filename by dropping every dotted suffix."""
    return file_path.name.split(".")[0]


@logger.catch
def calculate_cutoff_statistics(df: pd.DataFrame, initial_time_point: str, cutoff: float) -> dict[str, float | int]:
    """Compute row and count retention after applying the cutoff to the initial time point."""
    if initial_time_point not in df.columns:
        raise ValueError(f"Initial time point column '{initial_time_point}' not found in {list(df.columns)}")

    if not pd.api.types.is_numeric_dtype(df[initial_time_point]):
        raise ValueError(f"Initial time point column '{initial_time_point}' is not numeric")

    original_rows = len(df)
    original_counts = float(df[initial_time_point].sum())

    kept = df.loc[df[initial_time_point] >= cutoff, initial_time_point]
    rows_kept = len(kept)
    counts_kept = float(kept.sum())

    return {
        "original_rows": original_rows,
        "rows_kept": rows_kept,
        "pct_rows_kept": (rows_kept / original_rows) * 100.0 if original_rows > 0 else 0.0,
        "original_counts": original_counts,
        "counts_kept": counts_kept,
        "pct_counts_kept": (counts_kept / original_counts) * 100.0 if original_counts > 0 else 0.0,
    }


@logger.catch
def compute_binned_distribution(df: pd.DataFrame, sample: str, bins: int) -> pd.DataFrame:
    """Bin log10 of strictly positive values for every numeric column into a long-format frame."""
    numeric_cols = df.select_dtypes(include=np.number).columns
    if numeric_cols.empty:
        raise ValueError("No numeric columns found")

    frames: list[pd.DataFrame] = []

    for col_name in numeric_cols:
        col_data = df[col_name].dropna()
        positive = col_data[col_data > 0]

        # Marker row keeps the group present in the output so the renderer can
        # still allocate a "No valid data" panel for it.
        if positive.empty:
            logger.warning(f"  {sample} / {col_name}: no positive values, emitting empty marker row")
            frames.append(
                pd.DataFrame(
                    {
                        "sample": [sample],
                        "timepoint": [str(col_name)],
                        "bin_left": [np.nan],
                        "bin_right": [np.nan],
                        "count": [0],
                    }
                )
            )
            continue

        counts, edges = np.histogram(np.log10(positive.to_numpy()), bins=bins)

        frames.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "timepoint": str(col_name),
                    "bin_left": edges[:-1],
                    "bin_right": edges[1:],
                    "count": counts,
                }
            )
        )
        logger.info(f"  {sample} / {col_name}: {len(positive)} positive values across {bins} bins")

    return pd.concat(frames, ignore_index=True)[DISTRIBUTION_COLUMNS]


@logger.catch
def log_summary_table(stats_df: pd.DataFrame) -> None:
    """Log the per-sample retention statistics as a plain-text table."""
    if stats_df.empty:
        logger.info("No statistics to display.")
        return

    logger.info("--- Processing Summary ---")
    for line in stats_df.to_string(index=False).split("\n"):
        logger.info(line)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Compute binned read count distributions and cutoff statistics from TSV files.")
    parser.add_argument("-i", "--input", nargs="+", type=Path, required=True, help="One or more input TSV files.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV path for the binned distribution.")
    parser.add_argument("-s", "--stats", type=Path, required=True, help="Output TSV path for the cutoff retention statistics.")
    parser.add_argument("-t", "--initial_time_point", required=True, type=str, help="Name of the column representing the initial time point for cutoff application.")
    parser.add_argument("-c", "--cutoff", required=True, type=float, help="Cutoff value to apply to the initial time point column (values >= cutoff are kept).")
    parser.add_argument("--bins", type=int, default=50, help="Number of bins for histograms (default: %(default)s).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Bin read count distributions, compute cutoff retention, and write both TSVs."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = ReadCountDistributionAnalysisConfig(
            input_files=args.input,
            output_path=args.output,
            stats_path=args.stats,
            initial_time_point=args.initial_time_point,
            cutoff=args.cutoff,
            bins=args.bins,
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    logger.info("=== Read Count Distribution Computation ===")
    logger.info(f"Processing {len(config.input_files)} input files")
    logger.info(f"Initial time point column: '{config.initial_time_point}'")
    logger.info(f"Cutoff value: {config.cutoff}")
    logger.info(f"Histogram bins: {config.bins}")

    distribution_frames: list[pd.DataFrame] = []
    stats_records: list[dict[str, float | int | str]] = []

    for file_path in sorted(config.input_files, key=lambda p: p.name):
        sample = parse_sample_name(file_path)
        logger.info(f"--- Processing file: {file_path.name} (sample: {sample}) ---")

        # Per-file control flow: skip a file that fails to parse, continue with the rest.
        try:
            df = load_and_validate_data(file_path)
            stats = calculate_cutoff_statistics(df, config.initial_time_point, config.cutoff)
            distribution_frames.append(compute_binned_distribution(df, sample, config.bins))
        except (ValueError, KeyError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error(f"Failed to process {file_path.name}: {e}. Skipping.")
            continue

        stats_records.append({"sample": sample, **stats})

    if not distribution_frames:
        logger.error("Error: No valid data found in any input file!")
        return 1

    distribution_df = pd.concat(distribution_frames, ignore_index=True)
    stats_df = pd.DataFrame(stats_records)[STATS_COLUMNS]

    logger.info(f"Writing {len(distribution_df)} bin rows to {config.output_path}...")
    distribution_df.to_csv(config.output_path, sep="\t", index=False, float_format="%.6f")

    logger.info(f"Writing {len(stats_df)} stats rows to {config.stats_path}...")
    stats_df.to_csv(config.stats_path, sep="\t", index=False, float_format="%.4f")

    log_summary_table(stats_df)
    logger.success(f"Computation complete! Outputs: {config.output_path}, {config.stats_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
