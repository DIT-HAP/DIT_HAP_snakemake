#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read Count Distribution Figure (single stage)
=============================================

One-stage QC figure: read per-(sample, condition) read-count tables, histogram
their log-spaced distributions, and mark the hard-filtering cutoff — all in
memory, with no intermediate binned TSV.

Each input file's numeric columns are one panel per timepoint; the cutoff
retention footer is computed from the raw initial-timepoint column via
``qc.read_counts.calculate_cutoff_statistics`` (same ``>=`` boundary as the
hard-filtering step). Rendering is delegated to
``figure_render.histogram.render_grouped_histogram_figure`` in log-scale mode:
bars are laid out evenly in log10 space exactly like the previous pre-binned
figure, while ticks show real read counts and the cutoff marker sits at the raw
cutoff value.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-27
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
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from sample_metadata import parse_sample_name  # noqa: E402
from qc.read_counts import STATS_COLUMNS, calculate_cutoff_statistics  # noqa: E402
from figure_render.histogram import render_grouped_histogram_figure  # noqa: E402


# =============================================================================
# CONSTANTS
# =============================================================================
X_LABEL = "Read count"
Y_LABEL = "Frequency"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV paths, output stem, and annotations."""
    input_files: list[Path]
    output_stem: Path
    initial_time_point: str
    cutoff: float
    bins: int = 30

    def __post_init__(self) -> None:
        """Validate that inputs exist, the cutoff is positive, and the output directory is present."""
        if not self.input_files:
            raise ValueError("At least one input file must be provided")
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        if self.cutoff <= 0:
            raise ValueError(f"Cutoff must be positive: {self.cutoff}")
        if self.bins < 5 or self.bins > 200:
            raise ValueError(f"Number of bins must be between 5 and 200: {self.bins}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_value_table(file_path: Path) -> pd.DataFrame | None:
    """Load one raw read-count table with its 4-level row index; None on parse failure."""
    logger.info(f"Reading TSV file: {file_path}")
    try:
        df = pd.read_csv(file_path, sep="\t", engine="python", index_col=[0, 1, 2, 3])
        if df.empty:
            raise ValueError("Empty DataFrame after reading")
    except (ValueError, KeyError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        logger.error(f"Failed to process {file_path.name}: {e}. Skipping.")
        return None
    return df


def format_retention_caption(sample: str, stats_row: dict[str, float | int | str]) -> str:
    """Format a one-line cutoff retention summary for a sample."""
    return (
        f"{sample}: {int(stats_row['rows_kept']):,}/{int(stats_row['original_rows']):,} rows kept "
        f"({stats_row['pct_rows_kept']:.1f}%), "
        f"{int(stats_row['counts_kept']):,}/{int(stats_row['original_counts']):,} counts kept "
        f"({stats_row['pct_counts_kept']:.1f}%)"
    )


@logger.catch(reraise=True)
def assemble_distribution(input_files: list[Path], initial_time_point: str, cutoff: float):
    """Build the long-format value frame and per-sample retention records from raw tables.

    Bad files are logged and skipped so one malformed input cannot kill the
    whole figure; raises when no file yields usable data.
    """
    value_frames: list[pd.DataFrame] = []
    retention_records: list[dict[str, float | int]] = []

    for file_path in sorted(input_files, key=lambda p: p.name):
        sample = parse_sample_name(file_path)
        df = load_value_table(file_path)
        if df is None:
            continue

        stats = calculate_cutoff_statistics(df, initial_time_point, cutoff)
        # @logger.catch swallows validation errors to None (e.g. a missing or
        # non-numeric initial timepoint column); treat that as a skipped file
        # rather than crashing on the ** below.
        if stats is None:
            continue

        numeric_cols = df.select_dtypes(include="number").columns
        for col in numeric_cols:
            value_frames.append(
                pd.DataFrame({"sample": [sample] * len(df), "timepoint": str(col), "value": df[col].to_numpy()})
            )
        retention_records.append({"sample": sample, **stats})

    if not value_frames:
        raise ValueError("No valid data found in any input file!")

    distribution_df = pd.concat(value_frames, ignore_index=True)
    retention_records.sort(key=lambda r: str(r["sample"]))
    log_lines = [
        format_retention_caption(r["sample"], r)
        for r in retention_records
    ]
    logger.info(f"Assembled {len(distribution_df)} values across {len(retention_records)} samples")
    return distribution_df, log_lines


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render read count distribution figure from raw read-count tables")
    parser.add_argument("-i", "--input", nargs="+", type=Path, required=True,
                        help="Input read-count TSV files ({sample}_{condition}.*.tsv)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-t", "--initial_time_point", type=str, required=True,
                        help="Initial time point column name for the cutoff annotation")
    parser.add_argument("-c", "--cutoff", type=float, required=True, help="Hard filtering cutoff value to annotate")
    parser.add_argument("--bins", type=int, default=50, help="Number of bins per panel (default: %(default)s)")
    parser.add_argument("-q", "--upper-quantile", type=float, default=None,
                        help="Drop values above this per-group quantile (e.g. 0.999999) to tame extreme outliers")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Compute retention in memory, histogram raw values, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            input_files=args.input,
            output_stem=output_stem,
            initial_time_point=args.initial_time_point,
            cutoff=args.cutoff,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Read Count Distribution Figure Rendering ===")

    try:
        df, retention_lines = assemble_distribution(config.input_files, config.initial_time_point, config.cutoff)

        render_grouped_histogram_figure(
            df,
            config.output_stem,
            value_column="value",
            row_key="sample",
            col_key="timepoint",
            bins=config.bins,
            log_scale=True,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            marker_value=config.cutoff,
            marker_label=f"Cutoff = {config.cutoff:.2g}",
            marker_on_col_value=config.initial_time_point,
            footer_lines=retention_lines,
            footer_header=f"Cutoff applied to '{config.initial_time_point}' (>= {config.cutoff:.2g}):",
            upper_quantile=args.upper_quantile
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
