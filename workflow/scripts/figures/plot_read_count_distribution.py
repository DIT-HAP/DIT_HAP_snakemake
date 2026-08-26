#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read Count Distribution Figure Renderer
=======================================

Render pre-binned log10 read-count histograms with the hard-filtering cutoff
marked on the initial-timepoint panel. Consumes the binned distribution TSV and
the cutoff statistics TSV produced by the computation layer; it never re-bins.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from loguru import logger
from matplotlib import use

use("Agg")

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figure_render.histogram import render_prebinned_histogram_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV paths, output stem, and annotations."""
    input_path: Path
    stats_path: Path
    output_stem: Path
    initial_time_point: str
    cutoff: float

    def __post_init__(self) -> None:
        """Validate that inputs exist, the cutoff is positive, and the output directory is present."""
        if not self.input_path.exists():
            raise ValueError(f"Input file does not exist: {self.input_path}")
        if not self.stats_path.exists():
            raise ValueError(f"Stats file does not exist: {self.stats_path}")
        if self.cutoff <= 0:
            raise ValueError(f"Cutoff must be positive: {self.cutoff}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CONSTANTS
# =============================================================================
DISTRIBUTION_COLUMNS = ["sample", "timepoint", "bin_left", "bin_right", "count"]
STATS_COLUMNS = [
    "sample", "original_rows", "rows_kept", "pct_rows_kept",
    "original_counts", "counts_kept", "pct_counts_kept",
]

X_LABEL = "log$_{10}$(read count)"
Y_LABEL = "Frequency"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_distribution_data(input_path: Path) -> pd.DataFrame:
    """Load the binned distribution TSV and validate its schema."""
    logger.info(f"Loading binned distribution from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, DISTRIBUTION_COLUMNS, context=f"distribution TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} bin rows")
    return df


@logger.catch(reraise=True)
def load_cutoff_stats(stats_path: Path) -> pd.DataFrame:
    """Load the cutoff retention statistics TSV and validate its schema."""
    logger.info(f"Loading cutoff statistics from {stats_path}...")
    df = pd.read_csv(stats_path, sep="\t")

    require_columns(df, STATS_COLUMNS, context=f"cutoff stats TSV {stats_path.name}")

    logger.info(f"Loaded statistics for {len(df)} samples")
    return df


def format_retention_caption(sample: str, stats_row: pd.Series) -> str:
    """Format a one-line cutoff retention summary for a sample."""
    return (
        f"{sample}: {int(stats_row['rows_kept']):,}/{int(stats_row['original_rows']):,} rows kept "
        f"({stats_row['pct_rows_kept']:.1f}%), "
        f"{int(stats_row['counts_kept']):,}/{int(stats_row['original_counts']):,} counts kept "
        f"({stats_row['pct_counts_kept']:.1f}%)"
    )


def build_retention_footer(df: pd.DataFrame, stats_df: pd.DataFrame) -> list[str]:
    """Build one retention line per sample present in the distribution data."""
    by_sample = stats_df.set_index("sample")
    present = set(df["sample"].unique())

    return [
        format_retention_caption(sample, row)
        for sample, row in by_sample.iterrows()
        if sample in present
    ]


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render read count distribution figure from binned TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input binned distribution TSV file")
    parser.add_argument("-s", "--stats", type=Path, required=True, help="Input cutoff statistics TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-t", "--initial_time_point", type=str, required=True, help="Initial time point column name for the cutoff annotation")
    parser.add_argument("-c", "--cutoff", type=float, required=True, help="Hard filtering cutoff value to annotate")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load binned distributions and stats, render the figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            input_path=args.input,
            stats_path=args.stats,
            output_stem=output_stem,
            initial_time_point=args.initial_time_point,
            cutoff=args.cutoff,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Read Count Distribution Figure Rendering ===")

    try:
        df = load_distribution_data(config.input_path)
        stats_df = load_cutoff_stats(config.stats_path)

        render_prebinned_histogram_figure(
            df,
            config.output_stem,
            row_key="sample",
            col_key="timepoint",
            left_column="bin_left",
            right_column="bin_right",
            count_column="count",
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            marker_value=float(np.log10(config.cutoff)),
            marker_label=f"Cutoff = {config.cutoff:.2g}",
            marker_on_col_value=config.initial_time_point,
            footer_lines=build_retention_footer(df, stats_df),
            footer_header=f"Cutoff applied to '{config.initial_time_point}' (>= {config.cutoff:.2g}):",
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
