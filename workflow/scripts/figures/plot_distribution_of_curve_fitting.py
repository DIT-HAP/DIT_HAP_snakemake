#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Distribution of Curve Fitting Results Figure Renderer
======================================================

Render histograms of curve fitting metrics (A, DR, DL, t50, R², RMSE, etc.)
from pre-computed fitting statistics TSV. One panel per metric in a grid layout.

Input
-----
- Fitting statistics TSV (``-i/--input``): curve fitting results with columns
  for model parameters and goodness-of-fit metrics.

Output
------
- ``<output>.pdf`` — journal-quality vector figure with histogram grid.
- ``<output>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_distribution_of_curve_fitting.py -i fitting_stats.tsv -o figures/fitting_distribution
    python plot_distribution_of_curve_fitting.py -i fitting_stats.tsv -o figures/fitting_distribution --bins 50 --verbose

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

import pandas as pd

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from figure_render.histogram import render_histogram_grid_figure  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV path and output stem."""
    fitting_stats_path: Path
    output_stem: Path
    bins: int = 30

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.fitting_stats_path.exists():
            raise ValueError(f"Fitting stats file does not exist: {self.fitting_stats_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)
        if self.bins < 1:
            raise ValueError(f"bins must be positive, got {self.bins}")


# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

SUCCESS_STATUS = "Success"

# Curve-fitting metrics worth histogramming, in display order. Columns absent
# from a given stats file are skipped.
METRIC_COLUMNS = [
    "A", "DR", "DL", "t10", "t50", "t90", "t_window", "t_inflection",
    "y_inflection", "auc", "R2", "RMSE", "normalized_RMSE", "AIC", "BIC",
]


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_fitting_stats(fitting_stats_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Load successful fits and return them with the metric columns actually present."""
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    df = pd.read_csv(fitting_stats_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)
    logger.info(f"Loaded {len(df)} rows")

    successful = df[df["Status"] == SUCCESS_STATUS].copy()
    logger.info(f"Found {len(successful)} successful fits")

    available = [column for column in METRIC_COLUMNS if column in successful.columns]
    logger.info(f"Found {len(available)} metric columns: {available}")

    return successful[available], available


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render distribution of curve fitting results")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input fitting statistics TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("--bins", type=int, default=30, help="Number of histogram bins (default: 30)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load data, render histograms, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            fitting_stats_path=args.input,
            output_stem=output_stem,
            bins=args.bins,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Distribution of Curve Fitting Results Rendering ===")

    try:
        df, metric_cols = load_fitting_stats(config.fitting_stats_path)

        render_histogram_grid_figure(
            df,
            config.output_stem,
            value_columns=metric_cols,
            bins=config.bins
        )

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
