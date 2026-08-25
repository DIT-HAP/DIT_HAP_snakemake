#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Insertion Orientation Figure Renderer
=====================================

Render log-log correlation plots of ``+`` versus ``-`` strand insertion counts
from the pre-computed strand pairs TSV. One panel per sample and timepoint, with
regression line and correlation statistics via cnsplots.regplot.

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

import numpy as np
import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from figure_render.scatter import render_grouped_regression_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV path and output stem."""
    input_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.input_path.exists():
            raise ValueError(f"Input file does not exist: {self.input_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["sample", "timepoint", "plus_count", "minus_count"]

X_COLUMN = "log10_plus_count"
Y_COLUMN = "log10_minus_count"
X_LABEL = "log$_{10}$ (+) strand"
Y_LABEL = "log$_{10}$ (-) strand"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load the strand pairs TSV, keep strictly-positive pairs, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"strand pairs TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} rows")

    # The computation layer already applies min > 0, but re-filtering keeps the
    # renderer safe against a hand-made input and guarantees finite log10 values.
    positive = df[(df["plus_count"] > 0) & (df["minus_count"] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(positive)} rows")

    if positive.empty:
        logger.warning("No valid data points after filtering!")
        return positive

    positive[X_COLUMN] = np.log10(positive["plus_count"])
    positive[Y_COLUMN] = np.log10(positive["minus_count"])

    return positive


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render insertion orientation figure from strand pairs TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input strand pairs TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load strand pairs, render the orientation figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            input_path=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Insertion Orientation Figure Rendering ===")

    try:
        df = load_and_prepare_data(config.input_path)

        render_grouped_regression_figure(
            df,
            config.output_stem,
            x=X_COLUMN,
            y=Y_COLUMN,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            row_key="sample",
            col_key="timepoint",
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
