#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PBL-PBR Correlation Figure Renderer
===================================

Render log-log correlation plots of PBL vs PBR from pre-computed pairs TSV.
Groups data by sample and timepoint, rendering one panel per group with
regression line and correlation statistics via cnsplots.regplot.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-13
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
# CONSTANTS
# =============================================================================
# PBL/PBR-specific knowledge lives here, not in the shared renderer.
REQUIRED_COLUMNS = ["sample", "timepoint", "pbl", "pbr"]

VALUE_COLUMNS = ["pbl", "pbr"]
X_COLUMN = "log10_pbl"
Y_COLUMN = "log10_pbr"
X_LABEL = "log$_{10}$ PBL"
Y_LABEL = "log$_{10}$ PBR"


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
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load the pairs TSV, keep strictly-positive pairs, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"pairs TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} rows")

    positive = df[(df["pbl"] > 0) & (df["pbr"] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(positive)} rows")

    if positive.empty:
        logger.warning("No valid data points after filtering!")
        return positive

    # regplot annotates r and P from the columns it is handed, so the log10
    # columns must be explicit: passing raw values would report raw-space
    # correlation (r ~ -0.03) instead of the log-space PCC (r = 0.85) that this
    # figure has always shown.
    positive[X_COLUMN] = np.log10(positive["pbl"])
    positive[Y_COLUMN] = np.log10(positive["pbr"])

    return positive


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render PBL-PBR correlation figure from pairs TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input pairs TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load pairs data, render correlation figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            input_path=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== PBL-PBR Correlation Figure Rendering ===")

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
