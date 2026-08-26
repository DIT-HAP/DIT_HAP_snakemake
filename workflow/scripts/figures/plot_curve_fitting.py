#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Curve Fitting Figure Renderer
==============================

Render observed vs fitted sigmoid curves for a random sample of genes.
Reads fitting statistics TSV and original LFC data, samples N curves
(default 32, deterministic seed 42), and draws observed points + fitted
sigmoid overlay in a grid layout.

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

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from depletion.curve_model import sigmoid_function  # noqa: E402
from figure_render.curves import render_fitted_curves_figure  # noqa: E402

import pandas as pd  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated paths and sampling parameters."""
    fitting_stats_path: Path
    lfc_path: Path
    output_stem: Path
    n_curves: int = 32
    random_seed: int = 42

    def __post_init__(self) -> None:
        """Validate that inputs exist and output directory can be created."""
        if not self.fitting_stats_path.exists():
            raise ValueError(f"Fitting stats file does not exist: {self.fitting_stats_path}")
        if not self.lfc_path.exists():
            raise ValueError(f"LFC file does not exist: {self.lfc_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)
        if self.n_curves < 1:
            raise ValueError(f"n_curves must be positive, got {self.n_curves}")


# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

# The sigmoid's parameters, in the order sigmoid_function() accepts them.
MODEL_PARAM_COLUMNS = ["A", "DR", "DL"]

# Printed in each panel's text box.
ANNOTATION_COLUMNS = ["R2", "RMSE"]

SUCCESS_STATUS = "Success"
X_LABEL = "Time"
Y_LABEL = "LFC"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_and_sample_data(
    fitting_stats_path: Path, lfc_path: Path, n_curves: int, random_seed: int
) -> tuple[pd.DataFrame, list[float], list[str]]:
    """Sample successful fits, join their observed LFC values, and return the x values.

    Timepoint columns are taken from the LFC table, not a name whitelist: the old
    renderer tested `col in ['YES0'..'YES4']`, which selected the wrong subset for
    any project with a different timepoint count or naming scheme.
    """
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    fitting_stats = pd.read_csv(fitting_stats_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)
    logger.info(f"Loaded {len(fitting_stats)} rows")

    successful = fitting_stats[fitting_stats["Status"] == SUCCESS_STATUS].copy()
    logger.info(f"Found {len(successful)} successful fits")

    if successful.empty:
        logger.warning("No successful fits found!")
        return successful, [], []

    sampled = successful.sample(n=min(n_curves, len(successful)), random_state=random_seed)
    logger.info(f"Sampled {len(sampled)} curves for plotting")

    logger.info(f"Loading LFC data from {lfc_path}...")
    lfc_data = pd.read_csv(lfc_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)

    timepoint_columns = list(lfc_data.columns)
    time_points = [float(t) for t in sampled["time_points"].iloc[0].split(",")]
    logger.info(f"Time points: {time_points} for columns {timepoint_columns}")

    # The stats table's own timepoint columns collide with the LFC ones, so the
    # LFC values are joined under a suffix and then renamed back.
    lfc_sampled = lfc_data.loc[sampled.index, timepoint_columns]
    joined = sampled.drop(columns=timepoint_columns, errors="ignore").join(lfc_sampled)

    return joined, time_points, timepoint_columns


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render curve fitting figure from stats and LFC data")
    parser.add_argument("-s", "--stats", type=Path, required=True, help="Input fitting statistics TSV file")
    parser.add_argument("-l", "--lfc", type=Path, required=True, help="Input LFC TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-n", "--n-curves", type=int, default=32, help="Number of curves to sample (default: 32)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling (default: 42)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load data, sample curves, render figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            fitting_stats_path=args.stats,
            lfc_path=args.lfc,
            output_stem=output_stem,
            n_curves=args.n_curves,
            random_seed=args.seed,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Curve Fitting Figure Rendering ===")

    try:
        joined, time_points, timepoint_columns = load_and_sample_data(
            config.fitting_stats_path, config.lfc_path, config.n_curves, config.random_seed
        )

        render_fitted_curves_figure(
            joined,
            config.output_stem,
            x_values=time_points,
            value_columns=timepoint_columns,
            model=sigmoid_function,
            model_params=MODEL_PARAM_COLUMNS,
            annotations=ANNOTATION_COLUMNS,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
        )

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
