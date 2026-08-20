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
from figure_render.curve_fitting import load_and_sample_data, render_curve_fitting_figure  # noqa: E402


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
        # Load and sample data
        sampled_stats, lfc_sampled, time_points = load_and_sample_data(
            config.fitting_stats_path, config.lfc_path, config.n_curves, config.random_seed
        )

        # Render figure
        render_curve_fitting_figure(sampled_stats, lfc_sampled, time_points, config.output_stem)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
