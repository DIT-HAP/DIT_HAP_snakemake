#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Distribution of Curve Fitting Results Figure Renderer
======================================================

Render histograms of curve fitting metrics (A, DR, DL, t50, R², RMSE, etc.)
from pre-computed fitting statistics TSV. One panel per metric in a grid layout.

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

import cnsplots as cns
import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import apply_house_style, save_dual  # noqa: E402


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
# CORE LOGIC
# =============================================================================
@logger.catch
def load_fitting_stats(fitting_stats_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Load fitting statistics and identify numeric columns for histogramming."""
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    df = pd.read_csv(fitting_stats_path, sep='\t', index_col=[0, 1, 2, 3])

    logger.info(f"Loaded {len(df)} rows")

    # Filter to successful fits
    successful = df[df['Status'] == 'Success'].copy()
    logger.info(f"Found {len(successful)} successful fits")

    # Select numeric columns of interest (exclude fitted/residual columns)
    metric_cols = ['A', 'DR', 'DL', 't10', 't50', 't90', 't_window', 't_inflection', 'y_inflection', 'auc', 'R2', 'RMSE', 'normalized_RMSE', 'AIC', 'BIC']
    available_cols = [col for col in metric_cols if col in successful.columns]

    logger.info(f"Found {len(available_cols)} metric columns: {available_cols}")

    return successful[available_cols], available_cols


@logger.catch
def render_distribution_figure(df: pd.DataFrame, metric_cols: list[str], output_stem: Path, bins: int) -> None:
    """Render multipanel histogram figure with one panel per metric."""
    logger.info("Rendering distribution figure...")

    if df.empty or not metric_cols:
        logger.warning("No data to plot!")
        return

    # Apply house style before creating figure
    apply_house_style()

    n_metrics = len(metric_cols)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols

    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_metrics} panels)...")

    # Create multipanel layout
    panel_width = 180
    panel_height = 150
    multipanel = cns.multipanel(max_width=panel_width * n_cols)

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_metrics)]

    for metric, label in zip(metric_cols, panel_labels):
        data = df[metric].dropna()

        logger.info(f"  Panel {label}: {metric} (n={len(data)})")

        # Create panel
        ax = multipanel.panel(label=label, width=panel_width, height=panel_height)

        if data.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(metric)
            continue

        # Histogram - use ax.hist directly since cns.histplot requires DataFrame
        ax.hist(
            data,
            bins=bins,
            alpha=0.8,
            edgecolor='white',
            linewidth=0.5
        )

        # Add statistics text
        stats_text = f'n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=8)

        ax.set_title(metric, fontsize=10)
        ax.set_xlabel('Value', fontsize=9)
        ax.set_ylabel('Frequency', fontsize=9)
        ax.tick_params(labelsize=8)

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


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
        # Load data
        df, metric_cols = load_fitting_stats(config.fitting_stats_path)

        # Render figure
        render_distribution_figure(df, metric_cols, config.output_stem, config.bins)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
