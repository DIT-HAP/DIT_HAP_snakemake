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

import cnsplots as cns
import numpy as np
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
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru to emit uncolorised, timestamped records to stdout."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


# =============================================================================
# CORE LOGIC
# =============================================================================
def sigmoid_function(x: np.ndarray, A: float, DR: float, DL: float) -> np.ndarray:
    """Calculate sigmoid function values using gompertz function."""
    if A == 0:
        return np.zeros_like(x)
    alpha = (DR * np.e) / A
    u = alpha * (DL - x) + 1
    exponent = np.clip(u, -700, 700)
    return A * np.exp(-np.exp(exponent))


@logger.catch
def load_and_sample_data(config: PlotConfig) -> tuple[pd.DataFrame, pd.DataFrame, list[float]]:
    """Load fitting stats and LFC data, sample successful fits, return aligned data and time points."""
    logger.info(f"Loading fitting statistics from {config.fitting_stats_path}...")
    fitting_stats = pd.read_csv(config.fitting_stats_path, sep='\t', index_col=[0, 1, 2, 3])

    logger.info(f"Loaded {len(fitting_stats)} rows")

    # Filter to successful fits
    successful = fitting_stats[fitting_stats['Status'] == 'Success'].copy()
    logger.info(f"Found {len(successful)} successful fits")

    if successful.empty:
        logger.warning("No successful fits found!")
        return successful, pd.DataFrame(), []

    # Sample n_curves
    n_sample = min(config.n_curves, len(successful))
    sampled = successful.sample(n=n_sample, random_state=config.random_seed)
    logger.info(f"Sampled {len(sampled)} curves for plotting")

    # Load LFC data
    logger.info(f"Loading LFC data from {config.lfc_path}...")
    lfc_data = pd.read_csv(config.lfc_path, sep='\t', index_col=[0, 1, 2, 3])

    # Align with sampled indices
    lfc_sampled = lfc_data.loc[sampled.index]

    # Extract time points from time_points column
    time_points_str = sampled['time_points'].iloc[0]
    time_points = [float(t) for t in time_points_str.split(',')]
    logger.info(f"Time points: {time_points}")

    return sampled, lfc_sampled, time_points


@logger.catch
def render_curve_fitting_figure(sampled_stats: pd.DataFrame, lfc_sampled: pd.DataFrame,
                                  time_points: list[float], output_stem: Path) -> None:
    """Render multipanel figure with observed + fitted curves."""
    logger.info("Rendering curve fitting figure...")

    if sampled_stats.empty:
        logger.warning("No data to plot!")
        return

    # Apply house style before creating figure
    apply_house_style()

    n_curves = len(sampled_stats)
    n_cols = 4
    n_rows = (n_curves + n_cols - 1) // n_cols

    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_curves} panels)...")

    # Create multipanel layout
    panel_width = 180
    panel_height = 150
    multipanel = cns.multipanel(max_width=panel_width * n_cols)

    # Time points as numpy array
    x_values = np.array(time_points)
    timepoint_cols = [col for col in lfc_sampled.columns if col in ['YES0', 'YES1', 'YES2', 'YES3', 'YES4']]

    for idx, (gene_idx, row) in enumerate(sampled_stats.iterrows()):
        # Get observed values
        y_observed = lfc_sampled.loc[gene_idx, timepoint_cols].values

        # Get fitted parameters
        A, DR, DL = row['A'], row['DR'], row['DL']
        R2, RMSE = row['R2'], row['RMSE']

        # Format gene ID for title
        gene_id = "=".join(map(str, gene_idx))
        panel_label = chr(65 + idx) if idx < 26 else f"A{idx-25}"

        # Create panel
        ax = multipanel.panel(label=panel_label, width=panel_width, height=panel_height)

        # Plot observed points with ax.scatter directly
        ax.scatter(
            x_values,
            y_observed,
            s=30,
            color='#1f77b4',
            alpha=0.8,
            edgecolor='white',
            linewidths=0.5
        )

        # Plot fitted curve
        x_smooth = np.linspace(x_values.min(), x_values.max(), 100)
        y_fit = sigmoid_function(x_smooth, A, DR, DL)
        ax.plot(x_smooth, y_fit, color='#ff7f0e', linewidth=2, label='Fitted')

        # Add constraint lines
        ax.axhline(y=A, color='gray', linestyle='--', alpha=0.3, linewidth=1)
        ax.axvline(x=DL, color='gray', linestyle='--', alpha=0.3, linewidth=1)

        # Add parameter text
        param_text = f'A={A:.2f} R²={R2:.3f}\nDR={DR:.2f} RMSE={RMSE:.3f}\nDL={DL:.2f}'
        ax.text(0.05, 0.95, param_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=8)

        ax.set_title(gene_id.replace("=", " "), fontsize=9)
        ax.set_ylim(-1.5, 8.5)
        ax.set_xlabel('Time', fontsize=9)
        ax.set_ylabel('LFC', fontsize=9)
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
        sampled_stats, lfc_sampled, time_points = load_and_sample_data(config)

        # Render figure
        render_curve_fitting_figure(sampled_stats, lfc_sampled, time_points, config.output_stem)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
