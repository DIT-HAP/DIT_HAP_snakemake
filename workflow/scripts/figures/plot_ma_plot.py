#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MA Plot Figure Renderer
=======================

Render MA plots (M vs A) from pre-computed TSV with one panel per timepoint.
M = log-fold change, A = average expression.

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
    ma_values_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.ma_values_path.exists():
            raise ValueError(f"MA values file does not exist: {self.ma_values_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


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
@logger.catch
def load_ma_data(ma_values_path: Path) -> pd.DataFrame:
    """Load MA values TSV and validate schema."""
    logger.info(f"Loading MA values from {ma_values_path}...")
    df = pd.read_csv(ma_values_path, sep='\t')

    required_cols = ['timepoint', 'M', 'A']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")
    return df


@logger.catch
def render_ma_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render MA plot with one panel per timepoint."""
    logger.info("Rendering MA plot figure...")

    if df.empty:
        logger.warning("No data to plot!")
        return

    # Apply house style before creating figure
    apply_house_style()

    # Group by timepoint
    timepoints = sorted(df['timepoint'].unique())
    n_panels = len(timepoints)

    logger.info(f"Creating figure with {n_panels} panels (one per timepoint)...")

    # Create multipanel layout (vertical stack)
    panel_width = 510
    panel_height = 200
    multipanel = cns.multipanel(max_width=panel_width)

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_panels)]

    # Scatter kwargs for many points
    scatter_kws = dict(
        s=10,
        facecolor="none",
        edgecolor="black",
        alpha=0.5,
        linewidths=0.5,
        rasterized=True
    )

    for timepoint, label in zip(timepoints, panel_labels):
        group_df = df[df['timepoint'] == timepoint]
        logger.info(f"  Panel {label}: {timepoint} (n={len(group_df)})")

        if group_df.empty:
            ax = multipanel.panel(label=label, width=panel_width, height=panel_height)
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
            ax.set_xlabel('M value')
            ax.set_ylabel('A value')
            ax.set_title(f'MA plot - {timepoint}')
            continue

        # Create panel
        ax = multipanel.panel(label=label, width=panel_width, height=panel_height)

        # Scatter plot M vs A
        ax.scatter(
            group_df['M'],
            group_df['A'],
            **scatter_kws
        )

        # Add reference line at M=0
        ax.axvline(0, color='red', linestyle='--', linewidth=2, alpha=0.5)

        ax.set_xlabel('M value')
        ax.set_ylabel('A value')
        ax.set_title(f'MA plot - {timepoint}')

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render MA plot figure from pre-computed MA values TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input MA values TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load MA data, render figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            ma_values_path=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== MA Plot Figure Rendering ===")

    try:
        # Load data
        df = load_ma_data(config.ma_values_path)

        # Render figure
        render_ma_figure(df, config.output_stem)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
