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

import cnsplots as cns
import numpy as np
import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual  # noqa: E402


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
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load pairs TSV, validate schema, filter positive values, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep='\t')

    # Validate required columns
    required_cols = ['sample', 'timepoint', 'condition', 'pbl', 'pbr']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")

    # Filter to positive values
    df_filtered = df[(df['pbl'] > 0) & (df['pbr'] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(df_filtered)} rows")

    if df_filtered.empty:
        logger.warning("No valid data points after filtering!")
        return df_filtered

    # Add log10 columns
    df_filtered['log10_pbl'] = np.log10(df_filtered['pbl'])
    df_filtered['log10_pbr'] = np.log10(df_filtered['pbr'])

    return df_filtered


@logger.catch
def render_correlation_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render multipanel correlation figure with one panel per sample-timepoint group."""
    logger.info("Rendering correlation figure...")

    # Apply house style before creating figure
    apply_house_style()

    # Group by sample and timepoint
    grouped = df.groupby(['sample', 'timepoint'], sort=True)
    n_panels = len(grouped)

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    # Create figure
    fig = cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_panels)]

    # Verified scatter kwargs for 131k points
    scatter_kws = dict(
        s=3,
        facecolor="none",
        edgecolor="gray",
        alpha=0.15,
        linewidths=0.25,
        rasterized=True
    )

    for (sample, timepoint), group_df in grouped:
        label = panel_labels.pop(0)
        logger.info(f"  Panel {label}: {sample} {timepoint} (n={len(group_df)})")

        if group_df.empty:
            # Handle empty data case
            ax = multipanel.panel(label=label)
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
            ax.set_xlabel('log$_{10}$ PBL')
            ax.set_ylabel('log$_{10}$ PBR')
            ax.set_title(f'{sample} {timepoint}')
            continue

        # Create panel and regplot
        ax = multipanel.panel(label=label)

        # regplot on log10 columns
        cns.regplot(
            data=group_df,
            x='log10_pbl',
            y='log10_pbr',
            scatter_kws=scatter_kws,
            ax=ax
        )

        # Set labels with proper subscript and title
        ax.set_xlabel('log$_{10}$ PBL')
        ax.set_ylabel('log$_{10}$ PBR')
        ax.set_title(f'{sample} {timepoint}')

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


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
        # Load and prepare data
        df = load_and_prepare_data(config.input_path)

        # Render figure
        render_correlation_figure(df, config.output_stem)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
