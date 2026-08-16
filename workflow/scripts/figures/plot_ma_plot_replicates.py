#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Replicate-Branch MA Plot Figure Renderer
========================================

Render DESeq2 MA plots from a pre-computed TSV with one panel per timepoint:
log2 fold change against the mean of normalized counts, points coloured
darkred where the adjusted p-value clears the significance threshold and gray
otherwise. This is the rendering half of pydeseq2's ``DeseqStats.plot_MA()``.

This is the biological-replicates branch, distinct from ``plot_ma_plot.py``:
that script handles the no-replicate branch, which computes M/A directly from
median-normalized counts, whereas this one carries DESeq2's ``baseMean`` /
``log2FoldChange`` / ``padj``.

Rows with NaN ``padj`` (pydeseq2's independent filtering removes them from
testing) are kept and rendered gray, matching pydeseq2, where ``nan < alpha``
evaluates False.

Input
-----
- MA-values figure-data TSV written by ``write_ma_values_tsv()`` in
  ``insertion_level_depletion_analysis_has_replicates.py``: tab-separated, no
  index, columns ``timepoint``, ``baseMean``, ``log2FoldChange``, ``padj``;
  one row per insertion per non-initial timepoint.

Output
------
- ``<stem>.pdf`` — journal-quality vector figure (scatter rasterized).
- ``<stem>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_ma_plot_replicates.py -i ma_values_replicates.tsv -o figures/ma_plot
    python plot_ma_plot_replicates.py -i ma_values_replicates.tsv -o figures/ma_plot --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-16
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
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
REQUIRED_COLUMNS = ("timepoint", "baseMean", "log2FoldChange", "padj")

# pydeseq2's DeseqStats default alpha; the project config sets no override, so
# this is the threshold the legacy figures actually used.
PADJ_THRESHOLD = 0.05

SIGNIFICANT_COLOR = "darkred"
NONSIGNIFICANT_COLOR = "gray"

# pydeseq2's lfc_null default; no alt_hypothesis override, so a single line at 0
LFC_NULL = 0.0

X_LABEL = "mean of normalized counts"
Y_LABEL = "log2 fold change"


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
    """Load replicate-branch MA values TSV and validate schema."""
    logger.info(f"Loading MA values from {ma_values_path}...")
    df = pd.read_csv(ma_values_path, sep='\t')

    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")
    return df


def significance_colors(padj: pd.Series) -> pd.Series:
    """Map adjusted p-values to point colours; NaN padj stays non-significant like pydeseq2."""
    # (padj < threshold) is False for NaN, so filtered-out insertions render gray
    return (padj < PADJ_THRESHOLD).map({True: SIGNIFICANT_COLOR, False: NONSIGNIFICANT_COLOR})


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

    # Create multipanel layout (vertical stack). The default 10 px bottom margin
    # is too tight for a stack: each panel's title would collide with the x-axis
    # label of the panel above it.
    panel_width = 510
    panel_height = 200
    panel_margin_bottom = 32
    multipanel = cns.multipanel(max_width=panel_width)

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_panels)]

    # Rasterize: ~94k points per panel would bloat the vector PDF
    scatter_kws = dict(
        s=1.5,
        alpha=0.4,
        linewidths=0,
        rasterized=True
    )

    for timepoint, label in zip(timepoints, panel_labels):
        group_df = df[df['timepoint'] == timepoint]
        logger.info(f"  Panel {label}: {timepoint} (n={len(group_df)})")

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            margin_bottom=panel_margin_bottom,
        )

        if group_df.empty:
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
            ax.set_xlabel(X_LABEL)
            ax.set_ylabel(Y_LABEL)
            ax.set_title(f'MA plot - {timepoint}')
            continue

        colors = significance_colors(group_df['padj'])
        n_significant = int((colors == SIGNIFICANT_COLOR).sum())
        logger.debug(f"    padj < {PADJ_THRESHOLD}: {n_significant} ({n_significant / len(group_df):.1%})")

        # Scatter log2FoldChange vs baseMean
        ax.scatter(
            group_df['baseMean'],
            group_df['log2FoldChange'],
            c=colors,
            **scatter_kws
        )

        # Log x only: baseMean spans several orders of magnitude, LFC is symmetric
        ax.set_xscale('log')

        # Add reference line at the null log fold change
        ax.axhline(LFC_NULL, color='red', alpha=0.5, linestyle='--', linewidth=1, zorder=3)

        ax.set_xlabel(X_LABEL)
        ax.set_ylabel(Y_LABEL)
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
    parser = argparse.ArgumentParser(description="Render replicate-branch MA plot figure from pre-computed MA values TSV")
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

    logger.info("=== Replicate-Branch MA Plot Figure Rendering ===")

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
