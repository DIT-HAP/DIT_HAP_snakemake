"""Replicate-branch MA plot figure rendering.

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

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

# =============================================================================
# CONSTANTS
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
