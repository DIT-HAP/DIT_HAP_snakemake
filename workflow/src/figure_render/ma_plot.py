"""MA plot figure rendering (unified, no-replicate branch).

Input
-----
- baseMean TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex, tab-separated,
  one column per timepoint.
- LFC TSV: same row MultiIndex and timepoint columns as the baseMean table.

Output
------
- `<stem>.pdf` / `<stem>.review.png` — vertical stack, one panel per timepoint,
  LFC on the x-axis and baseMean (log-scaled) on the y-axis.
- `<stem>_horizontal.pdf` / `<stem>_horizontal.review.png` — the same panels
  arranged in a single row, with baseMean (log-scaled) on the x-axis and LFC
  on the y-axis.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from enum import StrEnum
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

# =============================================================================
# CONSTANTS
# =============================================================================
class Orientation(StrEnum):
    """Panel arrangement and axis assignment for the MA plot."""

    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"

NONSIGNIFICANT_COLOR = "gray"
LFC_NULL = 0.0
LFC_LABEL = "log2 fold change"
BASEMEAN_LABEL = "mean of normalized counts"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_ma_data(basemean_path: Path, lfc_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the wide-format baseMean/LFC tables and validate their columns match."""
    logger.info(f"Loading baseMean from {basemean_path}...")
    basemean_df = pd.read_csv(basemean_path, sep="\t", index_col=[0, 1, 2, 3])

    logger.info(f"Loading LFC from {lfc_path}...")
    lfc_df = pd.read_csv(lfc_path, sep="\t", index_col=[0, 1, 2, 3])

    missing_cols = [col for col in lfc_df.columns if col not in basemean_df.columns]
    if missing_cols:
        raise ValueError(f"Timepoints present in LFC but missing from baseMean: {missing_cols}")

    logger.info(f"Loaded {len(lfc_df)} insertions across {len(lfc_df.columns)} timepoints")

    return basemean_df, lfc_df


@logger.catch
def render_ma_figure(
    basemean_df: pd.DataFrame,
    lfc_df: pd.DataFrame,
    output_stem: Path,
    orientation: Orientation = Orientation.VERTICAL,
) -> None:
    """Render MA plot with one panel per timepoint column, stacked vertically or in a single row."""
    logger.info(f"Rendering MA plot figure ({orientation})...")

    if lfc_df.empty:
        logger.warning("No data to plot!")
        return

    # Apply house style before creating figure
    apply_house_style()

    timepoints = list(lfc_df.columns)
    n_panels = len(timepoints)

    logger.info(f"Creating figure with {n_panels} panels (one per timepoint)...")

    # Square-ish panels, sharing both axes so panels are directly comparable.
    panel_size = 220
    match orientation:
        case Orientation.VERTICAL:
            # A small max_width (one panel wide) forces one panel per row. The
            # default 10 px bottom margin is too tight for a stack: each panel's
            # title would collide with the x-axis label of the panel above it.
            multipanel = cns.multipanel(max_width=panel_size)
            panel_margin_bottom = 32
        case Orientation.HORIZONTAL:
            # A max_width spanning all panels forces a single row instead of a stack.
            # Each panel's true rendered width exceeds `panel_size` (log-scale y-axis
            # tick labels add a measured ~30-40 px "left_reserve", plus the 10 px
            # margin_right); layout wraps to a new row once the running width sum
            # exceeds max_width, so budget a generous +80 px/panel headroom. An
            # oversized max_width is safe: the figure is always rendered at exactly
            # max_width regardless of how much of it the panels actually fill.
            multipanel = cns.multipanel(max_width=panel_size * n_panels + 80 * n_panels)
            panel_margin_bottom = 10

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_panels)]

    # Rasterize: tens of thousands of points per panel would bloat the vector PDF
    scatter_kws = dict(
        s=1.5,
        alpha=0.4,
        linewidths=0,
        rasterized=True,
    )

    first_ax = None
    for timepoint, label in zip(timepoints, panel_labels):
        logger.info(f"  Panel {label}: {timepoint} (n={len(lfc_df[timepoint])})")

        ax = multipanel.panel(
            label=label,
            width=panel_size,
            height=panel_size,
            margin_bottom=panel_margin_bottom,
        )

        if first_ax is None:
            first_ax = ax
        else:
            ax.sharex(first_ax)
            ax.sharey(first_ax)

        match orientation:
            case Orientation.VERTICAL:
                # LFC on x, baseMean (log) on y: reference line is vertical.
                ax.scatter(lfc_df[timepoint], basemean_df[timepoint], c=NONSIGNIFICANT_COLOR, **scatter_kws)
                ax.set_yscale("log")
                ax.axvline(LFC_NULL, color="red", alpha=0.5, linestyle="--", linewidth=1, zorder=3)
                ax.set_xlabel(LFC_LABEL)
                ax.set_ylabel(BASEMEAN_LABEL)
            case Orientation.HORIZONTAL:
                # baseMean (log) on x, LFC on y: reference line is horizontal.
                ax.scatter(basemean_df[timepoint], lfc_df[timepoint], c=NONSIGNIFICANT_COLOR, **scatter_kws)
                ax.set_xscale("log")
                ax.axhline(LFC_NULL, color="red", alpha=0.5, linestyle="--", linewidth=1, zorder=3)
                ax.set_xlabel(BASEMEAN_LABEL)
                ax.set_ylabel(LFC_LABEL)

        ax.set_title(f"MA plot - {timepoint}")

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
