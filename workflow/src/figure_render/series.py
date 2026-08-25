"""Generic multi-series scatter rendering.

Supersedes ``dispersions.py``, which plotted three named dispersion series
against a shared x column on log-log axes. Series columns, labels and colours
are supplied by the caller.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Rasterize: ~94k points per series would bloat the vector PDF. Alpha is kept low
# because overlapping clouds otherwise paint over each other -- at 0.5 the second
# series hides the first entirely.
SERIES_SCATTER_KWS: dict[str, object] = {
    "s": 1.0,
    "alpha": 0.12,
    "linewidths": 0,
    "rasterized": True,
}

# Legend handles inherit the scatter's tiny size and low alpha, so they are
# scaled up and forced opaque or the legend reads as blank.
LEGEND_MARKER_SCALE = 6


@dataclass(kw_only=True, slots=True, frozen=True)
class Series:
    """One y series drawn against the shared x column."""

    column: str
    label: str
    color: str


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_series_scatter_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    series: Sequence[Series],
    xlabel: str,
    ylabel: str,
    title: str,
    log_x: bool = True,
    log_y: bool = True,
) -> None:
    """Render every series against the shared x column on one axes."""
    logger.info(f"Rendering series scatter figure ({len(series)} series)...")

    require_columns(
        df,
        [x, *(item.column for item in series)],
        context="series scatter input",
    )

    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    ax = multipanel.panel()

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if df.empty:
        logger.warning("No data to plot!")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        save_dual(output_stem)
        return

    x_values = df[x]
    for item in series:
        logger.info(f"  Series {item.label}: {item.column} (n={df[item.column].notna().sum()})")
        ax.scatter(x_values, df[item.column], c=item.color, label=item.label, **SERIES_SCATTER_KWS)

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    legend = ax.legend(loc="best", frameon=False, markerscale=LEGEND_MARKER_SCALE)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
