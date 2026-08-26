"""Generic histogram grid rendering.

Supersedes ``distribution.py`` (one panel per metric, raw values) and
``read_counts.py`` (pre-binned counts replayed via weights, with a cutoff marker
and a figure-level footer). Both are grids of histograms, but their inputs
differ in form, so each has its own entry point here.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import grid_panel_size, panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Height reserved per line of the footer at the bottom of the figure.
FOOTER_LINE_PX = 9

# Frequency histograms carry wider decorations than log10 scatter panels, whose
# tick labels are single-digit.
HISTOGRAM_DECORATION_PX = 55

PANEL_WIDTH_PX = 180
PANEL_HEIGHT_PX = 150


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_histogram_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_columns: Sequence[str],
    bins: int,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    show_summary_stats: bool = True,
    n_cols: int = 4,
) -> None:
    """Render one histogram panel per named value column."""
    logger.info("Rendering histogram grid figure...")

    require_columns(df, value_columns, context="histogram input")

    if df.empty or not value_columns:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    n_panels = len(value_columns)
    n_rows = (n_panels + n_cols - 1) // n_cols
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels} panels)...")

    multipanel = cns.multipanel(max_width=PANEL_WIDTH_PX * n_cols)
    labels = panel_labels(n_panels)

    for label, column in zip(labels, value_columns, strict=True):
        data = df[column].dropna()
        logger.info(f"  Panel {label}: {column} (n={len(data)})")

        ax = multipanel.panel(label=label, width=PANEL_WIDTH_PX, height=PANEL_HEIGHT_PX)

        if data.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(column)
            continue

        # ax.hist rather than cns.histplot: histplot requires a DataFrame, and
        # this path already holds a plain Series.
        ax.hist(data, bins=bins, alpha=0.8, edgecolor="white", linewidth=0.5)

        if show_summary_stats:
            stats_text = f"n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}"
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment="top", fontsize=8)

        ax.set_title(column, fontsize=10)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


@logger.catch(reraise=True)
def render_prebinned_histogram_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    row_key: str,
    col_key: str,
    left_column: str,
    right_column: str,
    count_column: str,
    xlabel: str,
    ylabel: str,
    marker_value: float | None = None,
    marker_label: str = "",
    marker_on_col_value: str | None = None,
    footer_lines: Sequence[str] = (),
    footer_header: str = "",
) -> None:
    """Replay pre-computed bins as one histogram panel per row_key/col_key group."""
    logger.info("Rendering pre-binned histogram figure...")

    require_columns(
        df,
        [row_key, col_key, left_column, right_column, count_column],
        context="pre-binned histogram input",
    )

    if df.empty:
        logger.warning("No bins to plot!")
        return

    grouped = df.groupby([row_key, col_key], sort=True)
    if grouped.ngroups == 0:
        logger.warning("No groups to plot!")
        return

    apply_house_style()

    n_cols = df[col_key].nunique()
    n_rows = df[row_key].nunique()
    panel_width, panel_height = grid_panel_size(
        JOURNAL_WIDTH_PX, JOURNAL_HEIGHT_PX, n_cols, n_rows,
        decoration_px=HISTOGRAM_DECORATION_PX,
    )

    # multipanel resizes the figure to fit its panels, so the footer strip must be
    # reserved inside the layout via the last row's bottom margin. Shrinking the
    # requested figure height instead lands the footer on the last row's labels.
    footer_reserve_px = FOOTER_LINE_PX * (n_rows + 2) if footer_lines else 0
    last_row_index = n_rows - 1

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    labels = panel_labels(grouped.ngroups)

    for panel_index, ((row_value, col_value), group_df) in enumerate(grouped):
        label = labels[panel_index]
        is_last_row = (panel_index // max(n_cols, 1)) == last_row_index

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            pad_left=2,
            pad_top=2,
            margin_right=4,
            margin_bottom=20 + (footer_reserve_px if is_last_row else 0),
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{row_value} {col_value}")

        valid = group_df.dropna(subset=[left_column, right_column])
        total = float(valid[count_column].sum()) if not valid.empty else 0.0

        if valid.empty or total <= 0:
            logger.warning(f"  Panel {label}: {row_value} {col_value} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            continue

        logger.info(f"  Panel {label}: {row_value} {col_value} ({len(valid)} bins, n={int(total):,})")

        # Bin centres weighted by their counts, with bin count and range from the
        # stored edges, so nothing is re-binned. An explicit edge array cannot be
        # used: seaborn compares `bins == "auto"` before binning, which raises on
        # an array when weights are supplied.
        centers = (valid[left_column] + valid[right_column]) / 2.0
        plot_df = valid.assign(_bin_center=centers)
        binrange = (float(valid[left_column].min()), float(valid[right_column].max()))

        cns.histplot(
            data=plot_df,
            x="_bin_center",
            weights=count_column,
            bins=len(valid),
            binrange=binrange,
            ax=ax,
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{row_value} {col_value}")

        if marker_value is None:
            continue
        if marker_on_col_value is not None and col_value != marker_on_col_value:
            continue

        ax.axvline(marker_value, color="firebrick", linestyle="--", linewidth=0.8, label=marker_label)
        ax.legend(frameon=False, loc="upper right", fontsize=4)

    if footer_lines:
        footer = "\n".join([footer_header, *footer_lines]) if footer_header else "\n".join(footer_lines)
        plt.gcf().text(0.02, 0.004, footer, ha="left", va="bottom", fontsize=5, linespacing=1.4)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
