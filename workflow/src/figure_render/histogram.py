"""Generic histogram grid rendering.

Supersedes ``distribution.py`` (one panel per metric, raw values) and
``read_counts.py`` (pre-binned counts replayed via weights, with a cutoff marker
and a figure-level footer). Both are grids of histograms, but their inputs
differ in form, so each has its own entry point here.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
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
PANEL_HEIGHT_PX = 180


# =============================================================================
# CORE LOGIC — single-axes primitive
# =============================================================================
def draw_histogram_panel(
    ax: Axes,
    data: pd.Series,
    *,
    bins: int,
    log_scale: bool = False,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    title: str = "",
    bin_edges: np.ndarray | None = None,
) -> None:
    """Histogram one Series onto ax; the shared primitive behind both grid renderers.

    Handles positive-value filtering in log mode, the empty-data placeholder,
    and the house bar style. When ``bin_edges`` is given it is used as-is (for
    callers coordinating ranges across panels); otherwise edges come from the
    data itself, laid evenly in log10 space under ``log_scale``.
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    values = data.dropna()
    if log_scale:
        values = values[values > 0]

    if values.empty:
        logger.warning("Panel has no valid data")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        return

    if bin_edges is not None:
        edges = bin_edges
    elif log_scale:
        edges = np.logspace(np.log10(values.min()), np.log10(values.max()), bins + 1)
    else:
        edges = bins

    ax.hist(values.to_numpy(), bins=edges, alpha=0.8, edgecolor="white", linewidth=0.5)  # type: ignore[arg-type]

    if log_scale:
        ax.set_xscale("log")


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
    share_y_range: bool = True,
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
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels}) panels)...")

    multipanel = cns.multipanel(max_width=PANEL_WIDTH_PX * n_cols)
    labels = panel_labels(n_panels)

    # Metric columns here are different physical quantities, so sharing an x
    # axis is meaningless; frequency magnitudes remain comparable, hence the
    # shared y by default. Columns that are entirely empty must be skipped so
    # a NaN max cannot poison the shared limit.
    finite_maxes = [m for c in value_columns if (m := df[c].max()) == m]
    shared_y_max = max(finite_maxes, default=0.0) if share_y_range else None

    for label, column in zip(labels, value_columns, strict=True):
        data = df[column].dropna()
        logger.info(f"  Panel {label}: {column} (n={len(data)})")

        ax = multipanel.panel(label=label, width=PANEL_WIDTH_PX, height=PANEL_HEIGHT_PX)

        had_data = data.notna().any()
        draw_histogram_panel(ax, data, bins=bins, xlabel=xlabel, ylabel=ylabel, title=column)
        if had_data and show_summary_stats:
            stats_text = f"n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}"
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment="top", fontsize=8)
        if shared_y_max is not None and had_data:
            ax.set_ylim(0, shared_y_max)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")



@logger.catch(reraise=True)
def render_grouped_histogram_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_column: str,
    row_key: str,
    col_key: str,
    bins: int = 50,
    log_scale: bool = False,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    marker_value: float | None = None,
    marker_label: str = "",
    marker_on_col_value: str | None = None,
    footer_lines: Sequence[str] = (),
    footer_header: str = "",
    upper_quantile: float | None = None,
    share_x_range: bool = True,
    share_y_range: bool = True,
) -> None:
    """Histogram one value column into one panel per row_key/col_key group.

    With ``log_scale`` the values must be strictly positive: bin edges are laid
    evenly in log10 space via ``np.logspace`` and the axis becomes a true log
    scale, so bars are geometrically identical to pre-transformed linear binning
    while ticks show real read counts. The cutoff marker is then placed at the
    raw value (e.g. x=8), not at its log.

    ``upper_quantile`` (e.g. 0.999999) drops values above each group's own
    quantile before binning — QC read-count tables carry one or two astronomic
    outliers that otherwise stretch the log axis across empty decades.

    The row_key value rides in the first column's ylabel rather than the title:
    a full "{row} {col}" title is wider than the axes and pushes the measured
    panel width past the point where n_cols panels fit in a row, which silently
    reflows the grid.
    """
    require_columns(df, [row_key, col_key, value_column], context="grouped histogram input")

    if df.empty:
        logger.warning("No values to plot!")
        return

    grouped = df.groupby([row_key, col_key], sort=True)
    if grouped.ngroups == 0:
        logger.warning("No groups to plot!")
        return

    # Shared ranges must be computed AFTER outlier dropping so a lone extreme
    # value does not stretch every panel's common x axis.
    def _panel_values(group_df: pd.DataFrame) -> pd.Series:
        data = group_df[value_column].dropna()
        if log_scale:
            data = data[data > 0]
        if upper_quantile is not None and not data.empty:
            threshold = data.quantile(upper_quantile)
            dropped = int((data > threshold).sum())
            if dropped:
                logger.info(f"upper_quantile={upper_quantile}: dropping {dropped} outlier(s) above {threshold:.4g}")
            data = data[data <= threshold]
        return data

    panel_values = {key: _panel_values(group_df) for key, group_df in grouped}
    if all(data.empty for data in panel_values.values()):
        logger.warning("All panels are empty after filtering!")
        return

    apply_house_style()

    n_cols = df[col_key].nunique()
    n_rows = df[row_key].nunique()
    panel_width, panel_height = grid_panel_size(
        JOURNAL_WIDTH_PX, JOURNAL_HEIGHT_PX, n_cols, n_rows,
        decoration_px=HISTOGRAM_DECORATION_PX,
    )

    footer_reserve_px = FOOTER_LINE_PX * (n_rows + 2) if footer_lines else 0
    last_row_index = n_rows - 1

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    labels = panel_labels(grouped.ngroups)

    if share_x_range and log_scale:
        lo = min(d.min() for d in panel_values.values() if not d.empty)
        hi = max(d.max() for d in panel_values.values() if not d.empty)
        shared_edges = np.logspace(np.log10(lo), np.log10(hi), bins + 1)
    else:
        shared_edges = None

    for panel_index, ((row_value, col_value), data) in enumerate(panel_values.items()):
        label = labels[panel_index]
        column = panel_index % max(n_cols, 1)
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

        panel_ylabel = f"{row_value}\n{ylabel}" if column == 0 else ylabel
        draw_histogram_panel(
            ax, data,
            bins=bins,
            log_scale=log_scale,
            xlabel=xlabel,
            ylabel=panel_ylabel,
            title=str(col_value),
            bin_edges=shared_edges,
        )

        if share_y_range and not data.empty:
            # One global y top over every panel; computing it per iteration is
            # wasteful but panels are few and histogramming is cheap next to
            # the render itself.
            y_max = max(
                int(np.histogram(d.to_numpy(), bins=(shared_edges if shared_edges is not None else bins))[0].max())
                for d in panel_values.values()
                if not d.empty
            )
            ax.set_ylim(0, y_max * 1.05)

        logger.info(f"  Panel {label}: {row_value} {col_value} (n={len(data):,})")

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
