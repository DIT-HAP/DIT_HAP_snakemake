"""Generic scatter and regression figure rendering.

Supersedes ``correlation.py`` (PBL vs PBR) and ``orientation.py`` (plus vs minus
strand), which rendered the same grouped log-log regression grid for different
columns. Column names and axis labels are supplied by the caller.

Callers pass explicitly log-transformed columns when they want log-space
statistics: r and P are computed directly from whatever columns are handed to
``render_scatter_panel``, so passing raw values reports a different statistic
while looking equally plausible.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import cnsplots as cns
import num2tex
import numpy as np
import pandas as pd
from loguru import logger
from matplotlib.axes import Axes
from scipy.stats import pearsonr

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import PANEL_DECORATION_PX, panel_labels, square_panel_size
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Verified to render ~50k points per panel legibly: density stays visible instead
# of collapsing into a solid block. cns.scatterplot always forces edgecolor=None
# internally, so a filled gray marker (not regplot's hollow-circle look) is the
# closest equivalent that doesn't collide with that forced kwarg.
REGRESSION_PANEL_SCATTER_KWS: dict[str, object] = {
    "s": 3,
    "color": "gray",
    "alpha": 0.15,
    "rasterized": True,
}

# Verified for ~5k points with categorical colour-coding. Deliberately NOT shared
# with REGRESSION_PANEL_SCATTER_KWS: cns.scatterplot always passes edgecolor=None
# to seaborn internally, so supplying 'edgecolor' here raises on the duplicate.
SCATTERPLOT_KWS: dict[str, object] = {
    "s": 6,
    "alpha": 0.4,
    "rasterized": True,
}


@dataclass(kw_only=True, slots=True, frozen=True)
class ScatterPanel:
    """One panel of a scatter grid: which columns, how to label, which reference line."""

    x: str
    y: str
    xlabel: str
    ylabel: str
    title: str
    reference: Literal["identity", "zero", "unit_identity", "none"] = "none"
    log_scale: bool = False
    show_stats: bool = False


# =============================================================================
# CORE LOGIC
# =============================================================================
def _annotate_fit_stats(ax: Axes, df: pd.DataFrame, *, x: str, y: str) -> None:
    """Annotate n, r, and P at regplot's own text position/style, so toggling show_stats doesn't reflow the panel."""
    valid = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) < 2:
        logger.warning(f"Not enough finite pairs ({len(valid)}) to compute r/P for {x!r} vs {y!r}")
        return

    r, p_value = pearsonr(valid[x], valid[y])
    ax.text(
        0.05,
        0.95,
        rf"$n$={len(valid)}" "\n" rf"$r$={r:.2f}" "\n" rf"$P={num2tex.num2tex(p_value, precision=2):.2g}$",
        color="black",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )


def _draw_reference_line(ax: Axes, df: pd.DataFrame, panel: ScatterPanel) -> None:
    """Draw the panel's reference guide line, dashed gray to read as a guide rather than a fit.

    Fitted/regression lines (drawn elsewhere via cns.regplot or curves.py) stay
    solid so the two kinds of line remain visually distinct.
    """
    match panel.reference:
        case "identity":
            min_val = min(df[panel.x].min(), df[panel.y].min())
            max_val = max(df[panel.x].max(), df[panel.y].max())
            ax.plot([min_val, max_val], [min_val, max_val], color="gray", linestyle="--", alpha=0.6, linewidth=1)
        case "unit_identity":
            ax.plot([0, 1], [0, 1], color="gray", linestyle="--", alpha=0.6, linewidth=1)
        case "zero":
            ax.axhline(0, color="gray", linestyle="--", alpha=0.6, linewidth=1)
        case "none":
            pass


def render_scatter_panel(
    ax: Axes,
    df: pd.DataFrame,
    panel: ScatterPanel,
    *,
    hue: str | None = None,
    hue_order: Sequence[str] | None = None,
    scatter_kws: Mapping[str, object] | None = None,
) -> None:
    """Draw one ScatterPanel (or an empty-data placeholder) onto ax.

    The single-axes primitive behind both ``render_scatter_grid_figure`` and
    ``render_grouped_regression_figure``. Any other grid orchestrator that
    wants one cns.scatterplot panel with this project's hue/reference-line/
    log-scale/stats handling can call this directly instead of duplicating
    that sequence.

    ``panel.show_stats=True`` annotates n, r, and P (computed over all of df,
    ignoring hue) in the same text position/style cnsplots' own regplot used.
    """
    ax.set_xlabel(panel.xlabel)
    ax.set_ylabel(panel.ylabel)
    ax.set_title(panel.title)

    if df.empty:
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        return

    kws = dict(SCATTERPLOT_KWS if scatter_kws is None else scatter_kws)

    if hue is not None and hue in df.columns:
        # Levels absent from the data must be dropped, or a project lacking
        # one of them raises inside seaborn.
        observed = set(df[hue].unique())
        order = [level for level in (hue_order or sorted(observed)) if level in observed]
        cns.scatterplot(data=df, x=panel.x, y=panel.y, hue=hue, hue_order=order, ax=ax, **kws)
        ax.legend(fontsize=6, loc="best", frameon=False)
    else:
        cns.scatterplot(data=df, x=panel.x, y=panel.y, ax=ax, **kws)

    ax.set_xlabel(panel.xlabel)
    ax.set_ylabel(panel.ylabel)
    ax.set_title(panel.title)

    _draw_reference_line(ax, df, panel)

    if panel.show_stats:
        _annotate_fit_stats(ax, df, x=panel.x, y=panel.y)

    if panel.log_scale:
        ax.set_xscale("symlog")
        ax.set_yscale("symlog")


@logger.catch(reraise=True)
def render_grouped_regression_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    row_key: str = "sample",
    col_key: str = "timepoint",
    show_stats: bool = True,
    scatter_kws: Mapping[str, object] | None = None,
    panel_decoration_px: int = PANEL_DECORATION_PX,
) -> None:
    """Render a row_key x col_key grid of square regression panels for the x/y columns.

    show_stats controls the per-panel n/r/P annotation, forwarded unchanged to
    ``render_scatter_panel`` for every panel in the grid.
    """
    logger.info("Rendering grouped regression figure...")

    # Validated even when df is empty, so a typo'd column name is caught rather
    # than silently returning with no artifact and no error.
    require_columns(df, [row_key, col_key, x, y], context="scatter input")

    if df.empty:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    grouped = df.groupby([row_key, col_key], sort=True)
    n_panels = grouped.ngroups

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    # One row per row_key value and one column per col_key value, so a row stays
    # on one line instead of wrapping mid-group.
    n_cols = df[col_key].nunique()
    # A single edge from the width budget alone: cnsplots grows the figure's
    # total height to fit however many rows the grid needs regardless of what
    # height is requested here, so budgeting height from n_rows independently
    # (as grid_panel_size does) produces a rectangle whenever n_cols != n_rows.
    panel_size = square_panel_size(JOURNAL_WIDTH_PX, n_cols, decoration_px=panel_decoration_px)

    # multipanel wraps to a new row once a row's rendered panels (including
    # cnsplots' own measured y-axis tick/label reserve, ~30-40px per panel
    # beyond the requested axes width) exceed max_width. Budgeting max_width
    # from the requested panel_size alone therefore wraps one column early:
    # n_cols panels at panel_size each can measure wider than n_cols *
    # panel_size once decoration is added back in. Multiplying back in the
    # same panel_decoration_px keeps a full row within max_width regardless of
    # how much square_panel_size's min-size floor grew panel_size beyond the
    # width budget's natural per-column share.
    row_max_width = n_cols * (panel_size + panel_decoration_px)

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=row_max_width)

    labels = panel_labels(n_panels)

    for panel_index, ((row_value, col_value), group_df) in enumerate(grouped):
        label = labels[panel_index]
        column = panel_index % max(n_cols, 1)

        ax = multipanel.panel(
            label=label,
            width=panel_size,
            height=panel_size,
            pad_left=2,
            pad_top=2,
            margin_right=6,
            margin_bottom=18,
        )

        # The row_key value rides in the first column's ylabel rather than the
        # title: a full "{row} {col}" title is wider than the axes and pushes the
        # measured panel width past the point where n_cols panels fit in a row,
        # which silently reflows the grid.
        panel_ylabel = f"{row_value}\n{ylabel}" if column == 0 else ylabel

        if group_df.empty:
            logger.warning(f"  Panel {label}: {row_value} {col_value} has no valid data")
        else:
            logger.info(f"  Panel {label}: {row_value} {col_value} (n={len(group_df)})")

        panel = ScatterPanel(
            x=x, y=y, xlabel=xlabel, ylabel=panel_ylabel, title=str(col_value),
            reference="identity", show_stats=show_stats,
        )
        render_scatter_panel(
            ax, group_df, panel,
            scatter_kws=REGRESSION_PANEL_SCATTER_KWS if scatter_kws is None else scatter_kws,
        )

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


@logger.catch(reraise=True)
def render_scatter_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    panels: Sequence[ScatterPanel],
    hue: str | None = None,
    hue_order: Sequence[str] | None = None,
    scatter_kws: Mapping[str, object] | None = None,
) -> None:
    """Render one scatter panel per ScatterPanel spec, with optional hue colouring."""
    logger.info(f"Rendering scatter grid with {len(panels)} panels...")

    if not panels:
        logger.warning("No panels requested!")
        return

    # Validated even when df is empty, so a typo'd column name is caught rather
    # than silently rendering a placeholder panel.
    for panel in panels:
        require_columns(df, [panel.x, panel.y], context=f"scatter panel '{panel.title}'")

    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    labels = panel_labels(len(panels))

    for label, panel in zip(labels, panels, strict=True):
        ax = multipanel.panel(label=label)

        if df.empty:
            logger.warning(f"  Panel {label}: {panel.title} has no valid data")
        else:
            logger.info(f"  Panel {label}: {panel.title} (n={len(df)})")

        render_scatter_panel(ax, df, panel, hue=hue, hue_order=hue_order, scatter_kws=scatter_kws)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
