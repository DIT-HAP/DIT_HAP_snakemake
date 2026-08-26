"""Generic scatter and regression figure rendering.

Supersedes ``correlation.py`` (PBL vs PBR) and ``orientation.py`` (plus vs minus
strand), which rendered the same grouped log-log regression grid for different
columns. Column names and axis labels are supplied by the caller.

Callers pass explicitly log-transformed columns when they want log-space
statistics: cnsplots ``regplot`` computes and annotates r and P from whatever
columns it is handed, so passing raw values reports a different statistic while
looking equally plausible.

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
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import PANEL_DECORATION_PX, grid_panel_size, panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Verified to render ~50k points per panel legibly: density stays visible instead
# of collapsing into a solid block. Marker size must live inside these kwargs --
# regplot's own `s` argument is silently dropped when scatter_kws is supplied.
REGPLOT_SCATTER_KWS: dict[str, object] = {
    "s": 3,
    "facecolor": "none",
    "edgecolor": "gray",
    "alpha": 0.15,
    "linewidths": 0.25,
    "rasterized": True,
}

# Verified for ~5k points with categorical colour-coding. Deliberately NOT shared
# with REGPLOT_SCATTER_KWS: cns.scatterplot always passes edgecolor=None to
# seaborn internally, so supplying 'edgecolor' here raises on the duplicate.
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


# =============================================================================
# CORE LOGIC
# =============================================================================
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
    scatter_kws: Mapping[str, object] | None = None,
    panel_decoration_px: int = PANEL_DECORATION_PX,
) -> None:
    """Render a row_key x col_key grid of regression panels for the x/y columns."""
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
    n_rows = df[row_key].nunique()
    panel_width, panel_height = grid_panel_size(
        JOURNAL_WIDTH_PX, JOURNAL_HEIGHT_PX, n_cols, n_rows, decoration_px=panel_decoration_px
    )

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    labels = panel_labels(n_panels)
    kws = dict(REGPLOT_SCATTER_KWS if scatter_kws is None else scatter_kws)

    for panel_index, ((row_value, col_value), group_df) in enumerate(grouped):
        label = labels[panel_index]
        column = panel_index % max(n_cols, 1)

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
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

        def _apply_labels() -> None:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(panel_ylabel)
            ax.set_title(str(col_value))

        _apply_labels()

        if group_df.empty:
            logger.warning(f"  Panel {label}: {row_value} {col_value} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            continue

        logger.info(f"  Panel {label}: {row_value} {col_value} (n={len(group_df)})")

        cns.regplot(data=group_df, x=x, y=y, scatter_kws=kws, ax=ax)

        # regplot resets the axis labels it manages, so reapply after drawing.
        _apply_labels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


def _draw_reference_line(ax, df: pd.DataFrame, panel: ScatterPanel) -> None:
    """Draw the panel's reference line, sized to the data where the mode requires it."""
    match panel.reference:
        case "identity":
            max_val = max(df[panel.x].max(), df[panel.y].max(), 1)
            ax.plot([0, max_val], [0, max_val], color="red", linestyle="--", alpha=0.6, linewidth=1)
        case "unit_identity":
            ax.plot([0, 1], [0, 1], color="red", linestyle="--", alpha=0.6, linewidth=1)
        case "zero":
            ax.axhline(0, color="red", linestyle="--", alpha=0.6, linewidth=1)
        case "none":
            pass


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

    kws = dict(SCATTERPLOT_KWS if scatter_kws is None else scatter_kws)
    labels = panel_labels(len(panels))

    for label, panel in zip(labels, panels, strict=True):
        ax = multipanel.panel(label=label)

        if df.empty:
            logger.warning(f"  Panel {label}: {panel.title} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlabel(panel.xlabel)
            ax.set_ylabel(panel.ylabel)
            ax.set_title(panel.title)
            continue

        logger.info(f"  Panel {label}: {panel.title} (n={len(df)})")

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

        if panel.log_scale:
            ax.set_xscale("symlog")
            ax.set_yscale("symlog")

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
