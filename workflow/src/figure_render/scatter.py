"""Generic scatter and regression figure rendering.

Supersedes ``correlation.py`` (PBL vs PBR) and ``orientation.py`` (plus vs minus
strand), which rendered the same grouped log-log regression grid for different
columns. Column names and axis labels are supplied by the caller.

Callers pass explicitly log-transformed columns when they want log-space
statistics: r and P are computed directly from whatever columns are handed to
``render_scatter_panel``, so passing raw values reports a different statistic
while looking equally plausible. The same applies to ``ScatterPanel.density``:
a ``log_scale=True`` panel colours by linear-space density while displaying
symlog axes unless the columns were logged first.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import cnsplots as cns
import num2tex
import numpy as np
import pandas as pd
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.colors import Normalize
from scipy.stats import pearsonr

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import PANEL_DECORATION_PX, panel_labels, square_panel_size
from ._point_density import point_density
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

# Density panels bypass cns.scatterplot for ax.scatter (cns.scatterplot exposes
# only categorical hue and forces edgecolor=None, so a per-point colour array
# cannot be threaded through it). Higher alpha than
# REGRESSION_PANEL_SCATTER_KWS: at 0.15 the colour signal washes out, and
# overplotting is what the colour already encodes. No 'color'/'c' key -- those
# collide with the density array.
DENSITY_SCATTER_KWS: dict[str, Any] = {
    "s": 3,
    "alpha": 0.6,
    "linewidths": 0,
    "rasterized": True,
}

# Perceptually uniform, unlike cns.settings.palette_seq's 'gnuplot' default,
# so equal colour steps read as equal density steps.
DENSITY_CMAP = "viridis"

# Axes-fraction (x0, y0, width, height) for the inset colorbar: bottom-right,
# clear of the top-left n/r/P annotation. An inset axes contributes nothing to
# the panel width/height cnsplots measures for layout (Multipanel never walks
# ax.child_axes), so it cannot reflow the grid the way fig.colorbar(ax=...)
# would by taking space from the panel itself.
DENSITY_CBAR_BOUNDS = (0.88, 0.06, 0.035, 0.30)
DENSITY_CBAR_LABELSIZE = 5
DENSITY_CBAR_LABEL = "Density"


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
    density: bool = False


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
        rf"$n$={len(valid)}" "\n" rf"$r$={r:.2f}" "\n" rf"$P={num2tex.num2tex(p_value, precision=3):.2g}$",
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


def _add_density_colorbar(ax: Axes, mappable: PathCollection) -> None:
    """Inset a thin vertical colorbar labelled low/high inside the panel's bottom-right corner."""
    cax = ax.inset_axes(DENSITY_CBAR_BOUNDS)
    colorbar = ax.figure.colorbar(mappable, cax=cax, ticks=[0, 1])

    # Only the ordering is meaningful: density is renormalised per panel, so
    # numeric ticks would invite comparison between panels with different n.
    labels = colorbar.ax.set_yticklabels(["low", "high"])

    # End ticks anchor their label on the tick centre, pushing half the text
    # past the bar. Aligning each label inwards keeps both inside the panel.
    labels[0].set_verticalalignment("bottom")
    labels[-1].set_verticalalignment("top")

    # length=0 drops the tick marks but keeps their labels.
    colorbar.ax.tick_params(labelsize=DENSITY_CBAR_LABELSIZE, length=0, pad=1)
    colorbar.set_label(DENSITY_CBAR_LABEL, fontsize=DENSITY_CBAR_LABELSIZE, labelpad=2)
    colorbar.outline.set_linewidth(0.3)


def _draw_density_scatter(ax: Axes, df: pd.DataFrame, panel: ScatterPanel) -> bool:
    """Colour the panel's points by 2D density, returning False when the cloud is too degenerate.

    Deliberately ignores the caller's ``scatter_kws``: those are tuned for the
    cns.scatterplot path and carry a flat ``color`` that collides with the
    per-point colour array.
    """
    density = point_density(df[panel.x].to_numpy(), df[panel.y].to_numpy())

    if density is None:
        return False

    finite = np.isfinite(density)

    # Ascending so the crowded points draw last and are not buried under the
    # sparse ones.
    order = np.argsort(density[finite])
    x = df[panel.x].to_numpy()[finite][order]
    y = df[panel.y].to_numpy()[finite][order]

    # Rescaled to [0, 1] so the colorbar's low/high ticks sit at the bar ends.
    # A bare Normalize() autoscales to the raw density range instead, leaving
    # ticks at 0 and 1 out of range and clipped onto the edges.
    ranked = density[finite][order]
    span = np.ptp(ranked)
    scaled = (ranked - ranked.min()) / span if span > 0 else np.zeros_like(ranked)

    mappable = ax.scatter(
        x, y, c=scaled, cmap=DENSITY_CMAP, norm=Normalize(vmin=0, vmax=1), **DENSITY_SCATTER_KWS
    )
    _add_density_colorbar(ax, mappable)

    return True


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

    ``panel.density=True`` colours points by 2D kernel density with an inset
    colorbar. It is ignored when ``hue`` applies: a categorical hue and a
    continuous density both claim the marker colour.
    """
    ax.set_xlabel(panel.xlabel)
    ax.set_ylabel(panel.ylabel)
    ax.set_title(panel.title)

    if df.empty:
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        return

    kws = dict(SCATTERPLOT_KWS if scatter_kws is None else scatter_kws)

    if hue is not None and hue in df.columns:
        if panel.density:
            logger.warning(
                f"Panel '{panel.title}': hue {hue!r} and density both set; "
                f"colouring by hue and skipping density"
            )
        # Levels absent from the data must be dropped, or a project lacking
        # one of them raises inside seaborn.
        observed = set(df[hue].unique())
        order = [level for level in (hue_order or sorted(observed)) if level in observed]
        cns.scatterplot(data=df, x=panel.x, y=panel.y, hue=hue, hue_order=order, ax=ax, **kws)
        ax.legend(fontsize=6, loc="best", frameon=False)
    else:
        density_drawn = panel.density and _draw_density_scatter(ax, df, panel)
        if panel.density and not density_drawn:
            logger.warning(f"Panel '{panel.title}': falling back to flat colour, no density estimate")
        if not density_drawn:
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
    density: bool = False,
    scatter_kws: Mapping[str, object] | None = None,
    panel_decoration_px: int = PANEL_DECORATION_PX,
) -> None:
    """Render a row_key x col_key grid of square regression panels for the x/y columns.

    show_stats and density are set on every panel in the grid: this renderer
    builds its own ScatterPanel specs, so they cannot be varied per panel here.
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
            reference="identity", show_stats=show_stats, density=density,
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
