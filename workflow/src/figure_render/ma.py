"""Generic MA plot rendering.

Supersedes ``ma_plot.py`` (wide baseMean/LFC tables, one column per timepoint)
and ``ma_plot_replicates.py`` (one long table with padj-driven point colours).
Callers reshape either form into a sequence of (title, abundance, effect)
triples, so this module needs no knowledge of the input layout.

The two pipeline branches are selected by ``has_replicates`` and never run
together; both must keep working.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from enum import StrEnum
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

from ._layout import panel_labels

# =============================================================================
# CONSTANTS
# =============================================================================
class Orientation(StrEnum):
    """Panel arrangement and axis assignment for the MA plot."""

    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"


NONSIGNIFICANT_COLOR = "gray"

# Rasterize: tens of thousands of points per panel would bloat the vector PDF.
MA_SCATTER_KWS: dict[str, object] = {
    "s": 1.5,
    "alpha": 0.4,
    "linewidths": 0,
    "rasterized": True,
}

# A stack needs more bottom margin than the 10 px default: each panel's title
# would otherwise collide with the x-axis label of the panel above it.
STACK_MARGIN_BOTTOM = 32
ROW_MARGIN_BOTTOM = 10

# Each panel's rendered width exceeds `panel_width`: log-scale tick labels add a
# measured ~30-40 px left reserve plus the 10 px margin_right. Layout wraps once
# the running width sum exceeds max_width, so a single row needs generous
# headroom per panel. Oversizing max_width is safe -- the figure is rendered at
# exactly max_width regardless of how much the panels fill.
ROW_WIDTH_HEADROOM_PX = 80


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_ma_figure(
    panels: Sequence[tuple[str, pd.Series, pd.Series]],
    output_stem: Path,
    *,
    abundance_label: str,
    effect_label: str,
    title_prefix: str,
    orientation: Orientation = Orientation.VERTICAL,
    stack: bool | None = None,
    point_colors: Sequence[pd.Series | str] | None = None,
    panel_width: int = 220,
    panel_height: int = 220,
    share_axes: bool = True,
    null_effect: float = 0.0,
) -> None:
    """Render one MA panel per (title, abundance, effect) triple."""
    logger.info(f"Rendering MA plot figure ({orientation})...")

    if not panels:
        logger.warning("No data to plot!")
        return

    if point_colors is not None and len(point_colors) != len(panels):
        raise ValueError(
            f"point_colors must have one entry per panel: "
            f"got {len(point_colors)} colours for {len(panels)} panels"
        )

    apply_house_style()

    n_panels = len(panels)
    logger.info(f"Creating figure with {n_panels} panels...")

    # Arrangement is independent of axis assignment: the replicate branch draws
    # horizontal axes in a vertical stack. Defaulting stack from orientation keeps
    # the no-replicate branch's calls unchanged.
    stacked = (orientation is Orientation.VERTICAL) if stack is None else stack

    if stacked:
        # A one-panel-wide max_width forces one panel per row. The default 10 px
        # bottom margin is too tight for a stack: each panel's title would collide
        # with the x-axis label of the panel above it.
        multipanel = cns.multipanel(max_width=panel_width)
        margin_bottom = STACK_MARGIN_BOTTOM
    else:
        multipanel = cns.multipanel(
            max_width=(panel_width + ROW_WIDTH_HEADROOM_PX) * n_panels
        )
        margin_bottom = ROW_MARGIN_BOTTOM

    labels = panel_labels(n_panels)
    colors = point_colors if point_colors is not None else [NONSIGNIFICANT_COLOR] * n_panels

    first_ax = None
    for label, (title, abundance, effect), color in zip(labels, panels, colors, strict=True):
        logger.info(f"  Panel {label}: {title} (n={len(effect)})")

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            margin_bottom=margin_bottom,
        )

        if share_axes:
            if first_ax is None:
                first_ax = ax
            else:
                ax.sharex(first_ax)
                ax.sharey(first_ax)

        match orientation:
            case Orientation.VERTICAL:
                # Effect on x, abundance (log) on y: reference line is vertical.
                ax.scatter(effect, abundance, c=color, **MA_SCATTER_KWS)
                ax.set_yscale("log")
                ax.axvline(null_effect, color="red", alpha=0.5, linestyle="--", linewidth=1, zorder=3)
                ax.set_xlabel(effect_label)
                ax.set_ylabel(abundance_label)
            case Orientation.HORIZONTAL:
                # Abundance (log) on x, effect on y: reference line is horizontal.
                ax.scatter(abundance, effect, c=color, **MA_SCATTER_KWS)
                ax.set_xscale("log")
                ax.axhline(null_effect, color="red", alpha=0.5, linestyle="--", linewidth=1, zorder=3)
                ax.set_xlabel(abundance_label)
                ax.set_ylabel(effect_label)

        ax.set_title(f"{title_prefix} - {title}")

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
