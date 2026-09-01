"""Generic MA plot rendering.

Supersedes ``ma_plot.py`` (wide baseMean/LFC tables, one column per timepoint)
and ``ma_plot_replicates.py`` (one long table with padj-driven point colours).
Callers reshape either form into a sequence of (title, abundance, effect)
triples, so this module needs no knowledge of the input layout.

The two pipeline branches are selected by ``has_replicates`` and never run
together; both must keep working.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from enum import StrEnum
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

from figures import FURNITURE_COLOR, apply_house_style, panel_labels, save_dual


# =============================================================================
# CONSTANTS
# =============================================================================
class Orientation(StrEnum):
    """Panel arrangement and axis assignment for the MA plot."""

    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"


NONSIGNIFICANT_COLOR = FURNITURE_COLOR

# Rasterize: tens of thousands of points per panel would bloat the vector PDF.
MA_SCATTER_KWS: dict[str, object] = {
    "s": 1.5,
    "alpha": 0.4,
    "linewidths": 0,
    "rasterized": True,
}


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

    # Panels are identical in size and type, so a matplotlib grid applies: it
    # aligns their axes edges, which cns.multipanel cannot, and sharing is
    # declared up front instead of chained panel-to-panel afterwards.
    n_rows, n_cols = (n_panels, 1) if stacked else (1, n_panels)
    cns.figure(width=panel_width * n_cols, height=panel_height * n_rows)
    fig = plt.gcf()
    axes = fig.subplots(
        n_rows, n_cols, squeeze=False, sharex=share_axes, sharey=share_axes
    )

    labels = panel_labels(n_panels)
    colors = point_colors if point_colors is not None else [NONSIGNIFICANT_COLOR] * n_panels

    for index, (label, (title, abundance, effect), color) in enumerate(
        zip(labels, panels, colors, strict=True)
    ):
        logger.info(f"  Panel {label}: {title} (n={len(effect)})")

        ax = axes[index // n_cols][index % n_cols]

        match orientation:
            case Orientation.VERTICAL:
                # Effect on x, abundance (log) on y: reference line is vertical.
                ax.scatter(effect, abundance, c=color, **MA_SCATTER_KWS)
                ax.set_yscale("log")
                ax.axvline(null_effect, color=FURNITURE_COLOR, linestyle="--", zorder=3)
                ax.set_xlabel(effect_label)
                ax.set_ylabel(abundance_label)
            case Orientation.HORIZONTAL:
                # Abundance (log) on x, effect on y: reference line is horizontal.
                ax.scatter(abundance, effect, c=color, **MA_SCATTER_KWS)
                ax.set_xscale("log")
                ax.axhline(null_effect, color=FURNITURE_COLOR, linestyle="--", zorder=3)
                ax.set_xlabel(abundance_label)
                ax.set_ylabel(effect_label)

        ax.set_title(f"{title_prefix} - {title}")

        plt.sca(ax)
        cns.add_panel_label(label)

    fig.tight_layout()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
