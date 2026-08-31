"""Shared panel layout helpers for figure rendering.

Panel-letter generation and grid sizing were duplicated across five figure
modules in two incompatible variants; this module is the single implementation.

Two layout strategies live here. ``square_panel_size`` / ``grid_panel_size``
feed ``cns.multipanel``, which flows panels left-to-right and derives each
panel's origin from the *rendered* width of that panel's own y-axis
decorations, so panels in the same grid column start at different x whenever
their tick labels differ in width. ``aligned_grid_axes`` instead places axes on
a fixed pixel pitch, which is what a QC grid comparing the same quantity across
samples needs.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-30
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence

import cnsplots as cns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

# =============================================================================
# CONSTANTS
# =============================================================================
# Pixels each panel spends on its label, tick labels and axis title, on top of
# the axes box itself. Subtracted from the available width/height so a whole row
# of panels fits inside the journal width. Measured at 40 px for log10 scatter
# panels; frequency histograms carry wider decorations and pass 55 explicitly.
# Overestimating shrinks the axes enough for an extra panel to fit in a row,
# which silently reflows the grid.
PANEL_DECORATION_PX = 40

# Per-cell reserves for aligned_grid_axes, in cns.figure() layout pixels (the
# 72-dpi base cnsplots sizes figures in, not rendered device pixels). Every cell
# gets the same reserve — that uniformity is what makes columns and rows line
# up — so each must be wide/tall enough for the *worst* panel in the grid:
#   left:   panel label + a two-line rotated ylabel + y tick labels
#   top:    axes title + panel label
#   bottom: x tick labels + xlabel
# Undersizing does not shift any axes (the pitch is fixed); it only lets text
# spill past the figure edge, which savefig's tight bbox then absorbs.
ALIGNED_GRID_LEFT_PX = 44
ALIGNED_GRID_RIGHT_PX = 8
ALIGNED_GRID_TOP_PX = 20
ALIGNED_GRID_BOTTOM_PX = 26

# Offset from the axes left edge to the panel label's right edge, in layout
# pixels: far enough left to clear the y tick labels and ylabel, while staying
# inside ALIGNED_GRID_LEFT_PX.
_PANEL_LABEL_LEFT_PX = 36

# Offset from the axes top edge to the panel label's bottom edge, in layout
# pixels: clears the axes title, which sits at title_fontsize + axes_titlepad.
_PANEL_LABEL_TOP_PX = 11

# Uppercase A-Z, then A1, A2, ... A bare chr(65 + i) emits '[', '\', ']' past Z.
_ALPHABET_SIZE = 26


# =============================================================================
# CORE LOGIC
# =============================================================================
def panel_labels(n: int) -> list[str]:
    """Return n unique panel labels: A-Z, then A1, A2, ... for n > 26."""
    return [
        chr(65 + i) if i < _ALPHABET_SIZE else f"A{i - _ALPHABET_SIZE + 1}"
        for i in range(n)
    ]


def grid_panel_size(
    total_width_px: int,
    total_height_px: int,
    n_cols: int,
    n_rows: int,
    decoration_px: int = PANEL_DECORATION_PX,
    min_width: int = 55,
    min_height: int = 55,
) -> tuple[int, int]:
    """Return the (width, height) in pixels for one panel of an n_cols x n_rows grid."""
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")
    if n_rows < 1:
        raise ValueError(f"n_rows must be at least 1, got {n_rows}")

    width = max(min_width, int(total_width_px / n_cols) - decoration_px)
    height = max(min_height, int(total_height_px / n_rows) - decoration_px)

    return width, height


def square_panel_size(
    total_width_px: int,
    n_cols: int,
    decoration_px: int = PANEL_DECORATION_PX,
    min_size: int = 55,
) -> int:
    """Return one edge length for a square panel that fits n_cols across total_width_px.

    cnsplots grows the whole figure's height to fit however many rows a grid
    needs (``multipanel._create_or_update_figure`` sums row heights into the
    final figure height), so budgeting height from a fixed total and n_rows --
    as ``grid_panel_size`` does -- fights that behaviour and yields rectangular
    panels whenever n_cols != n_rows. Deriving one edge from the width budget
    alone keeps panels square regardless of how many rows the grid needs.
    """
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")

    return max(min_size, int(total_width_px / n_cols) - decoration_px)


def aligned_grid_axes(
    n_rows: int,
    n_cols: int,
    *,
    panel_px: int,
    labels: Sequence[str] | None = None,
    left_px: int = ALIGNED_GRID_LEFT_PX,
    right_px: int = ALIGNED_GRID_RIGHT_PX,
    top_px: int = ALIGNED_GRID_TOP_PX,
    bottom_px: int = ALIGNED_GRID_BOTTOM_PX,
) -> list[Axes]:
    """Create an n_rows x n_cols grid of square axes on a fixed pixel pitch, row-major.

    Replaces ``cns.multipanel`` for grids that must align. multipanel positions
    each panel at ``x + margin_left + left_reserve`` where ``left_reserve`` is
    the *measured* rendered width of that panel's own y tick labels, ylabel and
    panel label, so two panels in the same grid column land at different x
    whenever those differ in width (10-15 px in practice for log10 QC grids).
    Every cell here gets the same reserve instead, so column x and row y are
    identical by construction.

    Panels stay square because both edges are ``panel_px``: cnsplots' own
    figure-height growth (``multipanel._create_or_update_figure`` sums row
    heights) is not in play, the figure is sized to the grid up front.
    """
    if n_rows < 1:
        raise ValueError(f"n_rows must be at least 1, got {n_rows}")
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")
    if panel_px < 1:
        raise ValueError(f"panel_px must be at least 1, got {panel_px}")

    cell_w = left_px + panel_px + right_px
    cell_h = top_px + panel_px + bottom_px
    fig_w = n_cols * cell_w
    fig_h = n_rows * cell_h

    cns.figure(width=fig_w, height=fig_h)
    fig = plt.gcf()

    names = list(labels) if labels is not None else panel_labels(n_rows * n_cols)

    axes: list[Axes] = []
    for index in range(n_rows * n_cols):
        row, col = divmod(index, n_cols)
        left = col * cell_w + left_px
        # add_axes takes the bottom edge, while the grid is laid out top-down.
        top = row * cell_h + top_px
        ax = fig.add_axes((
            left / fig_w,
            (fig_h - top - panel_px) / fig_h,
            panel_px / fig_w,
            panel_px / fig_h,
        ))
        cns.setup_ax(ax)
        if index < len(names):
            # plt.gca() is what add_panel_label reads, and add_axes already made
            # ax current; sca keeps that true if a caller reorders these steps.
            plt.sca(ax)
            cns.add_panel_label(
                names[index], pad_left=_PANEL_LABEL_LEFT_PX, pad_top=_PANEL_LABEL_TOP_PX
            )
        axes.append(ax)

    return axes
