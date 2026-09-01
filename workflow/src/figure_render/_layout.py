"""Shared panel layout helpers for figure rendering.

Panel-letter generation and grid sizing were duplicated across five figure
modules in two incompatible variants; this module is the single implementation.

Two layout strategies live here. ``square_panel_size`` / ``grid_panel_size``
feed ``cns.multipanel``, which flows panels left-to-right and derives each
panel's origin from the *rendered* width of that panel's own y-axis
decorations, so panels in the same grid column start at different x whenever
their tick labels differ in width. ``aligned_grid_axes`` instead lays cells out
on a ``GridSpec``, which allots equal cells regardless of decoration width --
that is what a QC grid comparing the same quantity across samples needs.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence

import cnsplots as cns
import matplotlib.pyplot as plt
from loguru import logger
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
# Undersizing does not shift any axes (GridSpec cells are equal by
# construction); it only lets text spill past the figure edge, which savefig's
# tight bbox then absorbs.
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
    panel_px: int | None = None,
    page_width_px: int | None = None,
    labels: Sequence[str] | None = None,
    cells: Sequence[tuple[int, int]] | None = None,
    left_px: int = ALIGNED_GRID_LEFT_PX,
    right_px: int = ALIGNED_GRID_RIGHT_PX,
    top_px: int = ALIGNED_GRID_TOP_PX,
    bottom_px: int = ALIGNED_GRID_BOTTOM_PX,
    min_panel_px: int = 55,
) -> list[Axes]:
    """Create an n_rows x n_cols grid of square axes in a GridSpec, row-major.

    Replaces ``cns.multipanel`` for grids that must align. multipanel positions
    each panel at ``x + margin_left + left_reserve`` where ``left_reserve`` is
    the *measured* rendered width of that panel's own y tick labels, ylabel and
    panel label, so two panels in the same grid column land at different x
    whenever those differ in width (10-15 px in practice for log10 QC grids).
    ``GridSpec`` allots equal cells regardless of decoration width, so column x
    and row y are identical by construction, and nothing recomputes the layout
    on draw the way multipanel's draw handler does.

    Panels stay square because the figure is sized to the grid up front and the
    reserves are subtracted uniformly: each cell's axes box is panel_px on both
    edges.

    Pass either ``panel_px`` (explicit size) or ``page_width_px`` (derive size
    from page width and n_cols). When ``page_width_px`` is given, panel size is
    ``max(min_panel_px, (page_width_px // n_cols) - (left_px + right_px))``,
    and the function logs a warning if the resulting figure width exceeds
    ``page_width_px`` due to the min_panel_px floor.

    ``cells`` limits which grid positions get axes. Pass a sequence of (row, col)
    tuples (0-indexed); omitted cells are skipped. When None, all n_rows * n_cols
    cells are filled.
    """
    if n_rows < 1:
        raise ValueError(f"n_rows must be at least 1, got {n_rows}")
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")
    if panel_px is None and page_width_px is None:
        raise ValueError("Must provide either panel_px or page_width_px")
    if panel_px is not None and page_width_px is not None:
        raise ValueError("Cannot provide both panel_px and page_width_px")

    if page_width_px is not None:
        h_reserve = left_px + right_px
        panel_px = max(min_panel_px, (page_width_px // n_cols) - h_reserve)
        fig_width = n_cols * (panel_px + h_reserve)
        if fig_width > page_width_px:
            logger.warning(
                f"Grid figure width {fig_width} px exceeds page width {page_width_px} px "
                f"(n_cols={n_cols}, panel={panel_px}, min_panel_px={min_panel_px})"
            )

    assert panel_px is not None  # one of the two size arguments is always set above
    if panel_px < 1:
        raise ValueError(f"panel_px must be at least 1, got {panel_px}")

    cell_w = left_px + panel_px + right_px
    cell_h = top_px + panel_px + bottom_px
    fig_w = n_cols * cell_w
    fig_h = n_rows * cell_h

    cns.figure(width=fig_w, height=fig_h)
    fig = plt.gcf()

    # GridSpec spacing is expressed relative to the *axes* size, not the cell, so
    # the inter-cell gap is the reserve pair over panel_px. Outer bounds keep one
    # half-reserve outside the first/last axes, which makes every cell's total
    # footprint exactly cell_w x cell_h.
    gridspec = fig.add_gridspec(
        n_rows,
        n_cols,
        left=left_px / fig_w,
        right=1 - right_px / fig_w,
        top=1 - top_px / fig_h,
        bottom=bottom_px / fig_h,
        wspace=(left_px + right_px) / panel_px,
        hspace=(top_px + bottom_px) / panel_px,
    )

    names = list(labels) if labels is not None else panel_labels(n_rows * n_cols)

    if cells is None:
        cells = [(r, c) for r in range(n_rows) for c in range(n_cols)]

    cell_to_label = {(r, c): names[i] for i, (r, c) in enumerate(cells) if i < len(names)}

    axes: list[Axes] = []
    for row, col in cells:
        if row >= n_rows or col >= n_cols:
            raise ValueError(f"Cell ({row}, {col}) outside grid bounds ({n_rows}, {n_cols})")

        ax = fig.add_subplot(gridspec[row, col])

        label = cell_to_label.get((row, col))
        if label:
            plt.sca(ax)
            cns.add_panel_label(
                label, pad_left=_PANEL_LABEL_LEFT_PX, pad_top=_PANEL_LABEL_TOP_PX
            )
        axes.append(ax)

    return axes
