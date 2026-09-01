"""
Shared Rendering Library
========================

House style, panel layout and dual-artifact saving for cnsplots figures.
This is a library module (no main(), no CLI) per python-script-conventions.

Every uniform grid in the project goes through ``grid_axes``, which is the
cnsplots-documented matplotlib grid pattern: ``cns.figure()`` sizes the page,
``fig.subplots()`` allots equal GridSpec cells, ``set_box_aspect`` keeps them
square and ``tight_layout`` fits the decorations. Heterogeneous figures use
``cns.multipanel`` directly.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


# =============================================================================
# CONSTANTS
# =============================================================================
# Page width in cnsplots layout pixels (72 per inch). Taken from cnsplots' own
# multipanel default so figures built either way share one page width.
JOURNAL_WIDTH_PX = int(cns.settings.multipanel_max_width)
JOURNAL_HEIGHT_PX = 425  # ~150 mm, 4:3 aspect

# Panel-label offsets from the axes corner, in layout pixels. cnsplots defaults
# both to 0, which puts the letter on top of the y tick labels; measured over a
# 2x6 log grid, this is the smallest pair with no overlap against either the
# tick labels or the axes title.
PANEL_LABEL_PAD_LEFT_PX = 18
PANEL_LABEL_PAD_TOP_PX = 6

# Uppercase A-Z, then A1, A2, ... A bare chr(65 + i) emits '[', '\\', ']' past Z.
_ALPHABET_SIZE = 26

# One qualitative palette for every figure in the project. A reader assumes
# consistent colour means the same thing, so the palette name lives here alone
# and no module may name a colour of its own.
HOUSE_PALETTE = "Cell"

# Semantic colour roles, resolved from HOUSE_PALETTE on demand rather than
# written as literals, so re-pointing HOUSE_PALETTE re-colours every figure.
# Indices are stable positions in the palette, not arbitrary picks:
#   0, 1  the first two entries, for a two-series contrast
#   2,3,4 the greyscale-safe trio (min Rec. 601 luma gap 0.145, against 0.026
#         for 0/1/2) — figures get printed and photocopied
_OBSERVED_INDEX = 0
_FITTED_INDEX = 1
_SERIES_INDICES = (2, 3, 4)

# Reference lines, diagonals, cutoff markers, bar edges and annotations are
# furniture, not data. Keeping them off-palette leaves the palette to mean data.
FURNITURE_COLOR = cns.GRAY

# Density is a magnitude, so it needs a sequential map. Installed as
# cns.settings.palette_seq by apply_house_style(), which lands it in
# rcParams['image.cmap'], so unqualified scatter/imshow calls inherit it.
# Overrides cnsplots' 'gnuplot', which is non-monotonic in perceptual lightness:
# equal density steps must read as equal colour steps.
DENSITY_CMAP = "viridis"


# =============================================================================
# CORE LOGIC
# =============================================================================
def apply_house_style() -> None:
    """Configure cnsplots with project house style: Cell palette, viridis density, Arial font."""
    # Remove Helvetica to prevent findfont warnings; Arial is installed
    cns.settings.font_sans_serif = ("Arial", "DejaVu Sans")
    cns.settings.panel_label_fontname = "Arial"

    # Fonttype 42 embeds TrueType, editable in Illustrator
    cns.settings.pdf_fonttype = 42

    cns.settings.panel_pad_left = PANEL_LABEL_PAD_LEFT_PX
    cns.settings.panel_pad_top = PANEL_LABEL_PAD_TOP_PX

    # Both are needed: setup_matplotlib sets the rcParams cycle for bare
    # matplotlib calls, palette_qual is what cns.* plot functions and
    # multipanel.panel() read when no explicit color_cycle is passed.
    cns.settings.palette_qual = HOUSE_PALETTE
    cns.settings.palette_seq = DENSITY_CMAP
    cns.setup_matplotlib(color_cycle=HOUSE_PALETTE)


def panel_labels(n: int) -> list[str]:
    """Return n unique panel labels: A-Z, then A1, A2, ... for n > 26."""
    return [
        chr(65 + i) if i < _ALPHABET_SIZE else f"A{i - _ALPHABET_SIZE + 1}"
        for i in range(n)
    ]


def grid_axes(
    n_rows: int,
    n_cols: int,
    *,
    labels: Sequence[str] | None = None,
    width_px: int = JOURNAL_WIDTH_PX,
    square: bool = True,
) -> list[Axes]:
    """Create an n_rows x n_cols grid of aligned axes, row-major, and label each one.

    GridSpec allots equal cells regardless of how wide a panel's tick labels
    render, so column x and row y agree by construction -- which cns.multipanel
    cannot do, since it offsets every panel by its own measured decoration
    width. ``square`` additionally pins each axes box to a 1:1 aspect, so a grid
    whose row and column counts differ still holds square panels.

    Callers hide the cells they do not fill; ``tight_layout`` then reclaims the
    space. Run it after all panels are drawn, since it measures rendered text.
    """
    if n_rows < 1:
        raise ValueError(f"n_rows must be at least 1, got {n_rows}")
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")

    cell_px = width_px / n_cols
    cns.figure(width=width_px, height=int(cell_px * n_rows))
    fig = plt.gcf()
    grid = fig.subplots(n_rows, n_cols, squeeze=False)

    names = list(labels) if labels is not None else panel_labels(n_rows * n_cols)

    axes: list[Axes] = []
    for index in range(n_rows * n_cols):
        ax = grid[index // n_cols][index % n_cols]
        if square:
            ax.set_box_aspect(1)
        if index < len(names):
            plt.sca(ax)
            cns.add_panel_label(names[index])
        axes.append(ax)

    return axes


def house_colors(indices: Sequence[int]) -> list[str]:
    """Return the house-palette hex colors at the given positions."""
    return cns.get_hexcolors_from_apalette(list(indices), HOUSE_PALETTE)


def observed_fitted_colors() -> tuple[str, str]:
    """Return the (observed, fitted) pair for measurement-versus-model figures."""
    observed, fitted = house_colors((_OBSERVED_INDEX, _FITTED_INDEX))
    return observed, fitted


def series_colors(n: int) -> list[str]:
    """Return n greyscale-separated house colors for overlaid series on one axes."""
    if n > len(_SERIES_INDICES):
        # Beyond the curated trio, luma separation is no longer guaranteed, so
        # fall back to palette order and let the caller carry the distinction on
        # a second channel (marker, dash) if print legibility matters.
        return house_colors(range(n))
    return house_colors(_SERIES_INDICES[:n])


def save_dual(stem: Path | str) -> None:
    """Save the current figure as a journal PDF plus a review PNG, given a stem with no extension."""
    # Append rather than replace suffixes: Path.with_suffix()/.stem would truncate
    # at the first dot, so a stem like "HD1328-4.YES0_corr" would silently collapse
    # to "HD1328-4" and two samples would overwrite each other's figures.
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    cns.savefig(stem.parent / f"{stem.name}.pdf")
    cns.savefig(stem.parent / f"{stem.name}.review.png")
