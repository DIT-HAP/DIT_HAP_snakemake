"""Shared panel layout helpers for figure rendering.

Panel-letter generation and grid sizing were duplicated across five figure
modules in two incompatible variants; this module is the single implementation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================

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
