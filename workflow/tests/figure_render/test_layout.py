#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for shared figure layout helpers.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pytest

from figure_render._layout import PANEL_DECORATION_PX, grid_panel_size, panel_labels, square_panel_size


# =============================================================================
# TESTS
# =============================================================================
def test_panel_labels_are_letters_below_26() -> None:
    """Assert the first 26 labels are plain uppercase letters."""
    labels = panel_labels(26)

    assert labels[0] == "A"
    assert labels[25] == "Z"
    assert all(label.isalpha() and label.isupper() for label in labels)


def test_panel_labels_do_not_emit_punctuation_past_z() -> None:
    r"""Assert overflow past Z stays alphanumeric, guarding the chr(65+i) bug.

    The old bare ``chr(65 + i)`` produced '[', '\\' and ']' for i >= 26.
    curve_fitting.py renders 32 panels by default, so it was already emitting
    those characters as panel labels.
    """
    labels = panel_labels(32)

    assert len(labels) == 32
    for label in labels:
        assert label.isalnum(), f"Non-alphanumeric panel label: {label!r}"
    assert "[" not in labels
    assert "\\" not in labels
    assert "]" not in labels


def test_panel_labels_are_unique() -> None:
    """Assert no label repeats, so panels stay individually citable."""
    labels = panel_labels(60)

    assert len(set(labels)) == 60


def test_panel_labels_zero_returns_empty() -> None:
    """Assert a zero-panel figure yields no labels rather than raising."""
    assert panel_labels(0) == []


def test_grid_panel_size_subtracts_decoration() -> None:
    """Assert per-panel size is the share of the figure minus decoration."""
    width, height = grid_panel_size(510, 425, n_cols=5, n_rows=3, decoration_px=40)

    assert width == 510 // 5 - 40
    assert height == 425 // 3 - 40


def test_grid_panel_size_enforces_floor() -> None:
    """Assert crowded grids clamp to the minimum rather than going negative."""
    width, height = grid_panel_size(510, 425, n_cols=50, n_rows=50, decoration_px=40)

    assert width == 55
    assert height == 55


def test_grid_panel_size_rejects_zero_cols() -> None:
    """Assert a zero column count is refused instead of dividing by zero."""
    with pytest.raises(ValueError, match="n_cols"):
        grid_panel_size(510, 425, n_cols=0, n_rows=1)


def test_grid_panel_size_reproduces_orientation_values() -> None:
    """Assert the helper matches orientation.py's arithmetic for its real grid.

    orientation.py rendered a 5-timepoint by 3-sample grid at
    JOURNAL_WIDTH_PX=510 / JOURNAL_HEIGHT_PX=425 with PANEL_DECORATION_PX=40.
    Reproducing it exactly is what keeps that figure pixel-identical.
    """
    width, height = grid_panel_size(510, 425, n_cols=5, n_rows=3, decoration_px=PANEL_DECORATION_PX)

    assert width == max(45, int(510 / 5) - 40)
    assert height == max(55, int(425 / 3) - 40)


def test_square_panel_size_derives_from_width_only() -> None:
    """Assert the edge length depends only on n_cols, not on any row count."""
    edge = square_panel_size(510, n_cols=5, decoration_px=40)

    assert edge == 510 // 5 - 40


def test_square_panel_size_enforces_floor() -> None:
    """Assert a crowded row clamps to the minimum rather than going negative."""
    edge = square_panel_size(510, n_cols=50, decoration_px=40)

    assert edge == 55


def test_square_panel_size_rejects_zero_cols() -> None:
    """Assert a zero column count is refused instead of dividing by zero."""
    with pytest.raises(ValueError, match="n_cols"):
        square_panel_size(510, n_cols=0)
