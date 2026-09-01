#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the shared figure layout helpers in figures.py.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-01
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from figures import JOURNAL_WIDTH_PX, apply_house_style, grid_axes, panel_labels


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture(autouse=True)
def _house_style():
    """Apply the house style and close figures, so each test starts clean."""
    apply_house_style()
    yield
    plt.close("all")


def _overlap(first, second) -> float:
    """Return the smaller overlap dimension of two bboxes, 0 when they are disjoint."""
    dx = min(first.x1, second.x1) - max(first.x0, second.x0)
    dy = min(first.y1, second.y1) - max(first.y0, second.y0)
    return min(dx, dy) if dx > 0 and dy > 0 else 0.0


# =============================================================================
# TESTS — panel labels
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
    curves.py renders 32 panels by default, so it was already emitting those
    characters as panel labels.
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


# =============================================================================
# TESTS — grid_axes
# =============================================================================
def test_grid_axes_returns_row_major_cells() -> None:
    """Assert the returned list is row-major: the first n_cols entries share a row."""
    axes = grid_axes(2, 3)

    assert len(axes) == 6
    top_row = {round(ax.get_position().y0, 6) for ax in axes[:3]}
    bottom_row = {round(ax.get_position().y0, 6) for ax in axes[3:]}
    assert len(top_row) == 1
    assert len(bottom_row) == 1
    assert top_row != bottom_row


@pytest.mark.parametrize(("n_rows", "n_cols"), [(1, 3), (2, 6), (3, 4), (3, 6), (4, 4)])
def test_grid_axes_aligns_columns_and_rows(n_rows: int, n_cols: int) -> None:
    """Assert every panel in a column shares one x0, and every panel in a row one y0.

    A QC grid compares the same quantity across samples, so the axes boxes have
    to line up for the eye to compare them. cns.multipanel cannot do this: it
    offsets each panel by the rendered width of that panel's own y tick labels
    and ylabel. GridSpec allots equal cells, so alignment holds by construction.
    """
    axes = grid_axes(n_rows, n_cols)

    # Deliberately uneven y magnitudes: this is what drifts under multipanel.
    for index, ax in enumerate(axes):
        ax.plot([0, 1], [0, 10 ** ((index % n_cols) * 2)])

    fig = plt.gcf()
    fig.tight_layout()
    fig.canvas.draw()
    boxes = [ax.get_window_extent() for ax in axes]

    for col in range(n_cols):
        x0s = [boxes[row * n_cols + col].x0 for row in range(n_rows)]
        assert max(x0s) - min(x0s) == pytest.approx(0, abs=0.01), (
            f"Column {col} panels start at different x: {x0s}"
        )
    for row in range(n_rows):
        y0s = [boxes[row * n_cols + col].y0 for col in range(n_cols)]
        assert max(y0s) - min(y0s) == pytest.approx(0, abs=0.01), (
            f"Row {row} panels start at different y: {y0s}"
        )


@pytest.mark.parametrize(("n_rows", "n_cols"), [(2, 6), (3, 4), (1, 5)])
def test_grid_axes_panels_are_square(n_rows: int, n_cols: int) -> None:
    """Assert axes stay square even when the row and column counts differ.

    set_box_aspect(1) constrains the axes box itself, so tight_layout resizing
    the cell cannot turn a panel into a rectangle.
    """
    axes = grid_axes(n_rows, n_cols)
    for ax in axes:
        ax.plot([0, 1], [0, 1])

    fig = plt.gcf()
    fig.tight_layout()
    fig.canvas.draw()

    for ax in axes:
        box = ax.get_window_extent()
        assert box.width == pytest.approx(box.height, abs=1.0), (
            f"Panel is not square: width={box.width:.1f} height={box.height:.1f}"
        )


def test_grid_axes_square_false_lets_the_cell_set_the_shape() -> None:
    """Assert square=False leaves the axes box free to be wider than it is tall."""
    axes = grid_axes(1, 2, square=False)
    for ax in axes:
        ax.plot([0, 1], [0, 1])

    fig = plt.gcf()
    fig.tight_layout()
    fig.canvas.draw()

    box = axes[0].get_window_extent()
    assert box.width != pytest.approx(box.height, abs=1.0)


def test_grid_axes_figure_width_is_the_page_width() -> None:
    """Assert the figure is sized to the page, so a row fills the journal column."""
    grid_axes(2, 3)

    width_px = plt.gcf().get_size_inches()[0] * 72
    assert width_px == pytest.approx(JOURNAL_WIDTH_PX, abs=1.0)


def test_grid_axes_labels_every_panel_by_default() -> None:
    """Assert each cell carries a panel letter when no labels are passed."""
    axes = grid_axes(2, 2)

    assert [text.get_text() for ax in axes for text in ax.texts] == ["A", "B", "C", "D"]


def test_grid_axes_labels_only_as_many_panels_as_given() -> None:
    """Assert a short label list leaves the trailing cells unlabelled, for hiding."""
    axes = grid_axes(2, 2, labels=["A", "B", "C"])

    assert [text.get_text() for ax in axes for text in ax.texts] == ["A", "B", "C"]
    assert len(axes[3].texts) == 0


def test_grid_axes_panel_label_clears_the_y_tick_labels() -> None:
    """Assert the panel letter overlaps neither the y tick labels nor the title.

    cnsplots defaults panel_pad_left/top to 0, which puts the letter on top of
    the tick labels; apply_house_style() installs the measured offsets instead.
    """
    axes = grid_axes(2, 3)
    for index, ax in enumerate(axes):
        ax.plot([0, 1], [0, 10 ** ((index % 3) * 2)])
        ax.set_title("YES0")
        ax.set_ylabel("PBR")

    fig = plt.gcf()
    fig.tight_layout()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for ax in axes:
        label_box = ax.texts[0].get_window_extent(renderer)
        assert _overlap(label_box, ax.title.get_window_extent(renderer)) == 0.0
        for tick in ax.get_yticklabels():
            if tick.get_text():
                assert _overlap(label_box, tick.get_window_extent(renderer)) == 0.0


@pytest.mark.parametrize(("n_rows", "n_cols"), [(0, 1), (1, 0)])
def test_grid_axes_rejects_empty_grids(n_rows: int, n_cols: int) -> None:
    """Assert a zero row or column count is refused instead of dividing by zero."""
    with pytest.raises(ValueError, match="at least 1"):
        grid_axes(n_rows, n_cols)
