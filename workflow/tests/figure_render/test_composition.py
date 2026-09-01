#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic composition (bar plus donuts) renderer.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pandas as pd
import pytest

from figure_render.composition import render_composition_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def composition_frame() -> pd.DataFrame:
    """Three categories with complementary part/whole counts."""
    return pd.DataFrame(
        {
            "bucket": ["low", "mid", "high"],
            "present": [80, 50, 20],
            "absent": [20, 50, 80],
            "total": [100, 100, 100],
            "present_pct": [80.0, 50.0, 20.0],
        }
    )


@pytest.fixture
def render_kwargs() -> dict[str, str]:
    """The column mapping shared by most tests."""
    return dict(
        category_column="bucket",
        percentage_column="present_pct",
        part_column="present",
        whole_column="absent",
        part_label="Present",
        whole_label="Absent",
        xlabel="Present (%)",
        ylabel="Bucket",
        title="Presence by bucket",
    )


# =============================================================================
# TESTS
# =============================================================================
def _content_axes(fig) -> list:
    """Return the figure's data-bearing axes, dropping the donut grid's host rectangle.

    The donut block is a GridSpec inside an empty host panel, which is how its
    cells come out aligned. The host holds no data and cannot be removed —
    multipanel's draw handler still reaches for it — so it is skipped by gid.
    """
    from figure_render.composition import DONUT_GRID_HOST_GID

    return [ax for ax in fig.get_axes() if ax.get_gid() != DONUT_GRID_HOST_GID]


def test_panel_count_is_one_plus_categories(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert one bar overview panel plus one donut per category."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "comp", **render_kwargs)

    assert len(_content_axes(plt.gcf())) == 1 + len(composition_frame)
    assert (tmp_path / "comp.pdf").exists()


def test_bar_panel_is_percentage_scaled(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert the overview panel's x-axis (horizontal bars) is pinned to 0-100."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "bar", **render_kwargs)

    assert plt.gcf().get_axes()[0].get_xlim() == (0.0, 100.0)


def test_percentage_labels_are_annotated(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert each bar carries its percentage as text."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "annot", **render_kwargs)

    texts = [t.get_text() for t in plt.gcf().get_axes()[0].texts]
    assert "80.0%" in texts


def test_zero_total_category_renders_placeholder(
    render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert a category with no items renders a placeholder instead of an empty donut."""
    df = pd.DataFrame(
        {
            "bucket": ["empty"],
            "present": [0],
            "absent": [0],
            "total": [0],
            "present_pct": [0.0],
        }
    )

    render_composition_figure(df, tmp_path / "zero", **render_kwargs)

    assert (tmp_path / "zero.pdf").exists()


def test_empty_frame_renders_placeholder(
    render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert an empty frame still writes artifacts.

    coverage.py rendered a single placeholder panel and saved, so the Snakemake
    output always exists.
    """
    empty = pd.DataFrame(columns=["bucket", "present", "absent", "total", "present_pct"])

    render_composition_figure(empty, tmp_path / "empty", **render_kwargs)

    assert (tmp_path / "empty.pdf").exists(), "Snakemake requires the output to exist"


def test_category_order_is_preserved(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert donut panels follow the frame's row order, not alphabetical order.

    Categories live in the donut hole now that the title band is gone; the
    per-panel texts also hold the panel letter, so match by category name.
    """
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "order", **render_kwargs)

    donut_texts = [
        ", ".join(t.get_text() for t in ax.texts) for ax in _content_axes(plt.gcf())[1:]
    ]
    assert "low" in donut_texts[0]
    assert "mid" in donut_texts[1]
    assert "high" in donut_texts[2]


def test_missing_column_raises(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert an absent column is reported by name."""
    kwargs = render_kwargs | {"percentage_column": "absent_pct"}

    with pytest.raises(ValueError, match="absent_pct"):
        render_composition_figure(composition_frame, tmp_path / "bad", **kwargs)


@pytest.fixture
def four_category_frame() -> pd.DataFrame:
    """Four categories, so the donuts wrap into a 2x2 grid with a row below."""
    return pd.DataFrame(
        {
            "bucket": ["low", "mid", "high", "top"],
            "present": [80, 50, 20, 10],
            "absent": [20, 50, 80, 90],
            "total": [100, 100, 100, 100],
            "present_pct": [80.0, 50.0, 20.0, 10.0],
        }
    )


def test_rows_are_separated_by_the_configured_gap(
    four_category_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert consecutive rows are separated by at least ROW_GAP_PX of clear space.

    multipanel reserves layout space for measured *top* decorations only, so
    without an explicit margin_bottom the bar's xlabel rendered on top of the
    next row's panel labels and each donut's legend rendered inside the row
    below it. Asserting a positive gap rather than mere non-overlap keeps rows
    reading as separate blocks.
    """
    import matplotlib.pyplot as plt

    from figure_render.composition import ROW_GAP_PX

    render_composition_figure(four_category_frame, tmp_path / "spill", **render_kwargs)

    fig = plt.gcf()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes = _content_axes(fig)
    tights = [ax.get_tightbbox(renderer) for ax in axes]
    # get_tightbbox is in device px; multipanel sizes in 72-dpi layout px.
    scale = fig.dpi / 72

    rows: dict[float, list[int]] = {}
    for index, ax in enumerate(axes):
        rows.setdefault(round(ax.get_window_extent().y0, 1), []).append(index)
    ordered = [rows[y] for y in sorted(rows, reverse=True)]
    assert len(ordered) == 3, f"Expected bar row plus two donut rows, got {len(ordered)}"

    for upper, lower in zip(ordered, ordered[1:], strict=False):
        gap = (min(tights[i].y0 for i in upper) - max(tights[i].y1 for i in lower)) / scale
        assert gap >= ROW_GAP_PX - 1, (
            f"Rows are {gap:.1f} layout px apart, below the configured {ROW_GAP_PX} px gap"
        )

    # The last row's legend must also stay on the canvas, not hang off the bottom.
    assert min(tights[i].y0 for i in ordered[-1]) >= 0, "Bottom row decorations spill off the figure"


def test_donut_legend_sits_close_under_its_ring(
    four_category_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert each legend sits LEGEND_RING_GAP_PX below its own ring, not the axes box.

    cns.donutplot anchors a bottom legend relative to the axes box, but the ring
    stops at 90% of it, so the default offset left ~20 layout px of dead space
    and the legend read as belonging to the panel below.
    """
    import matplotlib.pyplot as plt

    from figure_render.donut import LEGEND_RING_GAP_PX

    render_composition_figure(four_category_frame, tmp_path / "legend_gap", **render_kwargs)

    fig = plt.gcf()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    scale = fig.dpi / 72

    donuts = [ax for ax in fig.get_axes() if ax.get_legend() is not None]
    assert len(donuts) == len(four_category_frame), "Every donut panel must carry a legend"

    for ax in donuts:
        ring_bottom = min(patch.get_window_extent(renderer).y0 for patch in ax.patches)
        legend_top = ax.get_legend().get_window_extent(renderer).y1
        gap = (ring_bottom - legend_top) / scale
        assert gap == pytest.approx(LEGEND_RING_GAP_PX, abs=1.0), (
            f"Legend sits {gap:.1f} layout px under the ring, expected {LEGEND_RING_GAP_PX}"
        )
