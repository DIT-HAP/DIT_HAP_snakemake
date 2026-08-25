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
        total_column="total",
        part_label="Present",
        whole_label="Absent",
        xlabel="Bucket",
        ylabel="Present (%)",
        title="Presence by bucket",
    )


# =============================================================================
# TESTS
# =============================================================================
def test_panel_count_is_one_plus_categories(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert one bar overview panel plus one donut per category."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "comp", **render_kwargs)

    assert len(plt.gcf().get_axes()) == 1 + len(composition_frame)
    assert (tmp_path / "comp.pdf").exists()


def test_bar_panel_is_percentage_scaled(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert the overview panel's y-axis is pinned to 0-100."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "bar", **render_kwargs)

    assert plt.gcf().get_axes()[0].get_ylim() == (0.0, 100.0)


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
    """Assert donut panels follow the frame's row order, not alphabetical order."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "order", **render_kwargs)

    donut_titles = [ax.get_title() for ax in plt.gcf().get_axes()[1:]]
    assert donut_titles[0].startswith("low")
    assert donut_titles[1].startswith("mid")
    assert donut_titles[2].startswith("high")


def test_missing_column_raises(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert an absent column is reported by name."""
    kwargs = render_kwargs | {"percentage_column": "absent_pct"}

    with pytest.raises(ValueError, match="absent_pct"):
        render_composition_figure(composition_frame, tmp_path / "bad", **kwargs)
