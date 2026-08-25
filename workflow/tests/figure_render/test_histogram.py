#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic histogram grid renderers.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.histogram import (
    render_histogram_grid_figure,
    render_prebinned_histogram_figure,
)


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def metric_frame() -> pd.DataFrame:
    """A frame of arbitrarily named numeric metrics."""
    rng = np.random.default_rng(3)
    return pd.DataFrame({f"metric_{i}": rng.normal(0, 1, 500) for i in range(6)})


@pytest.fixture
def prebinned_frame() -> pd.DataFrame:
    """A pre-binned frame with two groups and five bins each."""
    rows = []
    for sample in ("s1", "s2"):
        for stage in ("early", "late"):
            for edge in range(5):
                rows.append(
                    {
                        "sample": sample,
                        "stage": stage,
                        "bin_left": float(edge),
                        "bin_right": float(edge + 1),
                        "count": 10 * (edge + 1),
                    }
                )
    return pd.DataFrame(rows)


# =============================================================================
# TESTS — raw values mode
# =============================================================================
def test_one_panel_per_value_column(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert the panel count equals the number of value columns requested."""
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "metrics",
        value_columns=["metric_0", "metric_1", "metric_2"], bins=20,
    )

    assert len(plt.gcf().get_axes()) == 3
    assert (tmp_path / "metrics.pdf").exists()


def test_value_columns_are_caller_supplied(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert no metric-name whitelist is applied.

    distribution.py hardcoded a 15-name list ['A','DR','DL','t10',...]; a caller
    naming arbitrary columns must work.
    """
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "arb", value_columns=["metric_5"], bins=10,
    )

    assert len(plt.gcf().get_axes()) == 1


def test_summary_stats_can_be_disabled(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert the per-panel n/mean/std box is optional."""
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "nostats",
        value_columns=["metric_0"], bins=10, show_summary_stats=False,
    )
    # cnsplots' multipanel.panel() always draws the panel corner-letter (e.g.
    # "A") as an ax.text() artist regardless of show_summary_stats, so
    # ax.texts is never empty; assert specifically that the stats box itself
    # is absent instead.
    texts = [t.get_text() for ax in plt.gcf().get_axes() for t in ax.texts]
    assert not any("Mean" in text for text in texts)

    render_histogram_grid_figure(
        metric_frame, tmp_path / "stats",
        value_columns=["metric_0"], bins=10, show_summary_stats=True,
    )
    texts = [t.get_text() for ax in plt.gcf().get_axes() for t in ax.texts]
    assert any("Mean" in text for text in texts)


def test_all_nan_column_renders_placeholder(tmp_path: Path) -> None:
    """Assert a column with no finite values renders a placeholder, not a crash."""
    df = pd.DataFrame({"empty_metric": [np.nan, np.nan], "ok_metric": [1.0, 2.0]})

    render_histogram_grid_figure(
        df, tmp_path / "nan", value_columns=["empty_metric", "ok_metric"], bins=5,
    )

    assert (tmp_path / "nan.pdf").exists()


def test_missing_value_column_raises(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert an absent column is reported by name."""
    with pytest.raises(ValueError, match="absent"):
        render_histogram_grid_figure(
            metric_frame, tmp_path / "bad", value_columns=["absent"], bins=5,
        )


def test_empty_value_columns_does_not_crash(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert requesting no columns returns without raising."""
    render_histogram_grid_figure(metric_frame, tmp_path / "none", value_columns=[], bins=5)


# =============================================================================
# TESTS — pre-binned mode
# =============================================================================
def test_prebinned_bar_heights_equal_stored_counts(
    prebinned_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert stored counts are replayed exactly, with no re-binning.

    This is the core guarantee of the pre-binned mode: the computation layer
    already binned the data, and the renderer must not bin again.
    """
    import matplotlib.pyplot as plt

    render_prebinned_histogram_figure(
        prebinned_frame, tmp_path / "prebinned",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
    )

    first_group = prebinned_frame[
        (prebinned_frame["sample"] == "s1") & (prebinned_frame["stage"] == "early")
    ]
    heights = np.array([p.get_height() for p in plt.gcf().get_axes()[0].patches])
    expected = first_group["count"].to_numpy().astype(float)

    assert np.array_equal(np.sort(heights), np.sort(expected)), (
        f"Bar heights {heights} differ from stored counts {expected}"
    )


def test_marker_only_on_named_column_value(
    prebinned_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the cutoff marker appears only on panels matching marker_on_col_value.

    read_counts.py drew the cutoff line only on the initial-timepoint panel,
    because that is the timepoint the cutoff is applied to.
    """
    import matplotlib.pyplot as plt

    render_prebinned_histogram_figure(
        prebinned_frame, tmp_path / "marker",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
        marker_value=2.5, marker_label="Cutoff", marker_on_col_value="early",
    )

    axes = plt.gcf().get_axes()
    with_marker = sum(1 for ax in axes if ax.get_legend() is not None)

    assert with_marker == 2, f"Expected 2 'early' panels to carry the marker, got {with_marker}"


def test_footer_is_rendered(prebinned_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert figure-level footer lines reach the figure.

    Retention numbers are per sample, not per panel; read_counts.py moved them
    out of a cramped in-axes box that overlapped the histograms.
    """
    import matplotlib.pyplot as plt

    render_prebinned_histogram_figure(
        prebinned_frame, tmp_path / "footer",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
        footer_lines=["s1: 5/10 rows kept"], footer_header="Retention:",
    )

    figure_texts = [t.get_text() for t in plt.gcf().texts]
    assert any("5/10 rows kept" in text for text in figure_texts)


def test_prebinned_all_nan_group_renders_placeholder(tmp_path: Path) -> None:
    """Assert a group whose bin edges are all NaN renders a placeholder panel."""
    df = pd.DataFrame(
        {
            "sample": ["s1", "s1"],
            "stage": ["early", "late"],
            "bin_left": [0.0, np.nan],
            "bin_right": [1.0, np.nan],
            "count": [10, 0],
        }
    )

    render_prebinned_histogram_figure(
        df, tmp_path / "marker_rows",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
    )

    assert (tmp_path / "marker_rows.pdf").exists()


def test_prebinned_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty pre-binned frame returns without raising."""
    empty = pd.DataFrame(columns=["sample", "stage", "bin_left", "bin_right", "count"])

    render_prebinned_histogram_figure(
        empty, tmp_path / "empty",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
    )
