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
    render_grouped_histogram_figure,
    render_histogram_grid_figure,
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
def grouped_frame() -> pd.DataFrame:
    """A long-format frame with two row groups, two column groups, and separable values."""
    rows = []
    for sample in ("s1", "s2"):
        for stage in ("early", "late"):
            base = ({"s1": 1.0, "s2": 100.0}[sample]) * ({"early": 1.0, "late": 10.0}[stage])
            values = base * np.geomspace(0.5, 2.0, 50)
            rows.append(pd.DataFrame({"sample": sample, "stage": stage, "value": values}))
    return pd.concat(rows, ignore_index=True)


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
# TESTS — grouped mode
# =============================================================================
def test_grouped_one_panel_per_row_col_pair(
    grouped_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert panels are allocated per row_key/col_key pair."""
    import matplotlib.pyplot as plt

    render_grouped_histogram_figure(
        grouped_frame, tmp_path / "grouped",
        value_column="value", row_key="sample", col_key="stage",
        bins=20,
    )

    assert len(plt.gcf().get_axes()) == 4, "Two samples x two stages should yield four panels"
    assert (tmp_path / "grouped.pdf").exists()


def test_grouped_marker_only_on_named_col_value(
    grouped_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the cutoff marker appears only on panels matching marker_on_col_value.

    read_counts.py drew the cutoff line only on the initial-timepoint panel,
    because that is the timepoint the cutoff is applied to.
    """
    import matplotlib.pyplot as plt

    render_grouped_histogram_figure(
        grouped_frame, tmp_path / "marker",
        value_column="value", row_key="sample", col_key="stage",
        marker_value=5.0, marker_label="Cutoff", marker_on_col_value="early",
    )

    axes = plt.gcf().get_axes()
    with_marker = sum(1 for ax in axes if ax.get_legend() is not None)

    assert with_marker == 2, f"Expected 2 'early' panels to carry the marker, got {with_marker}"


def test_grouped_log_scale_sets_log_axis(
    grouped_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert log_scale mode lays a true log axis with logspace bin edges.

    Bin edges must be exactly np.logspace over each group's own value range,
    which on a log axis renders bars geometrically identical to equal-width
    log10 bins — the property that lets the single-stage figure replace the
    pre-binned one without changing bar shapes.
    """
    import matplotlib.pyplot as plt

    render_grouped_histogram_figure(
        grouped_frame, tmp_path / "logscale",
        value_column="value", row_key="sample", col_key="stage",
        bins=10, log_scale=True,
    )

    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log"

    group = grouped_frame[grouped_frame["sample"] == "s1"]
    expected_edges = np.logspace(np.log10(group["value"].min()), np.log10(group["value"].max()), 11)
    edges = ax.patches[0].get_x()
    assert np.isclose(edges, expected_edges[0]), (
        f"First bar edge {edges} does not match logspace start {expected_edges[0]}"
    )


def test_grouped_footer_is_rendered(grouped_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert figure-level footer lines reach the figure.

    Retention numbers are per sample, not per panel; read_counts.py moved them
    out of a cramped in-axes box that overlapped the histograms.
    """
    import matplotlib.pyplot as plt

    render_grouped_histogram_figure(
        grouped_frame, tmp_path / "footer",
        value_column="value", row_key="sample", col_key="stage",
        footer_lines=["s1: 5/10 rows kept"], footer_header="Retention:",
    )

    figure_texts = [t.get_text() for t in plt.gcf().texts]
    assert any("5/10 rows kept" in text for text in figure_texts)


def test_grouped_all_nonpositive_group_renders_placeholder(tmp_path: Path) -> None:
    """Assert a group whose values are all non-positive (log mode) renders a placeholder panel."""
    df = pd.DataFrame(
        {
            "sample": ["s1", "s1"],
            "stage": ["early", "late"],
            "value": [1.0, -3.0],
        }
    )
    df.loc[df.index[-1], "value"] = -7.0

    render_grouped_histogram_figure(
        df, tmp_path / "nonpositive",
        value_column="value", row_key="sample", col_key="stage",
        log_scale=True,
    )

    assert (tmp_path / "nonpositive.pdf").exists()


def test_grouped_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty frame returns without raising."""
    empty = pd.DataFrame(columns=["sample", "stage", "value"])

    render_grouped_histogram_figure(
        empty, tmp_path / "empty",
        value_column="value", row_key="sample", col_key="stage",
    )
