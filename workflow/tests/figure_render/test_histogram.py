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

    visible = [ax for ax in plt.gcf().get_axes() if ax.get_visible()]
    assert len(visible) == 3
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

    visible = [ax for ax in plt.gcf().get_axes() if ax.get_visible()]
    assert len(visible) == 1


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

    visible = [ax for ax in plt.gcf().get_axes() if ax.get_visible()]
    assert len(visible) == 4, "Two samples x two stages should yield four panels"
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


def test_grouped_row_label_in_first_column_ylabel(grouped_frame, tmp_path):
    """Assert the sample name rides the first column's ylabel and titles carry only the col value."""
    import matplotlib.pyplot as plt

    render_grouped_histogram_figure(
        grouped_frame, tmp_path / "rowlabel",
        value_column="value", row_key="sample", col_key="stage",
    )

    axes = plt.gcf().get_axes()
    first_col_axes = axes[0::2]
    second_col_axes = axes[1::2]

    for ax in first_col_axes:
        assert ax.get_ylabel().startswith(("s1\n", "s2\n")), f"First-column ylabel should embed the row label: {ax.get_ylabel()!r}"
        assert "Frequency" in ax.get_ylabel()
    for ax in second_col_axes:
        assert "\n" not in ax.get_ylabel(), f"Non-first columns keep a plain ylabel: {ax.get_ylabel()!r}"
    for ax, expected in zip(axes, ["early", "late", "early", "late"], strict=True):
        assert ax.get_title() == expected, f"Title should be only the col value, got {ax.get_title()!r}"


def test_grouped_upper_quantile_drops_outliers(tmp_path):
    """Assert values above each group's quantile are excluded from binning.

    A lone extreme value must not stretch its panel's bars or the shared x
    range; with quantile=0.9 on 10 ascending values, the max is dropped.
    """
    import matplotlib.pyplot as plt

    df = pd.DataFrame({
        "sample": ["s1"] * 10,
        "stage": ["T0"] * 10,
        "value": np.geomspace(1.0, 1_000_000.0, 10),  # last value = extreme outlier
    })

    render_grouped_histogram_figure(
        df, tmp_path / "clipped",
        value_column="value", row_key="sample", col_key="stage",
        bins=5, log_scale=True, upper_quantile=0.9,
    )

    ax = plt.gcf().get_axes()[0]
    # Rightmost bar edge must stop near p90 (≈129155), not at the outlier 1e6.
    right_edge = max(p.get_x() + p.get_width() for p in ax.patches if p.get_height() > 0)
    assert right_edge < 400_000, f"Outlier should be dropped, but bars extend to {right_edge:.3g}"
    assert (tmp_path / "clipped.pdf").exists()


def test_grouped_shared_x_range_aligns_panels(tmp_path):
    """Assert shared log edges span the union of groups, identical across panels."""
    import matplotlib.pyplot as plt

    df = pd.DataFrame({
        "sample": ["s1"] * 4 + ["s2"] * 4,
        "stage": ["T0"] * 8,
        "value": [1.0, 2.0, 4.0, 8.0] + [16.0, 32.0, 64.0, 128.0],
    })

    render_grouped_histogram_figure(
        df, tmp_path / "sharedx",
        value_column="value", row_key="sample", col_key="stage",
        bins=7, log_scale=True, share_x_range=True,
    )

    axes = plt.gcf().get_axes()
    lefts = {round(ax.get_xlim()[0], 9) for ax in axes}
    rights = {round(ax.get_xlim()[1], 9) for ax in axes}
    assert len(lefts) == 1 and len(rights) == 1, f"x limits should match across panels: {lefts}, {rights}"
    # Bar EDGES (not the padded axis limits) must span exactly the global range.
    edges = sorted({p.get_x() for ax in axes for p in ax.patches} | {p.get_x() + p.get_width() for ax in axes for p in ax.patches})
    assert np.isclose(edges[0], df["value"].min()), f"First bar edge should be the global min: {edges[0]}"
    assert np.isclose(edges[-1], df["value"].max()), f"Last bar edge should be the global max: {edges[-1]}"


def test_grouped_unshared_x_ranges_differ(tmp_path):
    """Assert share_x_range=False preserves per-group autoscaling."""
    import matplotlib.pyplot as plt

    df = pd.DataFrame({
        "sample": ["s1"] * 4 + ["s2"] * 4,
        "stage": ["T0"] * 8,
        "value": [1.0, 2.0, 4.0, 8.0] + [16.0, 32.0, 64.0, 128.0],
    })

    render_grouped_histogram_figure(
        df, tmp_path / "freex",
        value_column="value", row_key="sample", col_key="stage",
        bins=3, log_scale=True, share_x_range=False,
    )

    xlims = [(ax.get_xlim()) for ax in plt.gcf().get_axes()]
    assert not np.isclose(xlims[0][1], xlims[1][1]), "Unshared ranges should remain per-panel"


def test_grid_shared_y_applies_uniform_ylim(metric_frame, tmp_path):
    """Assert grid mode's shared y makes every panel's ylim top equal."""
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "sharedy",
        value_columns=["metric_0", "metric_1"], bins=20, share_y_range=True,
    )

    visible = [ax for ax in plt.gcf().get_axes() if ax.get_visible()]
    tops = [ax.get_ylim()[1] for ax in visible]
    assert len(set(np.round(tops, 6))) == 1, f"Shared-y panels should share one y top: {tops}"


def test_grid_y_is_per_panel_by_default(tmp_path):
    """Assert grid mode autoscales each metric's y independently.

    Curve-fitting metrics span unrelated magnitudes (R² in [0,1] vs AIC in the
    thousands), so one shared y top flattens most panels to nothing.
    """
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(11)
    dense = rng.normal(0, 1, 500)
    sparse = np.full(500, np.nan)
    sparse[:50] = rng.normal(0, 1, 50)
    df = pd.DataFrame({"narrow": dense, "wide": sparse})

    render_histogram_grid_figure(
        df, tmp_path / "freey", value_columns=["narrow", "wide"], bins=20,
    )

    visible = [ax for ax in plt.gcf().get_axes() if ax.get_visible()]
    tops = [ax.get_ylim()[1] for ax in visible]
    assert len(set(np.round(tops, 6))) == 2, f"Panels should autoscale independently: {tops}"


def test_grouped_y_is_per_panel_by_default(tmp_path):
    """Assert grouped mode autoscales each group's y independently.

    Read-count groups differ in depth by orders of magnitude, so a global y top
    flattens the shallow panels.
    """
    import matplotlib.pyplot as plt

    df = pd.DataFrame({
        "sample": ["s1"] * 100 + ["s2"] * 10,
        "stage": ["T0"] * 110,
        "value": list(np.full(100, 4.0)) + list(np.full(10, 4.0)),
    })

    render_grouped_histogram_figure(
        df, tmp_path / "groupedfreey",
        value_column="value", row_key="sample", col_key="stage",
        bins=5,
    )

    tops = [ax.get_ylim()[1] for ax in plt.gcf().get_axes()]
    assert len(set(np.round(tops, 6))) == 2, f"Panels should autoscale independently: {tops}"


def test_grouped_shared_y_applies_uniform_ylim(tmp_path):
    """Assert share_y_range=True still forces one common y top across groups."""
    import matplotlib.pyplot as plt

    df = pd.DataFrame({
        "sample": ["s1"] * 100 + ["s2"] * 10,
        "stage": ["T0"] * 110,
        "value": list(np.full(100, 4.0)) + list(np.full(10, 4.0)),
    })

    render_grouped_histogram_figure(
        df, tmp_path / "groupedsharedy",
        value_column="value", row_key="sample", col_key="stage",
        bins=5, share_y_range=True,
    )

    tops = [ax.get_ylim()[1] for ax in plt.gcf().get_axes()]
    assert len(set(np.round(tops, 6))) == 1, f"Shared-y panels should share one y top: {tops}"


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
