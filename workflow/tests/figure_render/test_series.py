#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic multi-series scatter renderer.

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

from figure_render.series import Series, render_series_scatter_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def series_frame() -> pd.DataFrame:
    """A frame with one shared x column and three y series."""
    rng = np.random.default_rng(4)
    n = 400
    return pd.DataFrame(
        {
            "shared_x": np.abs(rng.normal(100, 30, n)),
            "series_one": np.abs(rng.normal(0.5, 0.2, n)),
            "series_two": np.abs(rng.normal(0.4, 0.2, n)),
            "series_three": np.abs(rng.normal(0.3, 0.1, n)),
        }
    )


@pytest.fixture
def three_series() -> list[Series]:
    """Three series with explicit labels and colours."""
    return [
        Series(column="series_one", label="First"),
        Series(column="series_two", label="Second"),
        Series(column="series_three", label="Third"),
    ]


# =============================================================================
# TESTS
# =============================================================================
def test_renders_all_series_on_one_axes(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert every series lands on a single shared axes."""
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "series",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="Three series",
    )

    axes = plt.gcf().get_axes()
    assert len(axes) == 1
    assert len(axes[0].collections) == 3
    assert (tmp_path / "series.pdf").exists()


def test_legend_labels_match_series(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert legend text comes from the Series labels."""
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "legend",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )

    legend = plt.gcf().get_axes()[0].get_legend()
    assert legend is not None
    assert [t.get_text() for t in legend.get_texts()] == ["First", "Second", "Third"]


def test_legend_markers_are_scaled_and_opaque(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert legend handles are made visible despite tiny, low-alpha scatter points.

    dispersions.py used markerscale=6 and reset handle alpha to 1.0, because the
    handles otherwise inherit s=1.0 and alpha=0.12 and are invisible.
    """
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "handles",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )

    legend = plt.gcf().get_axes()[0].get_legend()
    assert all(handle.get_alpha() == 1.0 for handle in legend.legend_handles)


def test_log_scales_are_toggleable(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert log_x / log_y control the axis scales."""
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "loglog",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"

    render_series_scatter_figure(
        series_frame, tmp_path / "linear",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
        log_x=False, log_y=False,
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "linear" and ax.get_yscale() == "linear"


def test_empty_frame_renders_placeholder(three_series: list[Series], tmp_path: Path) -> None:
    """Assert an empty frame still writes artifacts with a placeholder panel.

    dispersions.py rendered the placeholder and saved, so the Snakemake output
    always exists.
    """
    empty = pd.DataFrame(columns=["shared_x", "series_one", "series_two", "series_three"])

    render_series_scatter_figure(
        empty, tmp_path / "empty",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )

    assert (tmp_path / "empty.pdf").exists(), "Snakemake requires the output to exist"


def test_missing_series_column_raises(series_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert an absent series column is reported by name rather than partially plotted."""
    bad = [Series(column="absent_series", label="X")]

    with pytest.raises(ValueError, match="absent_series"):
        render_series_scatter_figure(
            series_frame, tmp_path / "bad",
            x="shared_x", series=bad,
            xlabel="x", ylabel="y", title="t",
        )
