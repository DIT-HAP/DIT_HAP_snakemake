#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic multi-panel scatter grid renderer.

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

from figure_render.scatter import ScatterPanel, render_scatter_grid_figure, render_scatter_panel


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def arbitrary_frame() -> pd.DataFrame:
    """A frame with generic metric names and an optional grouping column."""
    rng = np.random.default_rng(1)
    n = 300
    return pd.DataFrame(
        {
            "metric_a": rng.uniform(0, 10, n),
            "metric_b": rng.uniform(0, 10, n),
            "delta": rng.normal(0, 1, n),
            "group": rng.choice(["alpha", "beta"], n),
        }
    )


@pytest.fixture
def two_panels() -> list[ScatterPanel]:
    """Two panels exercising an identity reference line and a zero line."""
    return [
        ScatterPanel(
            x="metric_a", y="metric_b",
            xlabel="A", ylabel="B", title="A vs B",
            reference="identity",
        ),
        ScatterPanel(
            x="metric_a", y="delta",
            xlabel="A", ylabel="delta", title="Delta vs A",
            reference="zero",
        ),
    ]


# =============================================================================
# TESTS
# =============================================================================
def test_renders_one_panel_per_spec(
    arbitrary_frame: pd.DataFrame, two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert the panel count matches the number of specs supplied."""
    import matplotlib.pyplot as plt

    render_scatter_grid_figure(arbitrary_frame, tmp_path / "grid", panels=two_panels)

    assert len(plt.gcf().get_axes()) == 2
    assert (tmp_path / "grid.pdf").exists()
    assert (tmp_path / "grid.review.png").exists()


def test_scatter_kws_omit_edgecolor() -> None:
    """Assert the defaults omit edgecolor, which cns.scatterplot supplies itself.

    cns.scatterplot always passes edgecolor=None to seaborn internally, so
    supplying it here raises on the duplicate keyword.
    """
    from figure_render.scatter import SCATTERPLOT_KWS

    assert "edgecolor" not in SCATTERPLOT_KWS
    assert SCATTERPLOT_KWS["rasterized"] is True


def test_hue_is_optional(
    arbitrary_frame: pd.DataFrame, two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert rendering works both with and without a hue column."""
    render_scatter_grid_figure(arbitrary_frame, tmp_path / "nohue", panels=two_panels)
    assert (tmp_path / "nohue.pdf").exists()

    render_scatter_grid_figure(
        arbitrary_frame, tmp_path / "hue", panels=two_panels,
        hue="group", hue_order=["alpha", "beta"],
    )
    assert (tmp_path / "hue.pdf").exists()


def test_hue_order_filters_absent_levels(
    arbitrary_frame: pd.DataFrame, two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert a hue_order naming absent levels does not raise.

    density.py filtered hue_order against the observed values; the generic
    renderer must keep doing so or a project missing a viability level crashes.
    """
    render_scatter_grid_figure(
        arbitrary_frame, tmp_path / "extra_levels", panels=two_panels,
        hue="group", hue_order=["alpha", "beta", "never_present"],
    )

    assert (tmp_path / "extra_levels.pdf").exists()


def test_log_scale_panel_sets_symlog(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert a log_scale panel applies symlog to both axes, as density panel C did."""
    import matplotlib.pyplot as plt

    panels = [
        ScatterPanel(
            x="metric_a", y="metric_b",
            xlabel="A", ylabel="B", title="log",
            log_scale=True,
        )
    ]
    render_scatter_grid_figure(arbitrary_frame, tmp_path / "log", panels=panels)

    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "symlog"
    assert ax.get_yscale() == "symlog"


def test_empty_frame_renders_placeholder(
    two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert an empty frame still produces artifacts with placeholder panels.

    density.py rendered 'No valid data' panels rather than returning early, so
    the Snakemake output file always exists.
    """
    empty = pd.DataFrame(columns=["metric_a", "metric_b", "delta", "group"])

    render_scatter_grid_figure(empty, tmp_path / "empty", panels=two_panels)

    assert (tmp_path / "empty.pdf").exists(), "Snakemake requires the output to exist"


def test_missing_panel_column_raises(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert a panel naming an absent column fails loudly."""
    panels = [
        ScatterPanel(x="absent", y="metric_b", xlabel="A", ylabel="B", title="bad")
    ]

    with pytest.raises(ValueError, match="absent"):
        render_scatter_grid_figure(arbitrary_frame, tmp_path / "bad", panels=panels)


def test_render_scatter_panel_draws_on_given_axes(arbitrary_frame: pd.DataFrame) -> None:
    """Assert the single-axes primitive draws a scatter and reference line onto ax."""
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="A vs B", reference="identity",
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, arbitrary_frame, panel)

    assert ax.get_xlabel() == "A"
    assert ax.get_ylabel() == "B"
    assert ax.get_title() == "A vs B"
    assert len(ax.collections) > 0, "No scatter points drawn"
    assert len(ax.lines) > 0, "No identity reference line drawn"
    plt.close(ax.figure)


def test_render_scatter_panel_handles_empty_data() -> None:
    """Assert an empty frame draws a placeholder instead of raising."""
    import matplotlib.pyplot as plt

    panel = ScatterPanel(x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Empty")
    empty = pd.DataFrame(columns=["metric_a", "metric_b"])

    _, ax = plt.subplots()
    render_scatter_panel(ax, empty, panel)

    assert ax.get_title() == "Empty"
    assert len(ax.texts) > 0, "No 'No valid data' placeholder drawn"
    plt.close(ax.figure)
