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


def test_density_colours_points_individually(arbitrary_frame: pd.DataFrame) -> None:
    """Assert density=True maps a per-point value array through a colormap.

    Asserted on the scalar array rather than get_facecolors(): matplotlib
    resolves colormapped facecolors lazily at draw time, so before a draw they
    still read as one flat colour.
    """
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Density", density=True,
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, arbitrary_frame, panel)

    collection = ax.collections[0]
    values = collection.get_array()
    assert values is not None, "No scalar array set; density was not applied"
    assert len(values) == len(arbitrary_frame)
    assert values.min() < values.max(), "Density is uniform across all points"

    collection.update_scalarmappable()
    colors = collection.get_facecolors()
    assert len(np.unique(colors, axis=0)) > 1, "Colormap did not produce distinct point colours"
    plt.close(ax.figure)


def test_density_adds_exactly_one_inset_colorbar(arbitrary_frame: pd.DataFrame) -> None:
    """Assert the colorbar is an inset child axes.

    An inset contributes nothing to the panel width/height cnsplots measures
    for layout, so it cannot reflow the grid. fig.colorbar(ax=...) would take
    space from the panel itself and fight square_panel_size's wrap arithmetic.
    """
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Density", density=True,
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, arbitrary_frame, panel)

    assert len(ax.child_axes) == 1, f"Expected one inset colorbar, got {len(ax.child_axes)}"
    assert ax.child_axes[0].get_ylabel() == "Density"
    assert [text.get_text() for text in ax.child_axes[0].get_yticklabels()] == ["low", "high"]
    plt.close(ax.figure)


def test_density_colorbar_labels_have_no_tick_marks(arbitrary_frame: pd.DataFrame) -> None:
    """Assert the low/high labels are kept but their tick marks are not drawn."""
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Density", density=True,
    )

    fig, ax = plt.subplots(dpi=200)
    render_scatter_panel(ax, arbitrary_frame, panel)
    fig.canvas.draw()

    cax = ax.child_axes[0]
    assert [label.get_text() for label in cax.get_yticklabels()] == ["low", "high"]
    for tick in cax.yaxis.get_major_ticks():
        assert tick.tick1line.get_markersize() == 0, "Tick mark is still drawn"
        assert tick.tick2line.get_markersize() == 0, "Tick mark is still drawn"
    plt.close(fig)


def test_density_colorbar_labels_fall_inside_the_bar(arbitrary_frame: pd.DataFrame) -> None:
    """Assert both end ticks lie within the norm range and their labels inside the bar.

    A bare Normalize() autoscales to the raw density range, leaving ticks
    requested at 0 and 1 out of range: matplotlib then clips them onto the
    edges, dropping one label entirely and offsetting the other. Asserting the
    label text alone does not catch that.
    """
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Density", density=True,
    )

    fig, ax = plt.subplots(dpi=200)
    render_scatter_panel(ax, arbitrary_frame, panel)
    fig.canvas.draw()

    norm = ax.collections[0].norm
    cax = ax.child_axes[0]

    ticks = cax.get_yticks()
    assert all(norm.vmin <= tick <= norm.vmax for tick in ticks), (
        f"Ticks {list(ticks)} fall outside the norm range {norm.vmin}-{norm.vmax}"
    )

    bar = cax.get_window_extent()
    for label in cax.get_yticklabels():
        box = label.get_window_extent()
        assert box.y0 >= bar.y0 - 1 and box.y1 <= bar.y1 + 1, (
            f"Label {label.get_text()!r} at y {box.y0:.1f}-{box.y1:.1f} "
            f"spills outside the bar at {bar.y0:.1f}-{bar.y1:.1f}"
        )
    plt.close(fig)


def test_density_gradient_spans_the_full_colormap(arbitrary_frame: pd.DataFrame) -> None:
    """Assert the per-panel rescaling maps to [0, 1], so the gradient is not a partial slice."""
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Density", density=True,
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, arbitrary_frame, panel)

    collection = ax.collections[0]
    assert [collection.norm.vmin, collection.norm.vmax] == [0, 1]
    values = collection.get_array()
    assert values is not None
    assert values.min() == pytest.approx(0) and values.max() == pytest.approx(1)
    plt.close(ax.figure)


def test_density_defaults_off(arbitrary_frame: pd.DataFrame) -> None:
    """Assert the default leaves the existing flat-colour path and adds no colorbar."""
    import matplotlib.pyplot as plt

    panel = ScatterPanel(x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Plain")

    _, ax = plt.subplots()
    render_scatter_panel(ax, arbitrary_frame, panel)

    assert ax.child_axes == [], "Default panel gained a colorbar"
    assert ax.collections[0].get_array() is None, "Default panel got a per-point value array"
    plt.close(ax.figure)


def test_hue_wins_over_density(arbitrary_frame: pd.DataFrame) -> None:
    """Assert hue takes the marker colour and density is skipped, since both claim it."""
    import matplotlib.pyplot as plt

    panel = ScatterPanel(
        x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Both", density=True,
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, arbitrary_frame, panel, hue="group", hue_order=["alpha", "beta"])

    assert ax.child_axes == [], "Density colorbar drawn despite hue being active"
    assert ax.get_legend() is not None, "Hue legend missing"
    plt.close(ax.figure)


def test_density_falls_back_on_degenerate_column(tmp_path: Path) -> None:
    """Assert a zero-variance column still renders rather than raising inside gaussian_kde."""
    constant = pd.DataFrame({"metric_a": np.full(50, 2.0), "metric_b": np.arange(50.0)})
    panels = [
        ScatterPanel(
            x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Constant", density=True,
        )
    ]

    render_scatter_grid_figure(constant, tmp_path / "degenerate", panels=panels)

    assert (tmp_path / "degenerate.pdf").exists(), "Snakemake requires the output to exist"


def test_density_reaches_grid_figure_panels(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert the ScatterPanel field is honoured through the grid orchestrator, not just directly."""
    import matplotlib.pyplot as plt

    panels = [
        ScatterPanel(
            x="metric_a", y="metric_b", xlabel="A", ylabel="B", title="Density", density=True,
        )
    ]

    render_scatter_grid_figure(arbitrary_frame, tmp_path / "grid_density", panels=panels)

    ax = plt.gcf().get_axes()[0]
    assert len(ax.child_axes) == 1, "Density colorbar missing from the grid-rendered panel"
