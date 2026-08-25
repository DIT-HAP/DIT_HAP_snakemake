#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic MA plot renderer.

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

from figure_render.ma import Orientation, render_ma_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def two_panels() -> list[tuple[str, pd.Series, pd.Series]]:
    """Two (title, abundance, effect) triples with positive abundances."""
    rng = np.random.default_rng(5)
    return [
        (
            f"stage{i}",
            pd.Series(np.abs(rng.normal(500, 200, 300)) + 1),
            pd.Series(rng.normal(0, 1.5, 300)),
        )
        for i in (1, 2)
    ]


# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize("orientation", [Orientation.VERTICAL, Orientation.HORIZONTAL])
def test_renders_in_both_orientations(
    two_panels: list[tuple[str, pd.Series, pd.Series]],
    orientation: Orientation,
    tmp_path: Path,
) -> None:
    """Assert both orientations produce artifacts."""
    stem = tmp_path / f"ma_{orientation.value}"

    render_ma_figure(
        two_panels, stem,
        abundance_label="mean of normalized counts",
        effect_label="log2 fold change",
        title_prefix="MA plot",
        orientation=orientation,
    )

    assert stem.with_name(f"{stem.name}.pdf").exists()


def test_orientation_controls_axis_assignment_and_log_axis(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert abundance goes on the log axis in whichever direction it is drawn.

    Vertical put effect on x and log-scaled y; horizontal put abundance on x and
    log-scaled x. Abundance spans orders of magnitude; the effect is symmetric.
    """
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "vert",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        orientation=Orientation.VERTICAL,
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_yscale() == "log"
    assert ax.get_xlabel() == "effect"
    assert ax.get_ylabel() == "abundance"

    render_ma_figure(
        two_panels, tmp_path / "horiz",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        orientation=Orientation.HORIZONTAL,
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log"
    assert ax.get_xlabel() == "abundance"
    assert ax.get_ylabel() == "effect"


def test_stack_is_independent_of_orientation(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert horizontal axes can be stacked vertically.

    ma_plot_replicates.py used max_width=510 (a vertical stack) while assigning
    axes horizontally: scatter(abundance, effect) with a log x-axis and axhline.
    A single enum cannot express that combination.
    """
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "horiz_stacked",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        orientation=Orientation.HORIZONTAL, stack=True,
    )

    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log", "Axis assignment must stay horizontal"
    assert ax.get_xlabel() == "abundance"
    assert (tmp_path / "horiz_stacked.pdf").exists()


def test_stack_defaults_from_orientation(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert omitting `stack` stacks for VERTICAL and rows for HORIZONTAL."""
    render_ma_figure(
        two_panels, tmp_path / "def_vert",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.VERTICAL,
    )
    assert (tmp_path / "def_vert.pdf").exists()

    render_ma_figure(
        two_panels, tmp_path / "def_horiz",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.HORIZONTAL,
    )
    assert (tmp_path / "def_horiz.pdf").exists()


def test_per_point_colors_are_accepted(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert a per-point colour Series works, as the replicate branch needs."""
    colors = [
        pd.Series(np.where(np.arange(len(abundance)) % 2 == 0, "darkred", "gray"))
        for _, abundance, _ in two_panels
    ]

    render_ma_figure(
        two_panels, tmp_path / "colored",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        point_colors=colors,
    )

    assert (tmp_path / "colored.pdf").exists()


def test_single_color_string_is_accepted(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert a scalar colour per panel works, as the no-replicate branch needs."""
    render_ma_figure(
        two_panels, tmp_path / "single",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        point_colors=["gray", "gray"],
    )

    assert (tmp_path / "single.pdf").exists()


def test_null_effect_reference_line_orientation(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert the reference line is vertical when effect is on x, horizontal otherwise."""
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "refv",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.VERTICAL,
    )
    assert len(plt.gcf().get_axes()[0].lines) >= 1

    render_ma_figure(
        two_panels, tmp_path / "refh",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.HORIZONTAL,
    )
    assert len(plt.gcf().get_axes()[0].lines) >= 1


def test_titles_use_prefix_and_panel_name(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert each panel's title combines the prefix with the panel name."""
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "titles",
        abundance_label="a", effect_label="e", title_prefix="MA plot",
    )

    titles = [ax.get_title() for ax in plt.gcf().get_axes()]
    assert "MA plot - stage1" in titles
    assert "MA plot - stage2" in titles


def test_panel_count_matches_input(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert one panel per supplied triple."""
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "count",
        abundance_label="a", effect_label="e", title_prefix="MA",
    )

    assert len(plt.gcf().get_axes()) == 2


def test_empty_panels_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty panel sequence returns without raising."""
    render_ma_figure(
        [], tmp_path / "empty",
        abundance_label="a", effect_label="e", title_prefix="MA",
    )


def test_mismatched_color_count_raises(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert a point_colors length mismatch is refused rather than silently zipped short."""
    with pytest.raises(ValueError, match="point_colors"):
        render_ma_figure(
            two_panels, tmp_path / "badcolors",
            abundance_label="a", effect_label="e", title_prefix="MA",
            point_colors=["gray"],
        )
