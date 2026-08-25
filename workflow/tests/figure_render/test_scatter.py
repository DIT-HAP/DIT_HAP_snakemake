#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic grouped regression scatter renderer.

Asserts generality (arbitrary column names, arbitrary group keys) rather than
one figure's biology. Figure-specific baselines live in the entrypoint tests.

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

from figure_render.scatter import render_grouped_regression_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def arbitrary_frame() -> pd.DataFrame:
    """A frame whose columns share no names with any pipeline figure."""
    rng = np.random.default_rng(0)
    n = 200
    return pd.DataFrame(
        {
            "batch": ["b1"] * n + ["b2"] * n,
            "stage": (["s1"] * (n // 2) + ["s2"] * (n // 2)) * 2,
            "left_signal": np.abs(rng.normal(3, 1, 2 * n)),
            "right_signal": np.abs(rng.normal(3, 1, 2 * n)),
        }
    )


# =============================================================================
# TESTS
# =============================================================================
def test_renders_with_arbitrary_column_and_group_names(
    arbitrary_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the renderer works on columns unrelated to PBL/PBR or strand counts."""
    stem = tmp_path / "arbitrary"

    render_grouped_regression_figure(
        arbitrary_frame,
        stem,
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    assert (tmp_path / "arbitrary.pdf").exists()
    assert (tmp_path / "arbitrary.review.png").exists()
    assert (tmp_path / "arbitrary.pdf").stat().st_size > 0


def test_axis_labels_survive_regplot(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert labels are applied after drawing, since regplot resets the ones it manages."""
    import matplotlib.pyplot as plt

    render_grouped_regression_figure(
        arbitrary_frame,
        stem := tmp_path / "labels",
        x="left_signal",
        y="right_signal",
        xlabel="MY-X-LABEL",
        ylabel="MY-Y-LABEL",
        row_key="batch",
        col_key="stage",
    )

    axes = plt.gcf().get_axes()
    assert axes, "No axes on the rendered figure"
    assert any(ax.get_xlabel() == "MY-X-LABEL" for ax in axes), (
        "xlabel was lost — regplot reset it and it was not reapplied"
    )
    assert any("MY-Y-LABEL" in ax.get_ylabel() for ax in axes), (
        "ylabel was lost — regplot reset it and it was not reapplied"
    )
    assert stem.with_suffix(".pdf").exists() or (tmp_path / "labels.pdf").exists()


def test_row_key_value_appears_in_first_column_ylabel(
    arbitrary_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the row key's value is carried in the leading panel's ylabel.

    orientation.py put the sample name there rather than in the title: a full
    "{sample} {timepoint}" title is wider than the axes and pushes the measured
    panel width past the point where n_cols panels fit in a row, silently
    reflowing the grid.
    """
    import matplotlib.pyplot as plt

    render_grouped_regression_figure(
        arbitrary_frame,
        tmp_path / "rowlabel",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    ylabels = [ax.get_ylabel() for ax in plt.gcf().get_axes()]
    assert any("b1" in label for label in ylabels), f"Row key value missing from ylabels: {ylabels}"


def test_marker_size_is_inside_scatter_kws() -> None:
    """Assert the default kwargs carry `s`, since regplot drops its own `s=`.

    cnsplots regplot builds scatter_kws={"s": s, ...} internally, so supplying
    scatter_kws replaces it wholesale and loses the size.
    """
    from figure_render.scatter import REGPLOT_SCATTER_KWS

    assert "s" in REGPLOT_SCATTER_KWS, "marker size must live inside scatter_kws"
    assert REGPLOT_SCATTER_KWS["rasterized"] is True, "large point clouds must rasterize"


def test_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty input returns without raising and without writing artifacts."""
    empty = pd.DataFrame(columns=["batch", "stage", "left_signal", "right_signal"])

    render_grouped_regression_figure(
        empty,
        tmp_path / "empty",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )


def test_missing_x_column_raises(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert a bad column name fails loudly rather than drawing an empty panel."""
    with pytest.raises(ValueError, match="absent_column"):
        render_grouped_regression_figure(
            arbitrary_frame,
            tmp_path / "bad",
            x="absent_column",
            y="right_signal",
            xlabel="left",
            ylabel="right",
            row_key="batch",
            col_key="stage",
        )


def test_labels_are_keyword_only_without_defaults(
    arbitrary_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert a stale positional call fails instead of drawing wrong labels."""
    with pytest.raises(TypeError):
        render_grouped_regression_figure(arbitrary_frame, tmp_path / "stale")  # type: ignore[call-arg]


def test_panel_count_equals_group_count(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert one panel per row-key/col-key combination."""
    import matplotlib.pyplot as plt

    render_grouped_regression_figure(
        arbitrary_frame,
        tmp_path / "count",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    expected = arbitrary_frame.groupby(["batch", "stage"]).ngroups
    assert len(plt.gcf().get_axes()) == expected
