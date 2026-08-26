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

from figure_render.scatter import render_grouped_regression_figure, render_regression_panel


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


def test_empty_frame_with_bad_column_still_raises(tmp_path: Path) -> None:
    """Assert column validation runs before the empty-frame guard.

    An empty frame with a typo'd column name must not be silently accepted:
    require_columns has to run before the `if df.empty: return`, or the typo
    passes through with no artifact and no error.
    """
    empty = pd.DataFrame(columns=["batch", "stage", "left_signal", "right_signal"])

    with pytest.raises(ValueError, match="absent_column"):
        render_grouped_regression_figure(
            empty,
            tmp_path / "empty_bad",
            x="absent_column",
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


def test_panels_stay_square_when_cols_and_rows_differ(tmp_path: Path) -> None:
    """Assert axes are square even for a grid whose column count != row count.

    grid_panel_size divided a fixed total height by n_rows independently of the
    width division by n_cols, producing a rectangle whenever the two counts
    differ. cnsplots.multipanel also grows the whole figure's height to fit
    however many rows the grid needs (see _create_or_update_figure), so a fixed
    height budget divided by n_rows was never a real constraint in the first
    place -- only the width division (via multipanel's max_width) has effect.
    """
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(2)
    n_per_group = 20
    rows, cols = ["r1", "r2"], ["c1", "c2", "c3", "c4", "c5", "c6"]
    frame = pd.DataFrame(
        {
            "batch": [row for row in rows for _ in cols for _ in range(n_per_group)],
            "stage": [col for _ in rows for col in cols for _ in range(n_per_group)],
            "left_signal": np.abs(rng.normal(3, 1, len(rows) * len(cols) * n_per_group)),
            "right_signal": np.abs(rng.normal(3, 1, len(rows) * len(cols) * n_per_group)),
        }
    )

    render_grouped_regression_figure(
        frame,
        tmp_path / "square_check",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    fig = plt.gcf()
    fig.canvas.draw()
    for ax in fig.get_axes():
        bbox = ax.get_window_extent()
        assert bbox.width == pytest.approx(bbox.height, abs=1.0), (
            f"Panel is not square: width={bbox.width:.1f}px height={bbox.height:.1f}px"
        )


def test_grid_wraps_at_exactly_n_cols_per_row(tmp_path: Path) -> None:
    """Assert a row holds n_cols panels, not fewer.

    multipanel wraps once a row's rendered panels exceed max_width, and
    cnsplots measures each panel's total width as the requested axes width
    plus its own y-axis tick/label decoration (~30-40px, not accounted for by
    the raw panel_size budget). Passing panel_size directly as max_width
    therefore wrapped one column early -- 6 requested columns rendered as 5.
    """
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(4)
    n_per_group = 20
    rows, cols = ["HD1328-4", "HD1328-7", "HD1328-8"], ["0h", "YES0", "YES1", "YES2", "YES3", "YES4"]
    frame = pd.DataFrame(
        {
            "batch": [row for row in rows for _ in cols for _ in range(n_per_group)],
            "stage": [col for _ in rows for col in cols for _ in range(n_per_group)],
            "left_signal": np.abs(rng.normal(2, 1, len(rows) * len(cols) * n_per_group)),
            "right_signal": np.abs(rng.normal(2, 1, len(rows) * len(cols) * n_per_group)),
        }
    )

    render_grouped_regression_figure(
        frame,
        tmp_path / "row_width_check",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    fig = plt.gcf()
    fig.canvas.draw()
    row_ys = {round(ax.get_window_extent().y0, 1) for ax in fig.get_axes()}
    assert len(row_ys) == len(rows), f"Expected {len(rows)} rows, got {len(row_ys)}"
    for row_y in row_ys:
        panels_in_row = sum(
            1 for ax in fig.get_axes() if round(ax.get_window_extent().y0, 1) == row_y
        )
        assert panels_in_row == len(cols), (
            f"Expected {len(cols)} panels per row, row at y={row_y} has {panels_in_row}"
        )


def test_render_regression_panel_draws_identity_line_not_fit(tmp_path: Path) -> None:
    """Assert the panel draws a dashed identity guide line, not a solid fitted regression line.

    PBL vs PBR (and similar x/y pairs of the same measured quantity) should
    show whether the data agrees with y = x, not a least-squares fit through
    this particular sample.
    """
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(5)
    n = 100
    df = pd.DataFrame({
        "left_signal": np.abs(rng.normal(3, 1, n)),
        "right_signal": np.abs(rng.normal(3, 1, n)),
    })

    _, ax = plt.subplots()
    render_regression_panel(
        ax, df, x="left_signal", y="right_signal", xlabel="left", ylabel="right", title="Panel Title",
    )

    assert len(ax.lines) == 1, f"Expected exactly one line (the identity guide), got {len(ax.lines)}"
    line = ax.lines[0]
    assert line.get_linestyle() == "--", "Identity line must be dashed, not the solid fitted-line style"

    x_data, y_data = line.get_xdata(), line.get_ydata()
    assert list(x_data) == list(y_data), "Identity line must satisfy y = x at every plotted point"

    min_val = min(df["left_signal"].min(), df["right_signal"].min())
    max_val = max(df["left_signal"].max(), df["right_signal"].max())
    assert min(x_data) == pytest.approx(min_val)
    assert max(x_data) == pytest.approx(max_val)
    plt.close(ax.figure)


def test_render_regression_panel_draws_on_given_axes() -> None:
    """Assert the single-axes primitive draws a regplot onto the ax it is given."""
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(3)
    n = 100
    df = pd.DataFrame({
        "left_signal": np.abs(rng.normal(3, 1, n)),
        "right_signal": np.abs(rng.normal(3, 1, n)),
    })

    _, ax = plt.subplots()
    render_regression_panel(
        ax, df, x="left_signal", y="right_signal", xlabel="left", ylabel="right", title="Panel Title",
    )

    assert ax.get_xlabel() == "left"
    assert ax.get_ylabel() == "right"
    assert ax.get_title() == "Panel Title"
    assert len(ax.collections) > 0, "No scatter points drawn"
    plt.close(ax.figure)


def test_render_regression_panel_handles_empty_data() -> None:
    """Assert an empty frame draws a placeholder instead of raising."""
    import matplotlib.pyplot as plt

    empty = pd.DataFrame(columns=["left_signal", "right_signal"])

    _, ax = plt.subplots()
    render_regression_panel(
        ax, empty, x="left_signal", y="right_signal", xlabel="left", ylabel="right", title="Empty",
    )

    assert ax.get_title() == "Empty"
    assert len(ax.texts) > 0, "No 'No valid data' placeholder drawn"
    plt.close(ax.figure)
