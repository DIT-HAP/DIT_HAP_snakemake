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

from figure_render.scatter import ScatterPanel, render_grouped_regression_figure, render_scatter_panel


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


def test_axis_labels_survive_scatterplot(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert labels are applied after drawing, since scatterplot resets the ones it manages."""
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
        "xlabel was lost — scatterplot reset it and it was not reapplied"
    )
    assert any("MY-Y-LABEL" in ax.get_ylabel() for ax in axes), (
        "ylabel was lost — scatterplot reset it and it was not reapplied"
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
    """Assert the default kwargs carry `s` and rasterize, since these ship as **kws to cns.scatterplot."""
    from figure_render.scatter import REGRESSION_PANEL_SCATTER_KWS

    assert "s" in REGRESSION_PANEL_SCATTER_KWS, "marker size must live inside scatter_kws"
    assert REGRESSION_PANEL_SCATTER_KWS["rasterized"] is True, "large point clouds must rasterize"


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


def test_render_scatter_panel_identity_draws_dashed_line_not_fit(tmp_path: Path) -> None:
    """Assert an identity-reference panel draws a dashed guide line, not a solid fitted regression line.

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
    panel = ScatterPanel(
        x="left_signal", y="right_signal", xlabel="left", ylabel="right", title="Panel Title",
        reference="identity",
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, df, panel)

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


def test_render_scatter_panel_show_stats_true_annotates_n_r_p() -> None:
    """Assert show_stats=True draws an n/r/P text annotation."""
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(6)
    n = 100
    df = pd.DataFrame({
        "left_signal": np.abs(rng.normal(3, 1, n)),
        "right_signal": np.abs(rng.normal(3, 1, n)),
    })
    panel = ScatterPanel(
        x="left_signal", y="right_signal", xlabel="left", ylabel="right", title="Panel Title",
        show_stats=True,
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, df, panel)

    stats_texts = [t for t in ax.texts if "$n$=" in t.get_text()]
    assert len(stats_texts) == 1, f"Expected exactly one n/r/P annotation, got {ax.texts}"
    assert f"$n$={n}" in stats_texts[0].get_text()
    plt.close(ax.figure)


def test_render_scatter_panel_show_stats_default_false_omits_annotation() -> None:
    """Assert ScatterPanel's default show_stats=False omits the n/r/P text, preserving density.py's current panels."""
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(7)
    n = 100
    df = pd.DataFrame({
        "left_signal": np.abs(rng.normal(3, 1, n)),
        "right_signal": np.abs(rng.normal(3, 1, n)),
    })
    panel = ScatterPanel(
        x="left_signal", y="right_signal", xlabel="left", ylabel="right", title="Panel Title",
        reference="identity",
    )

    _, ax = plt.subplots()
    render_scatter_panel(ax, df, panel)

    assert not any("$n$=" in t.get_text() for t in ax.texts), "show_stats=False must omit the n/r/P annotation"
    assert len(ax.collections) > 0, "Scatter points must still be drawn"
    assert len(ax.lines) == 1, "Identity line must still be drawn"
    plt.close(ax.figure)
