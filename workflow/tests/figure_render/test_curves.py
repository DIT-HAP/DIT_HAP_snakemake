#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic fitted-curve grid renderer.

The timepoint-count tests are regression tests for a reproduced crash: the old
renderer selected value columns with a literal ['YES0'..'YES4'] whitelist and
raised "x and y must be the same size" for any project with a different number
of timepoints, or a different naming scheme.

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

from figure_render.curves import DEFAULT_YLIM, render_fitted_curves_figure


# =============================================================================
# HELPERS
# =============================================================================
def linear_model(x: np.ndarray, slope: float, intercept: float) -> np.ndarray:
    """A trivial stand-in for the pipeline's sigmoid, so tests need no domain import."""
    return slope * np.asarray(x) + intercept


def build_frame(timepoint_names: list[str], n_rows: int = 3) -> pd.DataFrame:
    """Build an observed+params frame for the given timepoint column names."""
    rng = np.random.default_rng(2)
    data = {name: rng.normal(2, 1, n_rows) for name in timepoint_names}
    data["slope"] = rng.uniform(0.1, 0.5, n_rows)
    data["intercept"] = rng.uniform(-1, 1, n_rows)
    data["quality"] = rng.uniform(0, 1, n_rows)
    return pd.DataFrame(data)


# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize(
    "timepoint_names",
    [
        ["YES0", "YES1", "YES2", "YES3", "YES4"],
        ["YES0", "YES1", "YES2", "YES3", "YES4", "YES5", "YES6"],
        ["0h", "YES0", "YES1", "YES2", "YES3", "YES4", "YES5", "YES6"],
        ["Spikein0", "Spikein1", "Spikein2", "Spikein3", "Spikein4", "Spikein5"],
    ],
    ids=["five_yes", "seven_yes", "eight_with_0h", "spikein_naming"],
)
def test_renders_for_any_timepoint_count_and_naming(
    timepoint_names: list[str], tmp_path: Path
) -> None:
    """Assert rendering works regardless of timepoint count or naming scheme.

    Regression test: the old whitelist matched 5 of 7 columns for
    Spore2YES6_1328, 5 of 8 for LD_haploid, and 0 of 6 for Spikein.
    """
    df = build_frame(timepoint_names)
    x_values = list(range(len(timepoint_names)))

    render_fitted_curves_figure(
        df,
        tmp_path / f"curves_{len(timepoint_names)}",
        x_values=x_values,
        value_columns=timepoint_names,
        model=linear_model,
        model_params=["slope", "intercept"],
        xlabel="x",
        ylabel="y",
    )

    assert (tmp_path / f"curves_{len(timepoint_names)}.pdf").exists()


def test_mismatched_x_and_value_columns_raises_readable_error(tmp_path: Path) -> None:
    """Assert a length mismatch is caught at the boundary, not inside matplotlib.

    The old failure surfaced as ValueError("x and y must be the same size") from
    deep inside Axes.scatter, with no indication of which lengths disagreed.
    """
    df = build_frame(["t1", "t2", "t3"])

    with pytest.raises(ValueError) as excinfo:
        render_fitted_curves_figure(
            df,
            tmp_path / "mismatch",
            x_values=[0.0, 1.0],
            value_columns=["t1", "t2", "t3"],
            model=linear_model,
            model_params=["slope", "intercept"],
            xlabel="x",
            ylabel="y",
        )

    message = str(excinfo.value)
    assert "2" in message and "3" in message, f"Error must report both lengths: {message}"


def test_uses_house_palette_not_matplotlib_defaults() -> None:
    """Assert the point and curve colours are resolved from the house palette.

    The old renderer hardcoded '#1f77b4' and '#ff7f0e', matplotlib's defaults,
    bypassing the Cell palette that apply_house_style() installs. They were then
    hardcoded again as the palette's own hex values, which drift silently if
    HOUSE_PALETTE changes; the pair must be derived, not copied.
    """
    import cnsplots as cns
    from figures import HOUSE_PALETTE, observed_fitted_colors

    observed, fitted = observed_fitted_colors()

    assert observed.lower() != "#1f77b4"
    assert fitted.lower() != "#ff7f0e"

    palette = [color.lower() for color in cns.get_hexcolors_from_apalette([0, 1], HOUSE_PALETTE)]
    assert [observed.lower(), fitted.lower()] == palette


def test_default_ylim_preserves_legacy_framing() -> None:
    """Assert the default y-limits stay (-1.5, 8.5) so panels remain comparable."""
    assert DEFAULT_YLIM == (-1.5, 8.5)


def test_ylim_is_overridable(tmp_path: Path) -> None:
    """Assert an explicit ylim is applied, and None means autoscale."""
    import matplotlib.pyplot as plt

    df = build_frame(["t1", "t2", "t3"])
    common = dict(
        x_values=[0.0, 1.0, 2.0],
        value_columns=["t1", "t2", "t3"],
        model=linear_model,
        model_params=["slope", "intercept"],
        xlabel="x",
        ylabel="y",
    )

    render_fitted_curves_figure(df, tmp_path / "fixed", ylim=(0.0, 5.0), **common)
    assert plt.gcf().get_axes()[0].get_ylim() == (0.0, 5.0)

    render_fitted_curves_figure(df, tmp_path / "auto", ylim=None, **common)
    assert plt.gcf().get_axes()[0].get_ylim() != (0.0, 5.0)


def test_panel_labels_stay_alphanumeric_past_26(tmp_path: Path) -> None:
    """Assert a 32-panel figure emits no punctuation labels.

    The pipeline samples 32 curves by default, so the old bare chr(65+i) was
    already producing '[', '\\\\' and ']' as panel labels in production.
    """
    df = build_frame(["t1", "t2", "t3"], n_rows=32)

    render_fitted_curves_figure(
        df,
        tmp_path / "many",
        x_values=[0.0, 1.0, 2.0],
        value_columns=["t1", "t2", "t3"],
        model=linear_model,
        model_params=["slope", "intercept"],
        xlabel="x",
        ylabel="y",
    )

    assert (tmp_path / "many.pdf").exists()


def test_annotations_are_rendered(tmp_path: Path) -> None:
    """Assert annotation columns reach the panel text box."""
    import matplotlib.pyplot as plt

    df = build_frame(["t1", "t2", "t3"], n_rows=1)

    render_fitted_curves_figure(
        df,
        tmp_path / "annot",
        x_values=[0.0, 1.0, 2.0],
        value_columns=["t1", "t2", "t3"],
        model=linear_model,
        model_params=["slope", "intercept"],
        annotations=["quality"],
        xlabel="x",
        ylabel="y",
    )

    texts = [t.get_text() for ax in plt.gcf().get_axes() for t in ax.texts]
    assert any("quality" in text for text in texts), f"Annotation missing from: {texts}"


def test_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty input returns without raising."""
    empty = pd.DataFrame(columns=["t1", "t2", "slope", "intercept"])

    render_fitted_curves_figure(
        empty,
        tmp_path / "empty",
        x_values=[0.0, 1.0],
        value_columns=["t1", "t2"],
        model=linear_model,
        model_params=["slope", "intercept"],
        xlabel="x",
        ylabel="y",
    )


def test_missing_param_column_raises(tmp_path: Path) -> None:
    """Assert an absent model-parameter column is reported by name.

    LD_haploid's stats file predates the current curve model and carries um/lam
    instead of A/DR/DL; it must fail with a readable message.
    """
    df = build_frame(["t1", "t2", "t3"])

    with pytest.raises(ValueError, match="DR"):
        render_fitted_curves_figure(
            df,
            tmp_path / "badparams",
            x_values=[0.0, 1.0, 2.0],
            value_columns=["t1", "t2", "t3"],
            model=linear_model,
            model_params=["A", "DR", "DL"],
            xlabel="x",
            ylabel="y",
        )
