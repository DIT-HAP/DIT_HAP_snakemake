"""Generic observed-vs-fitted curve grid rendering.

Supersedes ``curve_fitting.py``. The curve model is injected by the caller, so
this module holds no specific model; the pipeline passes
``depletion.curve_model.sigmoid_function``.

Value columns come from the caller rather than a name whitelist. The old
renderer selected them with a literal ['YES0'..'YES4'] test, which silently
picked the wrong subset for any project with a different timepoint count and
raised "x and y must be the same size" from inside matplotlib.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Callable, Sequence
from pathlib import Path

import cnsplots as cns
import numpy as np
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

from ._layout import panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# House palette (Cell) indices 0 and 1. The old renderer hardcoded matplotlib's
# '#1f77b4' / '#ff7f0e', bypassing the palette apply_house_style() installs.
OBSERVED_COLOR = "#c84c3a"
FITTED_COLOR = "#2f7e8f"

# Preserved from the legacy renderer so panels stay comparable across projects.
# Callers may override, or pass None to autoscale.
DEFAULT_YLIM: tuple[float, float] = (-1.5, 8.5)

PANEL_WIDTH_PX = 180
PANEL_HEIGHT_PX = 150

# Points used to draw the fitted curve smoothly between the observed x values.
FITTED_CURVE_RESOLUTION = 100


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_fitted_curves_figure(
    observed: pd.DataFrame,
    output_stem: Path,
    *,
    x_values: Sequence[float],
    value_columns: Sequence[str],
    model: Callable[..., np.ndarray],
    model_params: Sequence[str],
    annotations: Sequence[str] = (),
    xlabel: str,
    ylabel: str,
    ylim: tuple[float, float] | None = DEFAULT_YLIM,
    n_cols: int = 4,
) -> None:
    """Render one panel per row: observed points plus the model curve from its parameters."""
    logger.info("Rendering fitted curves figure...")

    # Checked before the empty guard so a mismatch is reported even for an empty
    # frame, and before any column access so the message names the real problem.
    if len(x_values) != len(value_columns):
        raise ValueError(
            f"x_values and value_columns must have equal length: "
            f"got {len(x_values)} x values and {len(value_columns)} value columns "
            f"({list(value_columns)})"
        )

    require_columns(
        observed,
        [*value_columns, *model_params, *annotations],
        context="fitted curves input",
    )

    if observed.empty:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    n_panels = len(observed)
    n_rows = (n_panels + n_cols - 1) // n_cols
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels} panels)...")

    multipanel = cns.multipanel(max_width=PANEL_WIDTH_PX * n_cols)
    labels = panel_labels(n_panels)

    x_array = np.asarray(x_values, dtype=float)
    x_smooth = np.linspace(x_array.min(), x_array.max(), FITTED_CURVE_RESOLUTION)

    for label, (row_index, row) in zip(labels, observed.iterrows(), strict=True):
        ax = multipanel.panel(label=label, width=PANEL_WIDTH_PX, height=PANEL_HEIGHT_PX)

        y_observed = row[list(value_columns)].to_numpy(dtype=float)
        params = [float(row[name]) for name in model_params]

        ax.scatter(
            x_array,
            y_observed,
            s=30,
            color=OBSERVED_COLOR,
            alpha=0.8,
            edgecolor="white",
            linewidths=0.5,
        )
        ax.plot(x_smooth, model(x_smooth, *params), color=FITTED_COLOR, linewidth=2, label="Fitted")

        if annotations:
            text = "\n".join(f"{name}={row[name]:.3f}" for name in annotations)
            ax.text(0.05, 0.95, text, transform=ax.transAxes, verticalalignment="top", fontsize=8)

        title = row_index if isinstance(row_index, str) else " ".join(map(str, np.atleast_1d(row_index)))
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)

        if ylim is not None:
            ax.set_ylim(*ylim)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
