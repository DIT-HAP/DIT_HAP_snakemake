"""Generic composition figure rendering: a bar overview plus per-category donuts.

Supersedes ``coverage.py``, which showed gene coverage percentage per viability
category alongside covered/not-covered donuts. Column names and labels are
supplied by the caller.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
def _status_counts_frame(part: int, whole: int, part_label: str, whole_label: str) -> pd.DataFrame:
    """Expand two complementary counts into the long-format frame cns.donutplot expects."""
    return pd.DataFrame({"status": [part_label] * part + [whole_label] * whole})


@logger.catch(reraise=True)
def render_composition_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    category_column: str,
    percentage_column: str,
    part_column: str,
    whole_column: str,
    total_column: str,
    part_label: str,
    whole_label: str,
    xlabel: str,
    ylabel: str,
    title: str,
    donut_unit: str = "items",
) -> None:
    """Render a percentage bar chart per category plus one part/whole donut each."""
    logger.info("Rendering composition figure...")

    require_columns(
        df,
        [category_column, percentage_column, part_column, whole_column, total_column],
        context="composition input",
    )

    apply_house_style()

    if df.empty:
        logger.warning("No composition data to plot!")
        cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
        multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
        ax = multipanel.panel(label="A")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        save_dual(output_stem)
        return

    n_categories = len(df)
    logger.info(f"Creating figure with {1 + n_categories} panels ({n_categories} categories)...")

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    labels = panel_labels(1 + n_categories)

    # Panel A: percentage per category
    ax_bar = multipanel.panel(label=labels[0])
    cns.barplot(data=df, x=category_column, y=percentage_column, ax=ax_bar)
    ax_bar.set_ylabel(ylabel)
    ax_bar.set_xlabel(xlabel)
    ax_bar.set_title(title)
    ax_bar.set_ylim(0, 100)
    for index, row in df.reset_index(drop=True).iterrows():
        ax_bar.text(
            index, row[percentage_column], f"{row[percentage_column]:.1f}%",
            ha="center", va="bottom", fontsize=6,
        )

    # Panels B, C, ...: part vs whole per category, in the frame's row order
    for label, (_, row) in zip(labels[1:], df.iterrows(), strict=True):
        ax_donut = multipanel.panel(label=label)

        part, whole = int(row[part_column]), int(row[whole_column])
        if part + whole == 0:
            ax_donut.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax_donut.transAxes)
            ax_donut.set_title(str(row[category_column]))
            continue

        status_df = _status_counts_frame(part, whole, part_label, whole_label)
        cns.donutplot(data=status_df, x="status", order=[part_label, whole_label], ax=ax_donut)
        ax_donut.set_title(
            f"{row[category_column]}\n({part:,}/{int(row[total_column]):,} {donut_unit})"
        )

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
