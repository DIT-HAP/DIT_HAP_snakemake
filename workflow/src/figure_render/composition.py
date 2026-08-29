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
from .donut import draw_donut_panel

# =============================================================================
# CONSTANTS
# =============================================================================
# Bar overview panel: its total width (panel label + pad_left=45 + category
# names + axes + margin) must match the two-per-row donut grid beneath it so
# the figure reads as one aligned column. Measured from the rendered layout:
# the donut row spans ~247 px (labels + 110 + margin 8 twice); the bar's left
# decorations consume ~104 px, leaving ~136 px for the axes. multipanel's
# max_width is tightened from the journal width to ~250 so the bar keeps its
# own row (104+136+8+120 > 250) while the two donuts stay side by side
# (120 + 124 < 250); canvas height grows automatically.
BAR_PANEL_WIDTH_PX = 136
BAR_PANEL_HEIGHT_PX = 110
BAR_PANEL_PAD_LEFT_PX = 45  # gap for the category names, as in the gallery
BAR_PANEL_MARGIN_RIGHT_PX = 8  # match the donut row's right margin
MULTIPANEL_MAX_WIDTH_PX = 250
DONUT_PANEL_WIDTH_PX = 110
DONUT_PANEL_HEIGHT_PX = 110
DONUTS_PER_ROW = 2
DONUT_PANEL_MARGIN_RIGHT_PX = 8

PERCENTAGE_LABEL_FONTSIZE = 6
# Donut-hole text: 6 pt Arial runs ~2 px/char; the hole is ~53 px across so
# ~15 chars per line fits with headroom. Category names longer than that are
# hyphen-wrapped before the n-count line is appended.
CENTER_LINE_MAX_CHARS = 15


def _wrap_center_text(center_text: str) -> str:
    """Hyphen-wrap ``center_text`` so no line exceeds CENTER_LINE_MAX_CHARS."""
    if len(center_text) <= CENTER_LINE_MAX_CHARS:
        return center_text
    lines: list[str] = []
    current = ""
    for part in center_text.split("-"):
        candidate = f"{current}-{part}" if current else part
        if len(candidate) > CENTER_LINE_MAX_CHARS and current:
            lines.append(current)
            current = part
        else:
            current = candidate
    if current:
        lines.append(current)
    return "\n".join(lines)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_composition_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    category_column: str,
    percentage_column: str,
    part_column: str,
    whole_column: str,
    part_label: str,
    whole_label: str,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    """Render a horizontal percentage bar chart plus one part/whole donut each.

    Categories run down the y axis so long names stay readable; percentage
    labels sit at the bar ends. Each donut carries its category name in the
    hole instead of a title above, saving the vertical title band.
    """
    logger.info("Rendering composition figure...")

    require_columns(
        df,
        [category_column, percentage_column, part_column, whole_column],
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
    multipanel = cns.multipanel(max_width=MULTIPANEL_MAX_WIDTH_PX)
    labels = panel_labels(1 + n_categories)

    # Panel A: horizontal percentage bars, categories down the y axis so long
    # names stay readable; percentages annotate the bar ends.
    ax_bar = multipanel.panel(
        label=labels[0],
        width=BAR_PANEL_WIDTH_PX,
        height=BAR_PANEL_HEIGHT_PX,
        pad_left=BAR_PANEL_PAD_LEFT_PX,
        margin_right=BAR_PANEL_MARGIN_RIGHT_PX,
    )
    cns.barplot(data=df, y=category_column, x=percentage_column, ax=ax_bar)
    ax_bar.set_ylabel(ylabel)
    ax_bar.set_xlabel(xlabel)
    ax_bar.set_title(title)
    ax_bar.set_xlim(0, 100)
    for index, row in df.reset_index(drop=True).iterrows():
        ax_bar.text(
            row[percentage_column], index, f"{row[percentage_column]:.1f}%",
            ha="left", va="center", fontsize=PERCENTAGE_LABEL_FONTSIZE,
        )

    # Panels B, C, ...: part vs whole per category, in the frame's row order,
    # wrapped into a fixed two-per-row grid chained with below=.
    for panel_index, (label, (_, row)) in enumerate(zip(labels[1:], df.iterrows(), strict=True)):
        panel_kwargs: dict[str, object] = {
            "width": DONUT_PANEL_WIDTH_PX,
            "height": DONUT_PANEL_HEIGHT_PX,
            "margin_right": DONUT_PANEL_MARGIN_RIGHT_PX,
        }
        if panel_index >= DONUTS_PER_ROW:
            panel_kwargs["below"] = labels[1 + panel_index - DONUTS_PER_ROW]

        ax_donut = multipanel.panel(label=label, **panel_kwargs)

        category = str(row[category_column])
        part, whole = int(row[part_column]), int(row[whole_column])
        if part + whole == 0:
            ax_donut.text(0.5, 0.5, f"{category}\n(no data)", ha="center", va="center",
                          transform=ax_donut.transAxes, fontsize=PERCENTAGE_LABEL_FONTSIZE)
            continue

        draw_donut_panel(
            ax_donut, part, whole,
            part_label=part_label, whole_label=whole_label,
            center_text=f"{_wrap_center_text(category)}\n(n = {part + whole:,})",
        )

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
