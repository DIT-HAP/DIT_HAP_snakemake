"""Single-axes donut primitive built on ``cns.donutplot``.

``cns.donutplot`` hard-codes the column name as both the center annotation and
the legend title, which is meaningless for part/whole donuts driven by two raw
counts. This module owns the post-patching: center text replaced by caller
supplied text and the legend title dropped.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-28
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import matplotlib.axes
import matplotlib.text
import pandas as pd
from cnsplots import donutplot

# =============================================================================
# CONSTANTS
# =============================================================================
CENTER_FONTSIZE = 6

STATUS_COLUMN = "status"


# =============================================================================
# CORE LOGIC
# =============================================================================
def draw_donut_panel(
    ax: matplotlib.axes.Axes,
    part: int,
    whole: int,
    *,
    part_label: str,
    whole_label: str,
    center_text: str,
) -> matplotlib.text.Text:
    """Draw one part/whole donut on ``ax`` with ``center_text`` in the hole.

    Returns the center annotation so callers can style or assert on it. The
    ring starts at 12 o'clock and runs clockwise with ``part_label``'s wedge
    first, legend beneath. Counts of zero on both sides are handled by the
    caller — this primitive always draws a wedge ring, however degenerate.
    """
    counts = [part, whole]
    status_df = pd.DataFrame(
        {STATUS_COLUMN: [part_label] * part + [whole_label] * whole}
    )
    donutplot(
        data=status_df,
        x=STATUS_COLUMN,
        legend="bottom",
        order=[part_label, whole_label],
        ax=ax,
    )

    # cns.donutplot puts the column name ("status") at the center and as the
    # legend title; replace the former with the caller's text and drop the
    # latter, which would otherwise repeat the meaningless column name. The
    # panel-label Text may precede the center annotation in ax.texts, so pick
    # the annotation by its original column-name content instead of position.
    center = next(text for text in ax.texts if text.get_text() == STATUS_COLUMN)
    center.set_text(center_text)
    center.set_fontsize(CENTER_FONTSIZE)
    ax.get_legend().set_title(None)
    return center
