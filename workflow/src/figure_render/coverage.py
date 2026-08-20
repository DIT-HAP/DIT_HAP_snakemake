"""Gene coverage figure rendering.

Input
-----
- Coverage statistics TSV with columns ``category``, ``covered``,
  ``not_covered``, ``total``, ``coverage_pct``.

Output
------
- ``<stem>.pdf`` / ``<stem>.review.png`` — grouped bar chart of coverage
  percentage per viability category, followed by one donut chart per
  category (covered vs not-covered).

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
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

# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["category", "covered", "not_covered", "total", "coverage_pct"]

# Panel labels: A for the bar overview, B/C/D/... for per-category donuts.
DONUT_STATUS_ORDER = ["Covered", "Not covered"]


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_coverage_data(input_path: Path) -> pd.DataFrame:
    """Load the coverage statistics TSV and validate its schema."""
    logger.info(f"Loading coverage statistics from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} viability categories")
    return df


def _status_counts_frame(covered: int, not_covered: int) -> pd.DataFrame:
    """Expand covered/not-covered counts into a long-format frame for cns.donutplot."""
    return pd.DataFrame(
        {"status": ["Covered"] * covered + ["Not covered"] * not_covered}
    )


@logger.catch
def render_coverage_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render the coverage bar chart plus one donut chart per viability category."""
    logger.info("Rendering gene coverage figure...")

    apply_house_style()

    if df.empty:
        logger.warning("No coverage data to plot!")
        cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
        multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
        ax = multipanel.panel(label="A")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        save_dual(output_stem)
        return

    n_categories = len(df)
    n_panels = 1 + n_categories

    logger.info(f"Creating figure with {n_panels} panels ({n_categories} categories)...")

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    panel_labels = [chr(65 + i) for i in range(n_panels)]

    # Panel A: bar chart of coverage percentage per category
    label = panel_labels.pop(0)
    ax_bar = multipanel.panel(label=label)
    cns.barplot(data=df, x="category", y="coverage_pct", ax=ax_bar)
    ax_bar.set_ylabel("Coverage (%)")
    ax_bar.set_xlabel("Gene viability")
    ax_bar.set_title("Gene coverage by viability")
    ax_bar.set_ylim(0, 100)
    for idx, row in df.reset_index(drop=True).iterrows():
        ax_bar.text(idx, row["coverage_pct"], f"{row['coverage_pct']:.1f}%", ha="center", va="bottom", fontsize=6)

    # Panels B, C, ...: donut chart of covered vs not-covered genes per category
    for _, row in df.iterrows():
        label = panel_labels.pop(0)
        ax_donut = multipanel.panel(label=label)

        covered, not_covered = int(row["covered"]), int(row["not_covered"])
        if covered + not_covered == 0:
            ax_donut.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax_donut.transAxes)
            ax_donut.set_title(str(row["category"]))
            continue

        status_df = _status_counts_frame(covered, not_covered)
        cns.donutplot(data=status_df, x="status", order=DONUT_STATUS_ORDER, ax=ax_donut)
        ax_donut.set_title(f"{row['category']}\n({covered:,}/{int(row['total']):,} genes)")

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
