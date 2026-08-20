#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene Coverage Figure Renderer
==============================

Render a grouped bar chart of coverage percentage per viability category,
followed by one donut chart per category (covered vs not-covered), from the
pre-computed coverage statistics TSV.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
REQUIRED_COLUMNS = ["category", "covered", "not_covered", "total", "coverage_pct"]

# Panel labels: A for the bar overview, B/C/D/... for per-category donuts.
DONUT_STATUS_ORDER = ["Covered", "Not covered"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV path and output stem."""
    input_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.input_path.exists():
            raise ValueError(f"Input file does not exist: {self.input_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)




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


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render gene coverage figure from coverage statistics TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input coverage statistics TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load coverage statistics, render the figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            input_path=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Gene Coverage Figure Rendering ===")

    try:
        df = load_coverage_data(config.input_path)
        render_coverage_figure(df, config.output_stem)
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
