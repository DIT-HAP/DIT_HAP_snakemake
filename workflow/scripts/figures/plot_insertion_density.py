#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Insertion Density Figure Renderer
==================================

Render the initial-vs-final insertion density comparison from the
pre-computed density statistics TSV. Reproduces the four-panel scatter
comparison of the pre-refactor multi-page PDF (initial-vs-final density,
density depletion vs. initial coverage, initial-vs-final read depth,
initial-vs-final depth inequality), colored by PomBase gene viability
(FYPOviability) when available.

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
from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
VIABILITY_HUE_ORDER = ["viable", "inviable", "condition-dependent", "unknown"]

# Verified scatter kwargs for ~5k points with essentiality colour-coding.
# cns.scatterplot always passes edgecolor=None to seaborn internally, so
# 'edgecolor' must not be supplied here or seaborn raises on the duplicate.
SCATTER_KWS = dict(s=6, alpha=0.4, rasterized=True)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV path, output stem, and timepoint labels."""
    input_path: Path
    output_stem: Path
    initial_timepoint: str
    final_timepoint: str

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.input_path.exists():
            raise ValueError(f"Input file does not exist: {self.input_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)




# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_density_data(input_path: Path) -> pd.DataFrame:
    """Load the density statistics TSV and validate its schema."""
    logger.info(f"Loading density statistics from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    required_cols = [
        "insertion_density_per_kb_initial",
        "insertion_density_per_kb_final",
        "insertion_density_log2fc",
        "total_reads_initial",
        "total_reads_final",
        "gini_coefficient_of_depth_initial",
        "gini_coefficient_of_depth_final",
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} genes")
    return df


def _draw_panel_or_placeholder(ax, df: pd.DataFrame, x: str, y: str, xlabel: str, ylabel: str, title: str) -> None:
    """Draw an essentiality-coloured scatter panel, or a 'No valid data' placeholder if empty."""
    if df.empty:
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        return

    if "FYPOviability" in df.columns:
        hue_order = [level for level in VIABILITY_HUE_ORDER if level in df["FYPOviability"].unique()]
        cns.scatterplot(data=df, x=x, y=y, hue="FYPOviability", hue_order=hue_order, ax=ax, **SCATTER_KWS)
        ax.legend(fontsize=6, loc="best", frameon=False)
    else:
        cns.scatterplot(data=df, x=x, y=y, ax=ax, **SCATTER_KWS)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


@logger.catch
def render_density_figure(df: pd.DataFrame, output_stem: Path, initial_timepoint: str, final_timepoint: str) -> None:
    """Render the four-panel initial-vs-final insertion density comparison."""
    logger.info("Rendering insertion density figure...")

    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    # Panel A: initial vs. final density, with the y = x reference line
    ax_a = multipanel.panel(label="A")
    _draw_panel_or_placeholder(
        ax_a, df,
        "insertion_density_per_kb_initial", "insertion_density_per_kb_final",
        f"Insertion density per kb ({initial_timepoint})",
        f"Insertion density per kb ({final_timepoint})",
        "Initial vs. Final Insertion Density",
    )
    if not df.empty:
        max_val = max(
            df["insertion_density_per_kb_initial"].max(),
            df["insertion_density_per_kb_final"].max(),
            1,
        )
        ax_a.plot([0, max_val], [0, max_val], color="red", linestyle="--", alpha=0.6, linewidth=1)

    # Panel B: initial density vs. log2 fold-change, with the y = 0 reference line
    ax_b = multipanel.panel(label="B")
    _draw_panel_or_placeholder(
        ax_b, df,
        "insertion_density_per_kb_initial", "insertion_density_log2fc",
        f"Insertion density per kb ({initial_timepoint})",
        f"log2FC density ({final_timepoint} / {initial_timepoint})",
        "Density Depletion vs. Initial Coverage",
    )
    if not df.empty:
        ax_b.axhline(0, color="red", linestyle="--", alpha=0.6, linewidth=1)

    # Panel C: initial vs. final total reads (log scale)
    ax_c = multipanel.panel(label="C")
    _draw_panel_or_placeholder(
        ax_c, df,
        "total_reads_initial", "total_reads_final",
        f"Total reads ({initial_timepoint})",
        f"Total reads ({final_timepoint})",
        "Initial vs. Final Read Depth",
    )
    if not df.empty:
        ax_c.set_xscale("symlog")
        ax_c.set_yscale("symlog")

    # Panel D: initial vs. final Gini coefficient of read depth
    ax_d = multipanel.panel(label="D")
    _draw_panel_or_placeholder(
        ax_d, df,
        "gini_coefficient_of_depth_initial", "gini_coefficient_of_depth_final",
        f"Gini coefficient of depth ({initial_timepoint})",
        f"Gini coefficient of depth ({final_timepoint})",
        "Initial vs. Final Depth Inequality",
    )
    if not df.empty:
        ax_d.plot([0, 1], [0, 1], color="red", linestyle="--", alpha=0.6, linewidth=1)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render insertion density figure from density statistics TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input density statistics TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint column name")
    parser.add_argument("-f", "--final_timepoint", type=str, required=True, help="Final timepoint column name")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load density statistics, render the comparison figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            input_path=args.input,
            output_stem=output_stem,
            initial_timepoint=args.initial_timepoint,
            final_timepoint=args.final_timepoint,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Insertion Density Figure Rendering ===")

    try:
        df = load_density_data(config.input_path)
        render_density_figure(df, config.output_stem, config.initial_timepoint, config.final_timepoint)
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
