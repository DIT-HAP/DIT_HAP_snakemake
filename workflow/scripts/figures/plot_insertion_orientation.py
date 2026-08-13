#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Insertion Orientation Figure Renderer
=====================================

Render log-log correlation plots of ``+`` versus ``-`` strand insertion counts
from the pre-computed strand pairs TSV. One panel per sample and timepoint, with
regression line and correlation statistics via cnsplots.regplot.

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
import numpy as np
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
# Pixels each panel spends on its label, tick labels and axis title, on top of
# the axes box itself. Subtracted from the available width/height so a whole row
# of panels fits inside the journal width. Measured at 40 px for this figure:
# the log10 tick labels are single-digit, so panels carry narrower decorations
# than a frequency histogram. Overestimating shrinks the axes enough for an extra
# panel to fit in the row, which reflows the grid.
PANEL_DECORATION_PX = 40

# Verified to render ~50k points per panel legibly: density stays visible instead
# of collapsing into a solid block. Marker size must live inside scatter_kws --
# regplot's own `s` argument is silently dropped when scatter_kws is supplied.
SCATTER_KWS = {
    "s": 3,
    "facecolor": "none",
    "edgecolor": "gray",
    "alpha": 0.15,
    "linewidths": 0.25,
    "rasterized": True,
}


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
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru to emit uncolorised, timestamped records to stdout."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load strand pairs TSV, validate schema, filter positive values, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    required_cols = ["sample", "timepoint", "plus_count", "minus_count"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")

    # The computation layer already applies min > 0, but re-filtering keeps the
    # renderer safe against a hand-made input and guarantees finite log10 values.
    df_filtered = df[(df["plus_count"] > 0) & (df["minus_count"] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(df_filtered)} rows")

    if df_filtered.empty:
        logger.warning("No valid data points after filtering!")
        return df_filtered

    # regplot annotates r and P computed on the columns it is given, so the log10
    # columns must be explicit: passing raw counts would report a different
    # statistic than the log-log correlation the original figure showed.
    df_filtered["log10_plus_count"] = np.log10(df_filtered["plus_count"])
    df_filtered["log10_minus_count"] = np.log10(df_filtered["minus_count"])

    return df_filtered


@logger.catch
def render_orientation_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render a multipanel log-log strand correlation figure, one panel per sample-timepoint."""
    logger.info("Rendering insertion orientation figure...")

    apply_house_style()

    if df.empty:
        logger.warning("No data to plot!")
        return

    grouped = df.groupby(["sample", "timepoint"], sort=True)
    n_panels = len(grouped)

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    # Lay out one row per sample and one column per timepoint so samples stay on
    # their own rows rather than wrapping mid-sample.
    n_timepoints = df["timepoint"].nunique()
    n_samples = df["sample"].nunique()
    panel_width = max(45, int(JOURNAL_WIDTH_PX / max(n_timepoints, 1)) - PANEL_DECORATION_PX)
    panel_height = max(55, int(JOURNAL_HEIGHT_PX / max(n_samples, 1)) - PANEL_DECORATION_PX)

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    panel_index = 0

    for (sample, timepoint), group_df in grouped:
        label = chr(65 + panel_index) if panel_index < 26 else f"A{panel_index - 25}"
        column = panel_index % max(n_timepoints, 1)
        panel_index += 1

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            pad_left=2,
            pad_top=2,
            margin_right=6,
            margin_bottom=18,
        )

        # Titles carry the timepoint only, with the sample name as the row's
        # ylabel. A full "{sample} {timepoint}" title is wider than the axes and
        # pushes the measured panel width past the point where n_timepoints
        # panels fit in a row, which silently reflows the grid.
        ax.set_xlabel("log$_{10}$ (+) strand")
        ax.set_ylabel(f"{sample}\nlog$_{{10}}$ (-) strand" if column == 0 else "log$_{10}$ (-) strand")
        ax.set_title(timepoint)

        if group_df.empty:
            logger.warning(f"  Panel {label}: {sample} {timepoint} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            continue

        logger.info(f"  Panel {label}: {sample} {timepoint} (n={len(group_df)})")

        cns.regplot(
            data=group_df,
            x="log10_plus_count",
            y="log10_minus_count",
            scatter_kws=SCATTER_KWS,
            ax=ax,
        )

        # regplot resets the axis labels it manages, so reapply after drawing.
        ax.set_xlabel("log$_{10}$ (+) strand")
        ax.set_ylabel(f"{sample}\nlog$_{{10}}$ (-) strand" if column == 0 else "log$_{10}$ (-) strand")
        ax.set_title(timepoint)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render insertion orientation figure from strand pairs TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input strand pairs TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load strand pairs, render the orientation figure, and save dual artifacts."""
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

    logger.info("=== Insertion Orientation Figure Rendering ===")

    try:
        df = load_and_prepare_data(config.input_path)
        render_orientation_figure(df, config.output_stem)
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
