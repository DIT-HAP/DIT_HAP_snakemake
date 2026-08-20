#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read Count Distribution Figure Renderer
=======================================

Render pre-binned log10 read-count histograms with the hard-filtering cutoff
marked on the initial-timepoint panel. Consumes the binned distribution TSV and
the cutoff statistics TSV produced by the computation layer; it never re-bins.

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

import matplotlib.pyplot as plt  # noqa: E402

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual  # noqa: E402


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# Pixels each panel spends on its label, tick labels and axis title, on top of
# the axes box itself. Subtracted from the available width/height so a whole row
# of panels fits inside the journal width.
PANEL_DECORATION_PX = 55

# Height reserved per line of the retention footer at the bottom of the figure.
FOOTER_LINE_PX = 9


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV paths, output stem, and annotations."""
    input_path: Path
    stats_path: Path
    output_stem: Path
    initial_time_point: str
    cutoff: float

    def __post_init__(self) -> None:
        """Validate that inputs exist, the cutoff is positive, and the output directory is present."""
        if not self.input_path.exists():
            raise ValueError(f"Input file does not exist: {self.input_path}")
        if not self.stats_path.exists():
            raise ValueError(f"Stats file does not exist: {self.stats_path}")
        if self.cutoff <= 0:
            raise ValueError(f"Cutoff must be positive: {self.cutoff}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)




# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_distribution_data(input_path: Path) -> pd.DataFrame:
    """Load the binned distribution TSV, validate its schema, and add bin centres and widths."""
    logger.info(f"Loading binned distribution from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    required_cols = ["sample", "timepoint", "bin_left", "bin_right", "count"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} bin rows")

    if df.empty:
        return df

    # Bin centres and widths let histplot reproduce the pre-computed bins exactly
    # via weights, with no re-binning of underlying values.
    df["bin_center"] = (df["bin_left"] + df["bin_right"]) / 2.0
    df["bin_width"] = df["bin_right"] - df["bin_left"]

    return df


@logger.catch
def load_cutoff_stats(stats_path: Path) -> pd.DataFrame:
    """Load the cutoff retention statistics TSV and validate its schema."""
    logger.info(f"Loading cutoff statistics from {stats_path}...")
    df = pd.read_csv(stats_path, sep="\t")

    required_cols = [
        "sample",
        "original_rows",
        "rows_kept",
        "pct_rows_kept",
        "original_counts",
        "counts_kept",
        "pct_counts_kept",
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in stats: {missing_cols}")

    logger.info(f"Loaded statistics for {len(df)} samples")
    return df


def format_retention_caption(sample: str, stats_row: pd.Series | None) -> str:
    """Format a one-line cutoff retention summary for a sample, or an empty string when absent."""
    if stats_row is None:
        return ""

    return (
        f"{sample}: {int(stats_row['rows_kept']):,}/{int(stats_row['original_rows']):,} rows kept "
        f"({stats_row['pct_rows_kept']:.1f}%), "
        f"{int(stats_row['counts_kept']):,}/{int(stats_row['original_counts']):,} counts kept "
        f"({stats_row['pct_counts_kept']:.1f}%)"
    )


@logger.catch
def render_distribution_figure(
    df: pd.DataFrame,
    stats_df: pd.DataFrame,
    output_stem: Path,
    initial_time_point: str,
    cutoff: float,
) -> None:
    """Render one histogram panel per sample-timepoint group from pre-binned counts."""
    logger.info("Rendering read count distribution figure...")

    apply_house_style()

    if df.empty:
        logger.warning("No bins to plot!")
        return

    grouped = df.groupby(["sample", "timepoint"], sort=True)
    n_panels = len(grouped)

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    stats_by_sample = stats_df.set_index("sample")
    log_cutoff = float(np.log10(cutoff))

    # Lay out one row per sample and one column per timepoint. The axes width is
    # sized so exactly n_timepoints panels fit inside max_width once the label,
    # y-axis decorations and margins are accounted for, which keeps each sample
    # on its own row instead of wrapping mid-sample.
    n_timepoints = df["timepoint"].nunique()
    n_samples = df["sample"].nunique()
    panel_width = max(45, int(JOURNAL_WIDTH_PX / max(n_timepoints, 1)) - PANEL_DECORATION_PX)
    panel_height = max(55, int(JOURNAL_HEIGHT_PX / max(n_samples, 1)) - PANEL_DECORATION_PX)

    # multipanel resizes the figure to fit its panels, so a bottom strip must be
    # reserved inside the layout (via the last row's bottom margin) rather than by
    # shrinking the requested figure height, or the footer lands on top of the
    # last row's axis labels.
    footer_reserve_px = FOOTER_LINE_PX * (n_samples + 2)
    last_row_index = n_samples - 1

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    panel_index = 0
    retention_lines: list[str] = []

    for (sample, timepoint), group_df in grouped:
        label = chr(65 + panel_index) if panel_index < 26 else f"A{panel_index - 25}"
        is_last_row = (panel_index // max(n_timepoints, 1)) == last_row_index
        panel_index += 1

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            pad_left=2,
            pad_top=2,
            margin_right=4,
            margin_bottom=20 + (footer_reserve_px if is_last_row else 0),
        )

        valid = group_df.dropna(subset=["bin_left", "bin_right"])
        total_count = float(valid["count"].sum())

        if valid.empty or total_count <= 0:
            logger.warning(f"  Panel {label}: {sample} {timepoint} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlabel("log$_{10}$(read count)")
            ax.set_ylabel("Frequency")
            ax.set_title(f"{sample} {timepoint}")
            continue

        logger.info(f"  Panel {label}: {sample} {timepoint} ({len(valid)} bins, n={int(total_count):,})")

        # Replay the pre-computed bins: bin centres weighted by their counts, with
        # bin count and range taken from the stored edges so nothing is re-binned.
        # The computation layer always bins with a uniform integer bin count, so
        # (n_bins, binrange) reproduces the stored edges exactly. An explicit edge
        # array cannot be used here: seaborn compares `bins == "auto"` before
        # binning, which raises on an array when weights are supplied.
        binrange = (float(valid["bin_left"].min()), float(valid["bin_right"].max()))
        cns.histplot(
            data=valid,
            x="bin_center",
            weights="count",
            bins=len(valid),
            binrange=binrange,
            ax=ax,
        )

        ax.set_xlabel("log$_{10}$(read count)")
        ax.set_ylabel("Frequency")
        ax.set_title(f"{sample} {timepoint}")

        if timepoint != initial_time_point:
            continue

        # Cutoff annotation belongs only on the panel the cutoff is applied to.
        ax.axvline(log_cutoff, color="firebrick", linestyle="--", linewidth=0.8, label=f"Cutoff = {cutoff:.2g}")
        ax.legend(frameon=False, loc="upper right", fontsize=4)

        stats_row = stats_by_sample.loc[sample] if sample in stats_by_sample.index else None
        caption = format_retention_caption(sample, stats_row)
        if caption:
            retention_lines.append(caption)

    # Retention stats go in a figure-level footer rather than a per-panel textbox:
    # the numbers are per sample, not per panel, and a cramped in-axes box
    # overlapped the histograms in the original figure.
    if retention_lines:
        header = f"Cutoff applied to '{initial_time_point}' (>= {cutoff:.2g}):"
        footer = "\n".join([header, *retention_lines])
        plt.gcf().text(0.02, 0.004, footer, ha="left", va="bottom", fontsize=5, linespacing=1.4)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render read count distribution figure from binned TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input binned distribution TSV file")
    parser.add_argument("-s", "--stats", type=Path, required=True, help="Input cutoff statistics TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-t", "--initial_time_point", type=str, required=True, help="Initial time point column name for the cutoff annotation")
    parser.add_argument("-c", "--cutoff", type=float, required=True, help="Hard filtering cutoff value to annotate")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load binned distributions and stats, render the figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            input_path=args.input,
            stats_path=args.stats,
            output_stem=output_stem,
            initial_time_point=args.initial_time_point,
            cutoff=args.cutoff,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Read Count Distribution Figure Rendering ===")

    try:
        df = load_distribution_data(config.input_path)
        stats_df = load_cutoff_stats(config.stats_path)
        render_distribution_figure(
            df,
            stats_df,
            config.output_stem,
            config.initial_time_point,
            config.cutoff,
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
