#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DESeq2 Dispersion Estimates Figure Renderer
===========================================

Render the DESeq2 dispersion diagnostic from a pre-computed TSV: genewise
(estimated), MAP (final) and fitted dispersions scattered against the mean of
normalized counts on log-log axes. This is the rendering half of pydeseq2's
built-in ``DeseqDataSet.plot_dispersions()``, reproducing its visual semantics
(series order, legend labels, black/blue/red colours, both axes log-scaled).

Both axes are log-scaled because dispersion estimates bottom out at pydeseq2's
``min_disp`` floor of 1e-8 and span many orders of magnitude; on linear axes the
whole cloud collapses onto the x-axis.

Input
-----
- Dispersion figure-data TSV written by ``write_dispersion_data_tsv()`` in
  ``insertion_level_depletion_analysis_has_replicates.py``: tab-separated, with
  index columns ``Chr, Coordinate, Strand, Target`` plus value columns
  ``normed_mean``, ``genewise_dispersion``, ``MAP_dispersion``,
  ``fitted_dispersion``.

Output
------
- ``<stem>.pdf`` — journal-quality vector figure (scatter rasterized).
- ``<stem>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_dispersions.py -i dispersion_data.tsv -o figures/dispersions
    python plot_dispersions.py -i dispersion_data.tsv -o figures/dispersions --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-16
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
X_COLUMN = "normed_mean"

# Series order, legend labels and colours mirror pydeseq2's plot_dispersions():
# it passes [genewise, MAP, fitted] with labels ["Estimated", "Final", "Fitted"]
# and the matplotlib colour string "kbr" mapped positionally.
DISPERSION_SERIES: tuple[tuple[str, str, str], ...] = (
    ("genewise_dispersion", "Estimated", "k"),
    ("MAP_dispersion", "Final", "b"),
    ("fitted_dispersion", "Fitted", "r"),
)

REQUIRED_COLUMNS = (X_COLUMN, *(column for column, _, _ in DISPERSION_SERIES))

X_LABEL = "mean of normalized counts"
Y_LABEL = "dispersion"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV path and output stem."""
    dispersion_data_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.dispersion_data_path.exists():
            raise ValueError(f"Dispersion data file does not exist: {self.dispersion_data_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)




# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_dispersion_data(dispersion_data_path: Path) -> pd.DataFrame:
    """Load dispersion figure-data TSV and validate schema."""
    logger.info(f"Loading dispersion data from {dispersion_data_path}...")
    df = pd.read_csv(dispersion_data_path, sep='\t', index_col=[0, 1, 2, 3])

    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")
    return df


@logger.catch
def render_dispersion_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render the three dispersion series against normed mean on log-log axes."""
    logger.info("Rendering dispersion figure...")

    # Apply house style before creating figure
    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    ax = multipanel.panel()

    if df.empty:
        logger.warning("No data to plot!")
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel(X_LABEL)
        ax.set_ylabel(Y_LABEL)
        ax.set_title('DESeq2 dispersion estimates')
        save_dual(output_stem)
        return

    # Rasterize: ~93.6k points per series (3 series) would bloat the vector PDF.
    # Alpha is kept low because the three clouds overlap heavily — at pydeseq2's
    # own 0.5 the blue "Final" cloud paints over the black "Estimated" one.
    scatter_kws = dict(
        s=1.0,
        alpha=0.12,
        linewidths=0,
        rasterized=True,
    )

    x_values = df[X_COLUMN]
    for column, legend_label, color in DISPERSION_SERIES:
        logger.info(f"  Series {legend_label}: {column} (n={df[column].notna().sum()})")
        ax.scatter(x_values, df[column], c=color, label=legend_label, **scatter_kws)

    # Log both axes: dispersions floor at pydeseq2's min_disp of 1e-8
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_title('DESeq2 dispersion estimates')

    # Legend markers inherit the scatter's tiny size and low alpha, so scale them up
    legend = ax.legend(loc='best', frameon=False, markerscale=6)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render DESeq2 dispersion figure from pre-computed dispersion TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input dispersion data TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load dispersion data, render figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            dispersion_data_path=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== DESeq2 Dispersion Figure Rendering ===")

    try:
        # Load data
        df = load_dispersion_data(config.dispersion_data_path)

        # Render figure
        render_dispersion_figure(df, config.output_stem)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
