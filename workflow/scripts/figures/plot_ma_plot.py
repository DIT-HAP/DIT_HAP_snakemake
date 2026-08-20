#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MA Plot Figure Renderer (Unified)
==================================

Render MA plots directly from the wide-format ``baseMean.tsv`` and ``LFC.tsv``
tables produced by either insertion-level depletion analysis branch (with or
without replicates), with one panel per timepoint column. Both tables share
the same row index (Chr, Coordinate, Strand, Target) and the same timepoint
columns, so no branch-specific handling is needed.

Input
-----
- baseMean TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex, tab-separated,
  one column per timepoint.
- LFC TSV: same row MultiIndex and timepoint columns as the baseMean table.

Output
------
- `<stem>.pdf` / `<stem>.review.png` — vertical stack, one panel per timepoint,
  LFC on the x-axis and baseMean (log-scaled) on the y-axis.
- `<stem>_horizontal.pdf` / `<stem>_horizontal.review.png` — the same panels
  arranged in a single row, with baseMean (log-scaled) on the x-axis and LFC
  on the y-axis.

Usage
-----
    python plot_ma_plot.py -b baseMean.tsv -l LFC.tsv -o figures/ma_plot
    python plot_ma_plot.py -b baseMean.tsv -l LFC.tsv -o figures/ma_plot --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  4.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figure_render.ma_plot import Orientation, load_ma_data, render_ma_figure  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV paths and output stem."""
    basemean_path: Path
    lfc_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that inputs exist and output directory can be created."""
        for path in (self.basemean_path, self.lfc_path):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)




# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render unified MA plot figure from baseMean/LFC tables")
    parser.add_argument("-b", "--basemean", type=Path, required=True, help="Input baseMean TSV file")
    parser.add_argument("-l", "--lfc", type=Path, required=True, help="Input LFC TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load MA data, render figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix("")

    # Validate paths
    try:
        config = PlotConfig(
            basemean_path=args.basemean,
            lfc_path=args.lfc,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Unified MA Plot Figure Rendering ===")

    try:
        # Load data
        basemean_df, lfc_df = load_ma_data(config.basemean_path, config.lfc_path)

        # Render both orientations
        render_ma_figure(basemean_df, lfc_df, config.output_stem, Orientation.VERTICAL)
        horizontal_stem = config.output_stem.with_name(f"{config.output_stem.name}_horizontal")
        render_ma_figure(basemean_df, lfc_df, horizontal_stem, Orientation.HORIZONTAL)

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
