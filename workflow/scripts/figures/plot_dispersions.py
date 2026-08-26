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

import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from figure_render.series import Series, render_series_scatter_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


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
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

X_COLUMN = "normed_mean"

# Series order, labels and colours mirror pydeseq2's plot_dispersions(): it passes
# [genewise, MAP, fitted] with labels ["Estimated", "Final", "Fitted"] and the
# matplotlib colour string "kbr" mapped positionally.
DISPERSION_SERIES = [
    Series(column="genewise_dispersion", label="Estimated", color="k"),
    Series(column="MAP_dispersion", label="Final", color="b"),
    Series(column="fitted_dispersion", label="Fitted", color="r"),
]

REQUIRED_COLUMNS = [X_COLUMN, *(item.column for item in DISPERSION_SERIES)]

X_LABEL = "mean of normalized counts"
Y_LABEL = "dispersion"
TITLE = "DESeq2 dispersion estimates"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_dispersion_data(dispersion_data_path: Path) -> pd.DataFrame:
    """Load the dispersion figure-data TSV and validate its schema."""
    logger.info(f"Loading dispersion data from {dispersion_data_path}...")
    df = pd.read_csv(dispersion_data_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)

    require_columns(df, REQUIRED_COLUMNS, context=f"dispersion TSV {dispersion_data_path.name}")

    logger.info(f"Loaded {len(df)} rows")
    return df


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
        render_series_scatter_figure(
            df, config.output_stem,
            x=X_COLUMN, series=DISPERSION_SERIES,
            xlabel=X_LABEL, ylabel=Y_LABEL, title=TITLE,
        )

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
