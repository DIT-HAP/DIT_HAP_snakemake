#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Replicate-Branch MA Plot Figure Renderer
========================================

Render DESeq2 MA plots from a pre-computed TSV with one square panel per
timepoint in a single row: log2 fold change against the mean of normalized
counts, points coloured darkred where the adjusted p-value clears the
significance threshold and gray otherwise. This is the rendering half of
pydeseq2's ``DeseqStats.plot_MA()``.

This is the biological-replicates branch, distinct from ``plot_ma_plot.py``:
that script handles the no-replicate branch, which computes M/A directly from
median-normalized counts, whereas this one carries DESeq2's ``baseMean`` /
``log2FoldChange`` / ``padj``.

Rows with NaN ``padj`` (pydeseq2's independent filtering removes them from
testing) are kept and rendered gray, matching pydeseq2, where ``nan < alpha``
evaluates False.

Input
-----
- MA-values figure-data TSV written by ``write_ma_values_tsv()`` in
  ``insertion_level_depletion_analysis_has_replicates.py``: tab-separated, no
  index, columns ``timepoint``, ``baseMean``, ``log2FoldChange``, ``padj``;
  one row per insertion per non-initial timepoint.

Output
------
- ``<stem>.pdf`` — journal-quality vector figure (scatter rasterized).
- ``<stem>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_ma_plot_replicates.py -i ma_values_replicates.tsv -o figures/ma_plot
    python plot_ma_plot_replicates.py -i ma_values_replicates.tsv -o figures/ma_plot --verbose

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
from figure_render.ma import Orientation, render_ma_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV path and output stem."""
    ma_values_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that input exists and output directory can be created."""
        if not self.ma_values_path.exists():
            raise ValueError(f"MA values file does not exist: {self.ma_values_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["timepoint", "baseMean", "log2FoldChange", "padj"]

# pydeseq2's DeseqStats default alpha; the project config sets no override.
PADJ_THRESHOLD = 0.05
SIGNIFICANT_COLOR = "darkred"
NONSIGNIFICANT_COLOR = "gray"

ABUNDANCE_LABEL = "mean of normalized counts"
EFFECT_LABEL = "log2 fold change"
TITLE_PREFIX = "MA plot"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_ma_data(ma_values_path: Path) -> pd.DataFrame:
    """Load the long-format MA values TSV and validate its schema."""
    logger.info(f"Loading MA values from {ma_values_path}...")
    df = pd.read_csv(ma_values_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"MA values TSV {ma_values_path.name}")

    logger.info(f"Loaded {len(df)} rows")
    return df


def significance_colors(padj: pd.Series) -> pd.Series:
    """Map adjusted p-values to point colours; NaN padj stays non-significant like pydeseq2."""
    # (padj < threshold) is False for NaN, so filtered-out insertions render gray.
    return (padj < PADJ_THRESHOLD).map({True: SIGNIFICANT_COLOR, False: NONSIGNIFICANT_COLOR})


def build_ma_panels(
    df: pd.DataFrame,
) -> tuple[list[tuple[str, pd.Series, pd.Series]], list[pd.Series]]:
    """Split the long table into per-timepoint triples plus their point colours."""
    panels: list[tuple[str, pd.Series, pd.Series]] = []
    colors: list[pd.Series] = []

    for timepoint in sorted(df["timepoint"].unique()):
        group = df[df["timepoint"] == timepoint]
        panels.append((timepoint, group["baseMean"], group["log2FoldChange"]))
        colors.append(significance_colors(group["padj"]))

    return panels, colors


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render replicate-branch MA plot figure from pre-computed MA values TSV")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input MA values TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load MA data, render figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            ma_values_path=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Replicate-Branch MA Plot Figure Rendering ===")

    try:
        # Load data
        df = load_ma_data(config.ma_values_path)
        panels, colors = build_ma_panels(df)

        render_ma_figure(
            panels, config.output_stem,
            abundance_label=ABUNDANCE_LABEL, effect_label=EFFECT_LABEL,
            title_prefix=TITLE_PREFIX,
            orientation=Orientation.HORIZONTAL,
            point_colors=colors,
        )

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
