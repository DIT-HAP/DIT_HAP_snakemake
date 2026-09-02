#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Insertion Orientation Figure (extraction + rendering, single stage)
===================================================================

One-stage QC figure: read the strand-indexed filtered insertion tables
produced by ``hard_filtering`` and render log-log correlation plots of ``+``
versus ``-`` strand insertion counts directly — no intermediate pairs TSV.

Data assembly (in-script): each input file is read with a 4-level row index
and a 2-level column header (Sample, Timepoint), both column levels are
stacked, ``Strand`` is pivoted back out, and rows lacking one strand or with a
non-positive count are dropped — the same pairing semantics the former
``qc.orientation.extract_strand_pairs`` applied before the intermediate table
was retired. Rendering is delegated to the shared grouped-regression renderer:
one panel per sample x timepoint with an identity guide line and n/r/P stats.

Input
-----
- One or more filtered insertion TSV files (``-i/--input``): strand-indexed
  tables with 4-level row index and 2-level column header (Sample, Timepoint).

Output
------
- ``<output>.pdf`` — journal-quality vector figure with regression grid.
- ``<output>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_insertion_orientation.py -i filtered/*.tsv -o figures/insertion_orientation
    python plot_insertion_orientation.py -i filtered/*.tsv -o figures/insertion_orientation --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-28
Version:  2.0.0
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
from io_tables import read_table  # noqa: E402
from figure_render.scatter import render_grouped_regression_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


# =============================================================================
# CONSTANTS
# =============================================================================
READER_KWARGS = {"index_col": [0, 1, 2, 3], "header": [0, 1]}
REQUIRED_COLUMNS = ["sample", "timepoint", "plus_count", "minus_count"]

X_COLUMN = "plus_count"
Y_COLUMN = "minus_count"
X_LABEL = "(+) strand"
Y_LABEL = "(−) strand"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding the validated input TSV path and output stem."""
    input_path: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that the input exists and output directory can be created."""
        if not self.input_path.exists():
            raise ValueError(f"Input file does not exist: {self.input_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def assemble_strand_pairs(input_path: Path) -> pd.DataFrame:
    """Read the strand-indexed filtered insertion table and assemble the long-format pairs frame.

    The pipeline always feeds the single aggregated ``hard_filtering`` output
    (all samples x timepoints ride its (Sample, Timepoint) column MultiIndex),
    so no filename parsing is needed. Raises on malformed or empty input.
    """
    logger.info(f"--- Processing file: {input_path.name} ---")
    df = read_table(input_path, **READER_KWARGS)
    pairs = _extract_strand_pairs(df)
    if pairs.empty:
        raise ValueError(f"No valid strand pairs in {input_path.name}!")
    logger.info(f"Extracted {len(pairs)} pairs from {input_path.name}")
    return pairs


def _extract_strand_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Pair +/- strand counts per Sample x Timepoint, keeping rows with both strands positive.

    Formerly ``qc.orientation.extract_strand_pairs``; inlined here because its
    only remaining purpose is feeding this figure's in-memory assembly.
    """
    required_levels = {"Strand", "Sample", "Timepoint"}
    available_levels = set(df.index.names) | set(df.columns.names)
    missing_levels = required_levels - available_levels
    if missing_levels:
        raise ValueError(f"Missing required index/column levels: {sorted(missing_levels)}")

    # Stack both column levels then pivot Strand back out, so each row holds the
    # +/- pair. dropna(axis=0) removes keys lacking one strand, exactly as before.
    plus_minus_pair = (
        df.stack(future_stack=True).stack(future_stack=True).unstack("Strand").dropna(axis=0)
    )
    logger.info(f"Paired {len(plus_minus_pair)} strand rows before positivity filtering")

    missing_strands = [strand for strand in ("+", "-") if strand not in plus_minus_pair.columns]
    if missing_strands:
        raise ValueError(f"Missing strand columns after unstacking: {missing_strands}")

    # Both strands strictly positive, which is what makes the log-log
    # correlation well defined.
    filtered = plus_minus_pair[plus_minus_pair.min(axis=1) > 0]
    logger.info(f"Retained {len(filtered)} rows with both strands strictly positive")

    if filtered.empty:
        logger.warning("No valid strand pairs after positivity filtering!")
        return pd.DataFrame(columns=REQUIRED_COLUMNS)

    pairs = (
        filtered[["+", "-"]]
        .rename(columns={"+": "plus_count", "-": "minus_count"})
        .reset_index()
        .rename(columns={"Sample": "sample", "Timepoint": "timepoint"})
    )
    return pairs[REQUIRED_COLUMNS]


def _prepare(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
    """Validate schema and return the frame; renderer will filter positive values for scale="log"."""
    require_columns(df, REQUIRED_COLUMNS, context=context)
    logger.info(f"Loaded {len(df)} rows ({context})")
    return df


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render insertion orientation figure from strand-indexed insertion tables")
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Filtered insertion TSV with Sample/Timepoint column headers")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Assemble strand pairs in memory, render the orientation figure, save dual artifacts."""
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
        pairs_df = assemble_strand_pairs(config.input_path)
        df = _prepare(pairs_df, context="in-memory strand pairs")

        render_grouped_regression_figure(
            df,
            config.output_stem,
            x=X_COLUMN,
            y=Y_COLUMN,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            row_key="sample",
            col_key="timepoint",
            density=True,
            scale="log",
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
