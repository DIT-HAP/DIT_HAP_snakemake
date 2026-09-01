#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PBL-PBR Correlation Figure (extraction + rendering, single stage)
================================================================

One-stage QC figure: read the per-sample/timepoint/condition merged insertion
tables produced by ``merge_strand_insertions`` and render log-log correlation
plots of PBL vs PBR directly — no intermediate pairs TSV.

Data assembly (in-script): each input file is read with a three-level row
index, filtered to strictly positive PBL/PBR pairs, and tagged with
sample/timepoint/condition parsed from its filename stem. Rendering is
delegated to the shared grouped-regression renderer; groups are one panel per
sample x timepoint with an identity guide line and n/r/P correlation stats.

Points are coloured by 2D kernel density: each panel carries ~10^5 insertions,
where a flat colour collapses the diagonal ridge into a solid block. Density is
estimated on the log10 columns, so it matches the axes the panel displays.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-27
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
from sample_metadata import parse_filename  # noqa: E402
from figure_render.scatter import render_grouped_regression_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


# =============================================================================
# CONSTANTS
# =============================================================================
# PBL/PBR-specific knowledge lives here, not in the shared renderer.
REQUIRED_COLUMNS = ["sample", "timepoint", "condition", "pbl", "pbr"]
PAIRS_COLUMNS = ["sample", "timepoint", "condition", "chr", "coordinate", "strand", "pbl", "pbr"]

VALUE_COLUMNS = ["pbl", "pbr"]
X_COLUMN = "pbl"
Y_COLUMN = "pbr"
X_LABEL = "PBL"
Y_LABEL = "PBR"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input TSV paths and output stem."""
    input_files: list[Path]
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that every input exists and output directory can be created."""
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load one long-format pairs TSV, keep strictly-positive pairs, add log10 columns.

    Kept for tests and pixel-baseline tooling that operate on archived pairs
    tables; the pipeline main() path builds the same frame in-memory via
    ``extract_pairs`` instead.
    """
    df = pd.read_csv(input_path, sep="\t")
    return _prepare(df, context=f"pairs TSV {input_path.name}")


@logger.catch(reraise=True)
def extract_pairs(input_files: list[Path]) -> pd.DataFrame:
    """Read merged insertion tables and assemble the long-format pairs frame.

    Each file contributes strictly-positive PBL/PBR rows with metadata from its
    filename; invalid files are skipped so one malformed output cannot kill the
    whole figure. Raises when no file yields usable data.
    """
    frames: list[pd.DataFrame] = []
    for file_path in sorted(input_files, key=lambda p: p.name):
        metadata = parse_filename(file_path)
        if metadata is None:
            continue

        sample_df = _read_positive_pairs(file_path)
        if sample_df is None:
            continue

        sample, timepoint, condition = metadata
        pairs_df = sample_df.rename(columns={
            "Chr": "chr", "Coordinate": "coordinate", "Strand": "strand",
            "PBL": "pbl", "PBR": "pbr",
        }).assign(sample=sample, timepoint=timepoint, condition=condition)

        pairs_df = pairs_df[PAIRS_COLUMNS]
        logger.info(f"Extracted {len(pairs_df)} pairs from {file_path.name}")
        frames.append(pairs_df)

    if not frames:
        raise ValueError("No valid data found in any input file!")

    combined = pd.concat(frames, ignore_index=True)
    logger.info(f"Total pairs extracted: {len(combined)}")
    return combined


def _read_positive_pairs(file_path: Path) -> pd.DataFrame | None:
    """Read one merged table, keeping Chr/Coordinate/Strand plus strictly-positive PBL/PBR pairs.

    Returns None (with a warning) for files missing required columns or holding
    no valid data points.
    """
    logger.info(f"Reading TSV file: {file_path}")

    df = pd.read_csv(file_path, sep="\t", index_col=[0, 1, 2])

    required_cols = ["PBL", "PBR"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.warning(f"Missing columns {missing_cols} in {file_path}")
        return None

    df_clean = df[required_cols].dropna()
    df_clean = df_clean[(df_clean["PBL"] > 0) & (df_clean["PBR"] > 0)]

    if df_clean.empty:
        logger.warning(f"No valid data points in {file_path}")
        return None

    # Restore Chr/Coordinate/Strand from the row index as plain columns
    return df_clean.reset_index()


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
    parser = argparse.ArgumentParser(description="Render PBL-PBR correlation figure from merged insertion tables")
    parser.add_argument("-i", "--input", type=Path, nargs="+", required=True,
                        help="Input merged insertion TSV files ({sample}_{timepoint}_{condition}.tsv)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Extract pairs from inputs in-memory, render correlation figure, save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            input_files=args.input,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== PBL-PBR Correlation Figure Rendering ===")

    try:
        pairs_df = extract_pairs(config.input_files)
        df = _prepare(pairs_df, context="in-memory pairs")

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
