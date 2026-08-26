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

import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figure_render.composition import render_composition_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


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
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["category", "covered", "not_covered", "total", "coverage_pct"]

COVERED_LABEL = "Covered"
NOT_COVERED_LABEL = "Not covered"

X_LABEL = "Gene viability"
Y_LABEL = "Coverage (%)"
TITLE = "Gene coverage by viability"
DONUT_UNIT = "genes"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_coverage_data(input_path: Path) -> pd.DataFrame:
    """Load the coverage statistics TSV and validate its schema."""
    logger.info(f"Loading coverage statistics from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"coverage TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} viability categories")
    return df


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

        render_composition_figure(
            df, config.output_stem,
            category_column="category", percentage_column="coverage_pct",
            part_column="covered", whole_column="not_covered", total_column="total",
            part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
            xlabel=X_LABEL, ylabel=Y_LABEL, title=TITLE, donut_unit=DONUT_UNIT,
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
