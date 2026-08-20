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

from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figure_render.density import load_density_data, render_density_figure  # noqa: E402


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
