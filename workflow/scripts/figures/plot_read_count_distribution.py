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

from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from figure_render.read_counts import (  # noqa: E402
    load_cutoff_stats,
    load_distribution_data,
    render_distribution_figure,
)


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
