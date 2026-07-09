#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "numpy",
#     "matplotlib",
#     "loguru",
# ]
# ///

"""
Insertion Orientation Analysis
==============================

Analyse insertion strand orientation (+/-) pairs from one or more TSV files
that use multi-level row and column indexing. For every input file the script
extracts the positive/negative strand pair for each (Sample, Timepoint) group
and renders a log-log scatter plot with correlation annotation.

Each Sample produces one figure containing one subplot per Timepoint (n rows x
1 column). Rows whose minimum strand value is not strictly positive are dropped
before plotting so both axes remain valid on a log scale. All figures are
written to a single multi-page PDF report.

Input
-----
- One or more TSV/TXT files with a 4-level row index (``index_col=[0, 1, 2, 3]``)
  and a 2-level column header (``header=[0, 1]``) in which the ``Strand`` level
  holds the ``+`` / ``-`` orientation.

Output
------
- A single multi-page PDF report: one page per Sample, one subplot per
  Timepoint, each a log-log scatter plot of positive vs negative strand counts.

Usage
-----
    python insertion_orientation_analysis.py -i file1.tsv file2.tsv -o orientation_analysis.pdf
    python insertion_orientation_analysis.py -i file1.tsv -o out.pdf --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# 4. Local application imports (require the project ``src`` dir on sys.path)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from plot import create_scatter_correlation_plot  # noqa: E402
from utils import read_file  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
STYLE_PATH = str((SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(STYLE_PATH)
AX_WIDTH, AX_HEIGHT = plt.rcParams["figure.figsize"]
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InsertionOrientationAnalysisConfig:
    """Validated input/output paths for the insertion orientation analysis."""

    input_files: list[Path]
    output_path: Path

    def __post_init__(self) -> None:
        """Validate that inputs exist with a TSV suffix and prepare the PDF output directory."""
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
            if file_path.suffix.lower() not in [".tsv", ".txt"]:
                raise ValueError(f"Input file must be a TSV file: {file_path}")
        if self.output_path.suffix.lower() != ".pdf":
            raise ValueError(f"Output file must be a PDF: {self.output_path}")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger to write to stdout at the given level."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


setup_logger()


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def create_file_comparison_figure(df: pd.DataFrame, filename: str, output_path: Path) -> None:
    """Create a multi-page figure comparing +/- strand values for every Sample/Timepoint."""
    plus_minus_pair = df.stack(future_stack=True).stack(future_stack=True).unstack("Strand").dropna(axis=0)
    timepoints = plus_minus_pair.index.get_level_values("Timepoint").unique()

    # Create figure with subplots: n rows x 1 column (one row per timepoint)
    n_rows = len(timepoints)
    fig_width = AX_WIDTH
    fig_height = n_rows * AX_HEIGHT

    total_figures = 0
    with PdfPages(output_path) as pdf:
        for sample, sample_df in plus_minus_pair.groupby(level="Sample"):
            fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height))
            if n_rows == 1:
                axes = [axes]  # Ensure axes is always a list

            # Process each timepoint
            for row_idx, (timepoint, sub_df) in enumerate(sample_df.groupby(level="Timepoint")):
                ax = axes[row_idx]

                filtered_sub_df = sub_df[sub_df.min(axis=1) > 0]

                pos_array = filtered_sub_df["+"].to_numpy()
                neg_array = filtered_sub_df["-"].to_numpy()

                create_scatter_correlation_plot(
                    x=pos_array,
                    y=neg_array,
                    ax=ax,
                    xscale="log",
                    yscale="log",
                )

                # Customize subplot
                ax.set_xlabel("Positive Strand (+)")
                if row_idx == 0:  # Only label y-axis on leftmost subplot
                    ax.set_ylabel("Negative Strand (-)")
                ax.set_title(f"Sample: {sample}\nTimepoint: {timepoint}")
                ax.grid(True)

            plt.tight_layout()
            # Save figure to PDF
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            total_figures += 1

    logger.info(f"Generated {total_figures} figures")


@logger.catch
def analyze_multiple_files(input_files: list[Path], output_path: Path) -> None:
    """Analyse strand orientations across multiple files and write a combined PDF report."""
    logger.info("Starting multi-file strand orientation analysis...")

    # Sort files by name for consistent processing order
    sorted_files = sorted(input_files, key=lambda p: p.name)
    logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")

    # Generate plots and save to PDF
    for file_path in sorted_files:
        filename = file_path.name
        logger.info(f"--- Processing file: {filename} ---")

        # Per-file control flow: skip a file that fails to process, continue with the rest.
        try:
            df = read_file(file_path, **{"index_col": [0, 1, 2, 3], "header": [0, 1]})
            create_file_comparison_figure(df, filename, output_path)
            logger.info(f"Generated figure for {filename}")
        except Exception as e:
            logger.error(f"Failed to process {filename}: {e}", exc_info=True)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Analyze insertion orientation (+/-) strand pairs from multiple TSV files.")
    parser.add_argument("-i", "--input", nargs="+", type=Path, required=True, help="One or more input TSV files with multi-level indexing.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path for the plots.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Orchestrate CLI parsing, validation and multi-file orientation analysis."""
    args = parse_args()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(log_level)

    try:
        config = InsertionOrientationAnalysisConfig(
            input_files=args.input,
            output_path=args.output,
        )

        logger.info("=== Insertion Orientation Analysis ===")
        logger.info(f"Processing {len(config.input_files)} input files...")

        start_time = time.time()
        analyze_multiple_files(config.input_files, config.output_path)
        end_time = time.time()
        total_time = end_time - start_time
        logger.info(f"Completed insertion orientation analysis in {total_time:.2f} seconds.")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
