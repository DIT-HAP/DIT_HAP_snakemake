#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
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
PBL-PBR Correlation Analysis
============================

Analyse the correlation between PBL and PBR counts across one or more TSV
files produced by the DIT-HAP pipeline. Each input file is read with a
three-level row index, filtered to strictly positive PBL/PBR pairs (so the
values are valid on a log axis), and rendered as a log-log scatter plot with
a regression/diagonal reference line via ``create_scatter_correlation_plot``.

Files are processed in alphabetical order by filename, and one plot per file
is written as a separate page into a single multi-page PDF. A short run
summary (total data points, files processed) is logged at the end.

Input
-----
- One or more TSV files (``-i``/``--input``), tab-separated, with a
  three-column MultiIndex (``index_col=[0, 1, 2]``) and ``PBL`` and ``PBR``
  data columns.

Output
------
- A single PDF file (``-o``/``--output``) with one correlation plot per
  valid input file, rasterised scatter points on log-log axes.

Usage
-----
    python PBL_PBR_correlation_analysis.py -i file1.tsv file2.tsv -o results.pdf
    python PBL_PBR_correlation_analysis.py -i file1.tsv -o results.pdf --verbose

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
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# 4. Local module import (requires runtime injection of workflow/src on sys.path)
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from plot import create_scatter_correlation_plot  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
STYLE_PATH = str((SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(STYLE_PATH)
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PBLPBRCorrelationConfig:
    """Immutable config holding validated input TSV paths and the output PDF path."""
    input_files: list[Path]
    output_path: Path

    def __post_init__(self) -> None:
        """Validate that every input file exists and ensure the output directory is present."""
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru to emit uncolorised, timestamped records to stdout."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def read_tsv_file(file_path: Path) -> pd.DataFrame | None:
    """Read a TSV file and return only strictly-positive PBL/PBR pairs, or None if invalid."""
    logger.info(f"Reading TSV file: {file_path}")

    df = pd.read_csv(file_path, sep='\t', index_col=[0, 1, 2])

    # Check if PBL and PBR columns exist
    required_cols = ['PBL', 'PBR']
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        logger.warning(f"Warning: Missing columns {missing_cols} in {file_path}")
        return None

    # Remove rows with missing values in PBL or PBR
    df_clean = df[['PBL', 'PBR']].dropna()

    # Remove zero or negative values for log scaling
    df_clean = df_clean[(df_clean['PBL'] > 0) & (df_clean['PBR'] > 0)]

    if df_clean.empty:
        logger.warning(f"Warning: No valid data points in {file_path}")
        return None

    return df_clean

@logger.catch
def create_correlation_plot(filename: str, df: pd.DataFrame) -> plt.Figure:
    """Create a single log-log PBL-vs-PBR correlation figure for one file."""
    fig, ax = plt.subplots(figsize=(AX_WIDTH, AX_HEIGHT))

    ax = create_scatter_correlation_plot(
        x=df['PBL'],
        y=df['PBR'],
        ax=ax,
        xscale='log',
        yscale='log'
    )

    # Customize the plot
    ax.set_xlabel('PBL (log scale)')
    ax.set_ylabel('PBR (log scale)')
    ax.set_title(f'PBL vs PBR Correlation Analysis\n{filename}')

    return fig

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Analyze correlation between PBL and PBR from multiple TSV files")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="Input TSV files (space-separated)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

def main() -> int:
    """Read PBL/PBR TSVs, build per-file correlation plots, and write them to a multi-page PDF."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Validate input and output paths via the config dataclass (raises ValueError on bad input)
    try:
        config = PBLPBRCorrelationConfig(
            input_files=args.input,
            output_path=args.output,
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    logger.info("=== PBL-PBR Correlation Analysis ===")
    logger.info(f"Processing {len(config.input_files)} input files...")

    # Sort files by name
    sorted_files = sorted(config.input_files, key=lambda x: x.name)
    logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")

    # Read and process files
    data_dict: dict[str, pd.DataFrame] = {}
    for file_path in sorted_files:
        filename = file_path.name
        logger.info(f"Reading {filename}...")

        df = read_tsv_file(file_path)
        if df is not None:
            data_dict[filename] = df

    if not data_dict:
        logger.error("Error: No valid data found in any input file!")
        return 1

    # Create and save plots
    logger.info("Creating correlation plots...")

    # Save to PDF with rasterization
    logger.info(f"Saving plots to {config.output_path}...")
    try:
        with PdfPages(config.output_path) as pdf:
            for filename, df in data_dict.items():
                logger.info(f"  - Processing {filename}...")
                fig = create_correlation_plot(filename, df)
                pdf.savefig(fig)
                plt.close(fig)  # Close figure to free memory

        logger.success(f"Analysis complete! Output saved to: {config.output_path}")
        logger.info(f"Generated {len(data_dict)} correlation plots in PDF")

        # Print summary statistics
        logger.info("\n=== Summary Statistics ===")
        total_points = sum(len(df) for df in data_dict.values())
        logger.info(f"Total data points analyzed: {total_points}")
        logger.info(f"Files processed: {len(data_dict)}")

    except Exception as e:
        logger.error(f"Error saving plots: {str(e)}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
