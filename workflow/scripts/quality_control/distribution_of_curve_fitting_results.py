#!/usr/bin/env python3

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
Distribution of Curve-Fitting Results
======================================

Generate per-column histograms summarising the distribution of curve-fitting
result parameters (e.g. R-squared, inflection points, rate constants). The
script auto-detects every numeric column in the input table and plots one
histogram per column, so it adapts to whatever set of fitted parameters the
upstream fitting step produced.

Histograms are laid out on a fixed 4-column grid (rows grow with the number of
numeric columns) and rendered into a single multi-panel PDF via matplotlib's
``PdfPages`` using the project ``DIT_HAP.mplstyle`` theme. Each panel is
annotated with the sample count, mean, and standard deviation of the plotted
column. Columns that contain no non-null values render a "No data" placeholder
instead of a histogram.

Input
-----
- A tab-separated table (read with ``sep="\\t"``) whose numeric columns hold the
  fitted parameters to be summarised. Non-numeric columns are ignored.

Output
------
- A single PDF file containing the grid of histograms, one panel per numeric
  column, with an overall figure title "Distribution of Numeric Variables".

Usage
-----
    python distribution_of_curve_fitting_results.py --input results.tsv --output dist.pdf
    python distribution_of_curve_fitting_results.py -i results.tsv -o dist.pdf --bins 50 --verbose

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
from typing import Any

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# Matplotlib styling is configured once at import time so every figure inherits
# the project theme. AX_WIDTH / AX_HEIGHT / COLORS are read from the active
# rcParams for convenience; the plotting function re-reads them locally.
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths and histogram binning for a single run."""
    input_file: Path
    output_file: Path
    bins: int = 30

    def __post_init__(self) -> None:
        # Real validation: a missing input file is a fatal, caller-visible error.
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        # Side effect (no attribute assignment, frozen-safe): ensure the output
        # directory exists before matplotlib writes the PDF.
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Summary counts and timing for one distribution-plotting run."""
    total_rows: int
    total_columns: int
    numeric_columns: int
    plots_generated: int
    missing_values: int
    execution_time: float

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger to write to stdout at the given level."""
    logger.remove()  # Remove default handler
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
def load_and_analyze_data(input_file: Path) -> tuple[pd.DataFrame, list[str], dict[str, Any]]:
    """Load a TSV table and identify its numeric columns for analysis."""
    logger.info(f"Loading data from {input_file}")

    # Load data
    df = pd.read_csv(input_file, sep="\t")
    logger.info(f"Loaded {len(df)} rows and {len(df.columns)} columns")

    # Identify numeric columns
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    logger.info(f"Found {len(numeric_columns)} numeric columns: {numeric_columns}")

    # Generate basic statistics
    stats = {
        'total_rows': len(df),
        'total_columns': len(df.columns),
        'numeric_columns': len(numeric_columns),
        'missing_values': df.isnull().sum().sum()
    }

    return df, numeric_columns, stats


@logger.catch
def create_histogram_plots(df: pd.DataFrame, numeric_columns: list[str],
                           output_file: Path, bins: int = 30) -> dict[str, Any]:
    """Create histogram plots for numeric columns in a multi-panel PDF layout."""
    logger.info(f"Creating histogram plots with {bins} bins for {len(numeric_columns)} columns")

    # Calculate subplot layout (4 columns)
    n_cols = 4
    n_rows = (len(numeric_columns) + n_cols - 1) // n_cols  # Ceiling division

    logger.info(f"Creating {n_rows}x{n_cols} subplot layout")

    # Create figure
    plot_width, plot_height = plt.rcParams['figure.figsize']
    fig_width = plot_width * n_cols
    fig_height = plot_height * n_rows

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    with PdfPages(output_file) as pdf:
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))

        # Handle case where we have only one row or column
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)
        elif len(numeric_columns) == 1:
            axes = np.array([[axes]])

        plot_stats = {}

        for idx, column in enumerate(numeric_columns):
            row = idx // n_cols
            col = idx % n_cols
            ax = axes[row, col]

            # Get data and remove NaN values
            data = df[column].dropna()

            if len(data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(column)
                continue

            # Create histogram
            color = colors[idx % len(colors)]
            ax.hist(
                data,
                bins=bins,
                color=color,
                alpha=0.8,
                edgecolor='white',
                linewidth=0.5
            )

            # Customize plot
            ax.set_title(column)
            ax.set_xlabel('Value')
            ax.set_ylabel('Frequency')

            stats_text = f'n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                    verticalalignment='top')

            # Store statistics
            plot_stats[column] = {
                'count': len(data),
                'mean': data.mean(),
                'std': data.std(),
                'min': data.min(),
                'max': data.max(),
                'median': data.median()
            }

            logger.debug(f"Plotted {column}: {len(data)} values, range [{data.min():.3f}, {data.max():.3f}]")

        # Hide empty subplots
        for idx in range(len(numeric_columns), n_rows * n_cols):
            row = idx // n_cols
            col = idx % n_cols
            axes[row, col].set_visible(False)

        # Add overall title
        fig.suptitle('Distribution of Numeric Variables', y=1.02)

        # Save to PDF
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

    logger.info(f"Histogram plots saved to {output_file}")

    return plot_stats


@logger.catch
def analyze_curve_fitting_results(config: InputOutputConfig) -> AnalysisResult:
    """Load curve-fitting results, plot parameter distributions, and summarise the run."""
    start_time = time.time()

    # Load and analyze data
    df, numeric_columns, stats = load_and_analyze_data(config.input_file)

    if not numeric_columns:
        logger.warning("No numeric columns found in the input data")
        return AnalysisResult(
            total_rows=stats['total_rows'],
            total_columns=stats['total_columns'],
            numeric_columns=0,
            plots_generated=0,
            missing_values=stats['missing_values'],
            execution_time=time.time() - start_time
        )

    # Create histogram plots
    plot_stats = create_histogram_plots(df, numeric_columns, config.output_file, config.bins)

    # Create result object
    execution_time = time.time() - start_time
    result = AnalysisResult(
        total_rows=stats['total_rows'],
        total_columns=stats['total_columns'],
        numeric_columns=stats['numeric_columns'],
        plots_generated=len(plot_stats),
        missing_values=stats['missing_values'],
        execution_time=execution_time
    )

    logger.info(f"Analysis completed: {result.plots_generated} plots generated from {result.numeric_columns} numeric columns")
    logger.info(f"Data summary: {result.total_rows:,} rows, {result.total_columns} columns, {result.missing_values} missing values")
    logger.info(f"Execution time: {result.execution_time:.2f} seconds")

    return result

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Generate histogram plots for numeric columns in curve fitting results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input CSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output PDF file")
    parser.add_argument("--bins", type=int, default=30, help="Number of bins for histograms")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate paths, run the analysis, and report the outcome."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    # Validate input and output paths, then run the core analysis/logic.
    try:
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output,
            bins=args.bins
        )

        logger.info(f"Starting processing of {config.input_file}")

        results = analyze_curve_fitting_results(config)

        logger.success(f"Analysis complete. Results saved to {config.output_file}")
        logger.info(f"Generated {results.plots_generated} plots from {results.numeric_columns} numeric columns")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
