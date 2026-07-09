#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "loguru",
#     "tabulate",
# ]
# ///

"""
Read Count Distribution Analysis
================================

Analyzes read count distributions across one or more TSV files. For each input
file the script builds log10-transformed histograms of every numeric column,
draws the hard-filtering cutoff on the initial time point column, and reports how
many rows and counts survive that cutoff. All plots are collected into a single
multi-page PDF (one page per input file), and a formatted summary table of the
per-file retention statistics is written to the log.

The cutoff keeps rows whose initial-time-point value is greater than or equal to
the cutoff; histograms are drawn on a log10 scale using only strictly positive
values, while the cutoff line is placed at log10(cutoff).

Input
-----
- One or more TSV files with a 4-level multi-index (``index_col=[0, 1, 2, 3]``)
  and one numeric read-count column per time point.

Output
------
- A multi-page PDF (one distribution plot per input file), each page carrying a
  per-file statistics box, plus a summary table logged to stdout.

Usage
-----
    python read_count_distribution_analysis.py -i sample1.tsv sample2.tsv -t T0 -c 10 -o report.pdf
    python read_count_distribution_analysis.py -i *.tsv -t T0 -c 10 -o report.pdf --bins 80 --verbose

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
import matplotlib.pyplot as plt
from loguru import logger
from matplotlib.backends.backend_pdf import PdfPages
from tabulate import tabulate

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams["figure.figsize"]
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# Headers for summary table
SUMMARY_HEADERS = {
    "filename": "File Name",
    "original_rows": "Original Rows",
    "original_counts": "Original Counts",
    "rows_kept": "Rows Kept",
    "percentage_rows_kept": "% Rows Kept",
    "count_kept": "Counts Kept",
    "percentage_count_kept": "% Counts Kept",
}

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ReadCountDistributionAnalysisConfig:
    """Immutable, validated configuration for read count distribution analysis."""

    input_files: list[Path]
    output_path: Path
    initial_time_point: str
    cutoff: float
    bins: int = 50

    def __post_init__(self) -> None:
        """Validate configuration values and create the output directory."""
        if not self.input_files:
            raise ValueError("At least one input file must be provided")
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
            if file_path.suffix.lower() not in [".tsv", ".txt"]:
                raise ValueError(f"Input file must be a TSV file: {file_path}")
        if self.output_path.suffix.lower() != ".pdf":
            raise ValueError(f"Output file must be a PDF: {self.output_path}")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        if self.cutoff <= 0:
            raise ValueError("Cutoff value must be positive")
        if self.bins < 5 or self.bins > 200:
            raise ValueError("Number of bins must be between 5 and 200")


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Immutable container for the results of analyzing a single file."""

    filename: str
    original_rows: int
    original_counts: float
    rows_kept: int | str
    percentage_rows_kept: float | str
    count_kept: float | str
    percentage_count_kept: float | str

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
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
def load_and_validate_data(file_path: Path) -> pd.DataFrame:
    """Load TSV file with multi-index structure and validate data integrity."""
    df_reader_kwargs = {"sep": "\t", "index_col": [0, 1, 2, 3]}
    df = pd.read_csv(file_path, engine="python", **df_reader_kwargs)

    if df.empty:
        raise ValueError("Empty DataFrame after reading")

    return df

@logger.catch
def calculate_cutoff_statistics(df: pd.DataFrame, initial_time_point_col: str, cutoff_val: float) -> dict[str, Any]:
    """Calculate filtering statistics after applying cutoff to initial time point column."""
    stats: dict[str, Any] = {}
    original_rows = len(df)
    original_counts = df[initial_time_point_col].sum()

    stats["original_rows"] = original_rows
    stats["original_counts"] = original_counts

    if initial_time_point_col not in df.columns:
        stats.update({
            "rows_kept": "N/A",
            "percentage_rows_kept": "N/A",
            "count_kept": "N/A",
            "percentage_count_kept": "N/A",
        })
        return stats

    if not pd.api.types.is_numeric_dtype(df[initial_time_point_col]):
        stats.update({
            "rows_kept": "N/A (col not numeric)",
            "percentage_rows_kept": "N/A",
            "count_kept": "N/A",
            "percentage_count_kept": "N/A",
        })
        return stats

    kept_df = df[df[initial_time_point_col] >= cutoff_val]
    stats["rows_kept"] = len(kept_df)
    stats["percentage_rows_kept"] = (stats["rows_kept"] / original_rows) * 100.0 if original_rows > 0 else 0.0
    stats["count_kept"] = kept_df[initial_time_point_col].sum()
    stats["percentage_count_kept"] = (stats["count_kept"] / original_counts) * 100.0 if original_counts > 0 else 0.0

    return stats

@logger.catch
def create_distribution_plot(df: pd.DataFrame, filename: str, initial_time_point_col: str,
                           cutoff_val: float, bins: int, stats: dict[str, Any]) -> plt.Figure:
    """Create log-transformed distribution plots with original scale labels and cutoff line."""
    numeric_cols = df.select_dtypes(include=np.number).columns
    if not numeric_cols.any():
        raise ValueError("No numeric columns found")

    num_subplots = len(numeric_cols)
    fig, axes = plt.subplots(
        num_subplots, 1,
        figsize=(AX_WIDTH, num_subplots * AX_HEIGHT),
        sharex=True,
        sharey=True,
    )
    if num_subplots == 1:
        axes = [axes]

    fig.suptitle(
        f"Value Distributions: {filename}\nInitial Time Point: '{initial_time_point_col}', Cutoff >= {cutoff_val}", y=1.01
    )

    max_y_val = 0

    for i, col_name in enumerate(numeric_cols):
        ax = axes[i]
        col_data = df[col_name].dropna()
        positive_col_data = col_data[col_data > 0]
        ax.set_xlabel("Log10(Value)")
        ax.set_ylabel("Frequency")

        if positive_col_data.empty:
            ax.text(0.5, 0.5, "No positive data", ha="center", va="center", transform=ax.transAxes)
            log_col_data_for_plot = pd.Series(dtype=float)
        else:
            log_col_data_for_plot = np.log10(positive_col_data)
            hist_counts, _, _ = ax.hist(
                log_col_data_for_plot, bins=bins,
                edgecolor="black", alpha=0.9, rwidth=0.9
            )
            if hist_counts.size > 0:
                max_y_val = max(max_y_val, hist_counts.max())

        ax.set_title(col_name)
        ax.tick_params(axis="both", which="major", labelsize=10, labelbottom=True, labelleft=True, bottom=True, left=True)

        if col_name == initial_time_point_col:
            if cutoff_val > 0:
                log_cutoff = np.log10(cutoff_val)
                current_xlim = ax.get_xlim()
                if current_xlim[0] < current_xlim[1] and log_cutoff >= current_xlim[0] and log_cutoff <= current_xlim[1]:
                    ax.axvline(log_cutoff, color=COLORS[0], linestyle="--", label=f"Cutoff = {cutoff_val:.2g}")
                elif current_xlim[0] >= current_xlim[1]:
                    ax.axvline(log_cutoff, color=COLORS[0], linestyle="--", label=f"Cutoff = {cutoff_val:.2g}")
                else:
                    logger.info(f"Cutoff value {cutoff_val} (log10: {log_cutoff:.2f}) for {col_name} is outside plot x-limits")

                if ax.has_data():
                    ax.legend(frameon=False)
            else:
                logger.warning(f"Cutoff value ({cutoff_val}) for {initial_time_point_col} is not positive")

    # Add statistics text box
    stat_text_lines = [
        f"File: {filename}",
        f"Initial Time Point: '{initial_time_point_col}'",
        f"Cutoff Applied: >= {cutoff_val:.2g}",
        f"Original Rows: {stats['original_rows']:,}",
        f"Rows Kept: {stats['rows_kept'] if isinstance(stats['rows_kept'], str) else format(stats['rows_kept'], ',')} ({stats['percentage_rows_kept'] if isinstance(stats['percentage_rows_kept'], str) else format(stats['percentage_rows_kept'], '.1f')}%)",
        f"Original Counts: {stats['original_counts']:,}",
        f"Counts Kept: {stats['count_kept'] if isinstance(stats['count_kept'], str) else format(stats['count_kept'], ',')} ({stats['percentage_count_kept'] if isinstance(stats['percentage_count_kept'], str) else format(stats['percentage_count_kept'], '.1f')}%)",
    ]
    stat_text = "\n".join(stat_text_lines)

    fig.text(0.5, 0.95, stat_text, transform=fig.transFigure,
             verticalalignment="top", horizontalalignment="left",
             bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.7, ec="gray"))

    return fig

@logger.catch
def plot_distributions_and_calculate_stats(
    df: pd.DataFrame,
    filename: str,
    initial_time_point_col: str,
    cutoff_val: float,
    bins: int,
) -> tuple[plt.Figure, dict[str, Any]]:
    """Main analysis function to create plots and calculate filtering statistics."""
    stats = calculate_cutoff_statistics(df, initial_time_point_col, cutoff_val)
    fig = create_distribution_plot(df, filename, initial_time_point_col, cutoff_val, bins, stats)
    return fig, stats

@logger.catch
def display_summary_table(aggregated_stats: list[dict[str, Any]], total_time: float, headers: dict[str, str]) -> None:
    """Display formatted summary table of processing results and statistics."""
    if not aggregated_stats:
        logger.info("No statistics to display.")
        return

    table_data = []
    for file_stats in aggregated_stats:
        row = []
        if "error" in file_stats and len(file_stats) <= 2:
            row.append(file_stats.get("filename", "Unknown File"))
            error_message = file_stats.get("error", "Processing Error")
            for i, key in enumerate(headers.keys()):
                if i == 0:
                    continue
                if key == "original_rows":
                    row.append(error_message[:50])
                else:
                    row.append("N/A")
        else:
            for key in headers:
                value = file_stats.get(key, "N/A")
                if isinstance(value, float) and ("percentage" in key or "%" in headers[key]):
                    row.append(f"{value:.2f}%")
                elif isinstance(value, (int)) and not ("percentage" in key or "%" in headers[key]):
                    row.append(f"{value:,}")
                elif isinstance(value, float) and not ("percentage" in key or "%" in headers[key]):
                    row.append(f"{value:.2f}")
                else:
                    row.append(str(value))
        table_data.append(row)

    logger.info("\n--- Processing Summary ---")
    try:
        table_str = tabulate(table_data, headers=list(headers.values()), tablefmt="grid", stralign="left", numalign="right")
        for line in table_str.split("\n"):
            logger.info(line)
    except Exception as e:
        logger.error(f"Could not generate summary table with tabulate: {e}")
        logger.info("Raw aggregated stats:")
        for stat_item in aggregated_stats:
            logger.info(str(stat_item))

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Analyze read count distributions from TSV files and apply cutoffs.")
    parser.add_argument("-i", "--input", nargs="+", type=Path, required=True, help="One or more input TSV files.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path for the plots.")
    parser.add_argument("-t", "--initial_time_point", required=True, type=str, help="Name of the column representing the initial time point for cutoff application.")
    parser.add_argument("-c", "--cutoff", required=True, type=float, help="Cutoff value to apply to the initial time point column (values >= cutoff are kept).")
    parser.add_argument("--bins", type=int, default=50, help="Number of bins for histograms (default: %(default)s).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

def main() -> int:
    """Main entry point for read count distribution analysis with cutoff application."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = ReadCountDistributionAnalysisConfig(
            input_files=args.input,
            output_path=args.output,
            initial_time_point=args.initial_time_point,
            cutoff=args.cutoff,
            bins=args.bins,
        )

        logger.info(f"Starting processing of {len(config.input_files)} input files")
        logger.info(f"Initial time point column: '{config.initial_time_point}'")
        logger.info(f"Cutoff value: {config.cutoff}")
        logger.info(f"Histogram bins: {config.bins}")
        logger.info(f"Output PDF: {config.output_path}")

        start_time = time.time()
        all_file_stats: list[dict[str, Any]] = []

        with PdfPages(config.output_path) as pdf:
            for file_path in sorted(config.input_files, key=lambda p: p.name):
                filename = file_path.name
                logger.info(f"--- Processing file: {filename} ---")

                try:
                    df = load_and_validate_data(file_path)
                    fig, current_stats = plot_distributions_and_calculate_stats(
                        df, filename, config.initial_time_point, config.cutoff, config.bins
                    )

                    current_stats["filename"] = filename
                    all_file_stats.append(current_stats)

                    if fig:
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)
                        logger.info(f"Plot generated for {filename}.")
                        logger.debug(f"Detailed stats for {filename}: {current_stats}")
                    else:
                        logger.info(f"Plotting skipped for {filename} (no numeric data).")

                except FileNotFoundError:
                    logger.error(f"File {filename} not found. Skipping.")
                    all_file_stats.append({"filename": filename, "error": "File not found"})
                except pd.errors.EmptyDataError:
                    logger.warning(f"File {filename} is empty. Skipping.")
                    all_file_stats.append({"filename": filename, "error": "Empty file"})
                except pd.errors.ParserError as pe:
                    logger.error(f"ParserError for {filename}: {pe}. Skipping.")
                    all_file_stats.append({"filename": filename, "error": f"Parsing failed: {pe}"})
                except ValueError as ve:
                    logger.warning(f"ValueError for {filename}: {ve}. Skipping.")
                    all_file_stats.append({"filename": filename, "error": str(ve)})
                except Exception as e:
                    logger.error(f"An unexpected error occurred while processing {filename}: {e}")
                    all_file_stats.append({"filename": filename, "error": str(e)})

        end_time = time.time()
        total_processing_time = end_time - start_time

        logger.info(f"--- Analysis complete for all files. PDF saved to {config.output_path} ---")
        display_summary_table(all_file_stats, total_processing_time, SUMMARY_HEADERS)

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
