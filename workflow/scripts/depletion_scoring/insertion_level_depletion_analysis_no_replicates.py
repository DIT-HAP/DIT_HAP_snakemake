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
# ]
# ///

"""
Insertion-Level Depletion Analysis (No Replicates)
==================================================

Perform insertion-level depletion analysis on transposon sequencing data that
lacks replicates. The script loads insertion counts and a set of control
insertions, normalises counts against the control-insertion medians, computes
log-fold changes (LFC / M values) relative to an initial timepoint, and
generates MA plots for visual quality control.

Input
-----
- Counts TSV: 4 index columns and a two-level column header (sample, timepoint),
  read via ``--counts_file`` / ``-i``.
- Control insertions TSV: 4 index columns, read via
  ``--control_insertions_file`` / ``-c``.
- Initial timepoint label, passed via ``--init_timepoint`` / ``-t``.

Output
------
- LFC table (M values) written to ``--output_LFC_file`` / ``-o``.
- ``normed_counts.tsv`` and ``baseMean.tsv`` written alongside the LFC file.
- ``MA_plot.pdf`` with one MA panel per timepoint.

Usage
-----
    python insertion_level_depletion_analysis_no_replicates.py \
        -i counts.tsv -c controls.tsv -t T0 -o LFC.tsv
    python insertion_level_depletion_analysis_no_replicates.py \
        -i counts.tsv -c controls.tsv -t T0 -o LFC.tsv --verbose

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
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Input/output paths and parameters for insertion-level depletion analysis."""
    counts_file: Path
    control_insertions_file: Path
    init_timepoint: str
    output_LFC_file: Path

    def __post_init__(self) -> None:
        for path in (self.counts_file, self.control_insertions_file):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        self.output_LFC_file.parent.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Results container for insertion-level depletion analysis."""
    status: str
    message: str

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

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_and_preprocess_data(
    counts_file: Path, control_insertions_file: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load counts and control insertions data from TSV files."""
    counts_df = pd.read_csv(
        counts_file, sep="\t", index_col=[0, 1, 2, 3], header=[0, 1]
    )
    control_insertions_df = pd.read_csv(
        control_insertions_file, sep="\t", index_col=[0, 1, 2, 3]
    )
    return counts_df, control_insertions_df


@logger.catch
def perform_median_normalization(
    counts_df: pd.DataFrame, control_insertions_df: pd.DataFrame
) -> pd.DataFrame:
    """Normalize counts using median values from control insertions."""
    median_values = counts_df.loc[control_insertions_df.index].median()
    min_median_values = median_values.min()
    normalized_counts = counts_df.mul(min_median_values).div(median_values)
    return normalized_counts


@logger.catch
def calculate_MA_values(
    normalized_counts: pd.DataFrame, init_timepoint: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate M (log-fold change) and A (average abundance) values for MA plots."""
    M_values = (
        -(normalized_counts + 1)
        .div((normalized_counts.xs(init_timepoint, level=1, axis=1) + 1), axis=0)
        .map(np.log2)
    )

    A_values = (normalized_counts + 1).mul(
        (normalized_counts.xs(init_timepoint, level=1, axis=1) + 1), axis=0
    ).map(np.log2) * 0.5

    return M_values, A_values


@logger.catch
def write_ma_values_tsv(M_values: pd.DataFrame, A_values: pd.DataFrame, output_path: Path) -> None:
    """Write M and A values to a long-format TSV for the rendering layer."""
    # M_values and A_values have MultiIndex columns (Sample, Timepoint)
    # We need the timepoint level only
    timepoints = M_values.columns.get_level_values(1).unique().tolist()

    rows = []
    for timepoint in timepoints:
        # Select all columns with this timepoint
        M_cols = [col for col in M_values.columns if col[1] == timepoint]
        A_cols = [col for col in A_values.columns if col[1] == timepoint]

        for M_col, A_col in zip(M_cols, A_cols):
            M_data = M_values[M_col].values
            A_data = A_values[A_col].values

            for m_val, a_val in zip(M_data, A_data):
                rows.append({
                    'timepoint': timepoint,
                    'M': m_val,
                    'A': a_val,
                })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False, float_format='%.3f')


@logger.catch
def generate_MA_plots(M_values: pd.DataFrame, A_values: pd.DataFrame, output_path: Path):
    """Generate MA plots for each timepoint and save as PDF."""
    timepoints = M_values.columns.tolist()
    n_rows = len(timepoints)

    fig, ax = plt.subplots(n_rows, 1, figsize=(AX_WIDTH, AX_HEIGHT * n_rows), sharex=True, sharey=True)

    for row, timepoint in enumerate(timepoints):
        M_data = M_values[timepoint]
        A_data = A_values[timepoint]

        ax[row].scatter(
            M_data,
            A_data,
            s=10,
            facecolor="none",
            edgecolor="black",
            alpha=0.5,
            rasterized=True,
        )
        ax[row].axvline(0, c="r", ls="--", lw=2, alpha=0.5)
        ax[row].set_xlabel("M value")
        ax[row].set_ylabel("A value")
        ax[row].set_title(f"MA plot - {timepoint}")

    fig.savefig(output_path)
    plt.close()

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command line arguments for depletion analysis."""
    parser = argparse.ArgumentParser(description="Insertion-level depletion analysis for non-replicates data")
    parser.add_argument("-i", "--counts_file", type=Path, required=True, help="Path to the counts file")
    parser.add_argument("-c", "--control_insertions_file", type=Path, required=True, help="Path to the control insertions file")
    parser.add_argument("-t", "--init_timepoint", type=str, required=True, help="Initial timepoint")
    parser.add_argument("-o", "--output_LFC_file", type=Path, required=True, help="Path to the output LFC file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point: load, normalize, compute LFC, and render MA plots."""
    args = parse_args()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(log_level)

    try:
        config = InputOutputConfig(
            counts_file=args.counts_file,
            control_insertions_file=args.control_insertions_file,
            init_timepoint=args.init_timepoint,
            output_LFC_file=args.output_LFC_file,
        )

        logger.info(f"Starting insertion-level depletion analysis for {config.counts_file}")

        counts_df, control_insertions_df = load_and_preprocess_data(
            config.counts_file, config.control_insertions_file
        )
        logger.info("Data loaded successfully")

        normalized_counts = perform_median_normalization(counts_df, control_insertions_df)
        logger.info("Median-based normalization completed")

        M_values, A_values = calculate_MA_values(normalized_counts, config.init_timepoint)
        logger.info("M and A values calculated")

        M_values.droplevel(0, axis=1).to_csv(config.output_LFC_file, sep="\t", index=True, float_format="%.3f")
        logger.info(f"LFC results saved to {config.output_LFC_file}")

        normalized_counts.droplevel(0, axis=1).to_csv(config.output_LFC_file.parent / "normed_counts.tsv", sep="\t", index=True, float_format="%.3f")
        logger.info(f"Normalized counts saved to {config.output_LFC_file.parent / 'normed_counts.tsv'}")

        baseMean = normalized_counts.droplevel(0, axis=1).copy()
        for col in baseMean.columns:
            baseMean[col] = baseMean[config.init_timepoint]
        baseMean.to_csv(config.output_LFC_file.parent / "baseMean.tsv", sep="\t", index=True, float_format="%.3f")
        logger.info(f"Base mean saved to {config.output_LFC_file.parent / 'baseMean.tsv'}")

        # Write MA values TSV for rendering layer
        ma_values_path = config.output_LFC_file.parent / "ma_values.tsv"
        write_ma_values_tsv(M_values, A_values, ma_values_path)
        logger.info(f"MA values saved to {ma_values_path}")

        logger.success(f"Analysis complete. Results saved to {config.output_LFC_file}")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
