"""
Insertion-level depletion analysis for non-replicates data.

This script performs insertion-level depletion analysis on transposon sequencing data
without replicates. It loads insertion counts and control insertions, performs median
normalization, calculates log-fold changes (LFC), and generates MA plots for
visualization.

Typical Usage:
    python insertion_level_depletion_analysis_no_replicates.py \
        --counts_file counts.tsv \
        --control_insertions_file controls.tsv \
        --init_timepoint T0 \
        --output_LFC_file lfc_results.tsv

Input: TSV files containing insertion counts and control insertions
Output: TSV file with LFC values and PDF with MA plots
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================== Constants ===============================
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Configuration for insertion-level depletion analysis."""
    counts_file: Path = Field(..., description="Path to the counts file")
    control_insertions_file: Path = Field(..., description="Path to the control insertions file")
    init_timepoint: str = Field(..., description="Initial timepoint")
    output_LFC_file: Path = Field(..., description="Path to the output LFC file")

    @field_validator('counts_file', 'control_insertions_file')
    def validate_input_files(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_LFC_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True

class AnalysisResult(BaseModel):
    """Results container for insertion-level depletion analysis."""
    status: str = Field(..., description="Status of the analysis")
    message: str = Field(..., description="Completion message")


# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# =============================== Core Functions ===============================
@logger.catch
def load_and_preprocess_data(
    counts_file: Path, control_insertions_file: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
) -> Tuple[pd.DataFrame, pd.DataFrame]:
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

# =============================== Main Function ===============================
def parse_arguments():
    """Parse command line arguments for depletion analysis."""
    parser = argparse.ArgumentParser(description="Insertion-level depletion analysis for non-replicates data")
    parser.add_argument("-i", "--counts_file", type=Path, required=True, help="Path to the counts file")
    parser.add_argument("-c", "--control_insertions_file", type=Path, required=True, help="Path to the control insertions file")
    parser.add_argument("-t", "--init_timepoint", type=str, required=True, help="Initial timepoint")
    parser.add_argument("-o", "--output_LFC_file", type=Path, required=True, help="Path to the output LFC file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """
    Main entry point for insertion-level depletion analysis.
    
    Orchestrates the complete workflow: data loading, normalization, 
    LFC calculation, and visualization generation.
    """
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    try:
        config = InputOutputConfig(
            counts_file=args.counts_file,
            control_insertions_file=args.control_insertions_file,
            init_timepoint=args.init_timepoint,
            output_LFC_file=args.output_LFC_file
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

        MA_plot_path = config.output_LFC_file.parent / "MA_plot.pdf"
        generate_MA_plots(M_values, A_values, MA_plot_path)
        logger.info(f"MA plots saved to {MA_plot_path}")

        result = AnalysisResult(
            status="completed",
            message="Insertion-level depletion analysis completed successfully"
        )
        
        logger.success(f"Analysis complete. Results saved to {config.output_LFC_file}")
        logger.success(f"MA plots saved to {MA_plot_path}")
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()