"""
Insertion Orientation Analysis Script for the DIT-HAP project.

This script analyzes strand orientation (+/-) pairs from multiple TSV files with multi-level indexing.
For each file, it creates a single figure with subplots arranged in 1 row × n columns (where n = number 
of numeric columns). Each subplot shows all +/- strand pairs for that column using log-scale scatter plots. 
Results are saved as a multi-page PDF report with correlation statistics.

Typical Usage:
    python insertion_orientation_analysis.py --input file1.tsv file2.tsv --output orientation_analysis.pdf

Input: One or more TSV files with multi-level indexing where level 2 represents strand orientation (+/-)
Output: Multi-page PDF report with strand orientation analysis plots and correlation statistics TSV file
Other information: The script extracts +/- strand pairs and creates log-scale scatter plots with correlation analysis.
"""

# =============================== Imports ===============================
import sys
import argparse
import time
from pathlib import Path
from typing import List
import pandas as pd
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pydantic import BaseModel, Field, field_validator

SCRIPT_DIR = Path(__file__).parent.resolve()
TARGET_path = str((SCRIPT_DIR / "../../src").resolve())
sys.path.append(TARGET_path)
from plot import create_scatter_correlation_plot
from utils import read_file

# =============================== Constants ===============================
STYLE_path = str((SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(STYLE_path)
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================== Configuration & Models ===============================
class InsertionOrientationAnalysisConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths for insertion orientation analysis."""
    input_files: List[Path] = Field(..., description="List of input TSV files with multi-level indexing")
    output_path: Path = Field(..., description="Path for output PDF file")

    @field_validator('input_files')
    def validate_input_files(cls, v):
        if not v:
            raise ValueError("At least one input file must be provided")
        for file_path in v:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
            if not file_path.suffix.lower() in ['.tsv', '.txt']:
                raise ValueError(f"Input file must be a TSV file: {file_path}")
        return v
    
    @field_validator('output_path')
    def validate_output_path(cls, v):
        if not v.suffix.lower() == '.pdf':
            raise ValueError(f"Output file must be a PDF: {v}")
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True

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
def create_file_comparison_figure(
    df: pd.DataFrame,
    filename: str,
    output_path: Path,
) -> None:
    """Create a figure with subplots comparing +/- strand values for all columns."""
    
    plus_minus_pair = df.stack(future_stack=True).stack(future_stack=True).unstack("Strand").dropna(axis=0)
    timepoints = plus_minus_pair.index.get_level_values("Timepoint").unique()
    
    # Create figure with subplots: 1 row × n columns
    n_rows = len(timepoints)
    fig_width = AX_WIDTH
    fig_height = n_rows * AX_HEIGHT

    total_figures = 0
    with PdfPages(output_path) as pdf:
        for sample, sample_df in plus_minus_pair.groupby(level="Sample"):
    
            fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height))
            if n_rows == 1:
                axes = [axes]  # Ensure axes is always a list
            
            # Process each column
            for row_idx, (timepoint, sub_df) in enumerate(sample_df.groupby(level="Timepoint")):
                ax = axes[row_idx]

                filtered_sub_df = sub_df[sub_df.min(axis=1) > 0]

                pos_array = filtered_sub_df["+"].to_numpy()
                neg_array = filtered_sub_df["-"].to_numpy()

                create_scatter_correlation_plot(
                    x=pos_array,
                    y=neg_array,
                    ax=ax,
                    xscale='log',
                    yscale='log',
                )
                
                # Customize subplot
                ax.set_xlabel('Positive Strand (+)')
                if row_idx == 0:  # Only label y-axis on leftmost subplot
                    ax.set_ylabel('Negative Strand (-)')
                ax.set_title(f"Sample: {sample}\nTimepoint: {timepoint}")
                ax.grid(True)
            
            plt.tight_layout()
            # Save figure to PDF
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            total_figures += 1

    logger.info(f"Generated {total_figures} figures")

@logger.catch
def analyze_multiple_files(
    input_files: List[Path],
    output_path: Path,
) -> None:
    """Analyze strand orientations across multiple files and generate PDF report."""
    logger.info("Starting multi-file strand orientation analysis...")
    
    # Sort files by name for consistent processing order
    sorted_files = sorted(input_files, key=lambda p: p.name)
    logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")
    
    # Generate plots and save to PDF
    for file_path in sorted_files:
        filename = file_path.name
        logger.info(f"--- Processing file: {filename} ---")
        
        try:
            # Read the file
            df = read_file(file_path, **{"index_col": [0,1,2,3], "header": [0,1]})
            
            # Create figure for this file
            create_file_comparison_figure(
                df, filename, output_path
            )
            
            logger.info(f"Generated figure for {filename}")
            
        except Exception as e:
            logger.error(f"Failed to process {filename}: {e}", exc_info=True)

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="Analyze insertion orientation (+/-) strand pairs from multiple TSV files.")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="One or more input TSV files with multi-level indexing.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path for the plots.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script."""
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    config = InsertionOrientationAnalysisConfig(
        input_files=args.input,
        output_path=args.output
    )

    logger.info("=== Insertion Orientation Analysis ===")
    logger.info(f"Processing {len(config.input_files)} input files...")
    
    start_time = time.time()
    
    # Perform multi-file analysis
    analyze_multiple_files(config.input_files, config.output_path)
    
    end_time = time.time()
    total_time = end_time - start_time
    logger.info(f"Completed insertion orientation analysis in {total_time:.2f} seconds.")

if __name__ == "__main__":
    main()