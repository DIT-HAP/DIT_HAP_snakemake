"""
PBL-PBR correlation analysis script for the DIT-HAP project.

This script analyzes the correlation between PBL and PBR values from multiple TSV files,
generating scatter plots with regression lines and statistical summaries.

Typical Usage:
    python PBL_PBR_correlation_analysis.py --input file1.tsv file2.tsv --output results.pdf

Input: TSV files containing PBL and PBR columns with multi-index structure
Output: PDF file containing correlation plots and statistical analysis
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

SCRIPT_DIR = Path(__file__).parent.resolve()
TARGET_path = str((SCRIPT_DIR / "../../src").resolve())
sys.path.append(TARGET_path)
from plot import create_scatter_correlation_plot

# =============================== Constants ===============================
STYLE_path = str((SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle").resolve())
plt.style.use(STYLE_path)
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================== Configuration & Models ===============================
class PBL_PBR_CorrelationAnalysisConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths for PBL-PBR correlation analysis."""
    input_files: List[Path] = Field(..., description="List of input TSV files")
    output_path: Path = Field(..., description="Path for output PDF file")

    @field_validator('input_files')
    def validate_input_files(cls, v):
        """Validate that all input files exist."""
        for file_path in v:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        return v
    
    @field_validator('output_path')
    def validate_output_path(cls, v):
        """Validate output path and create directory if needed."""
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
def read_tsv_file(file_path: Path) -> Optional[pd.DataFrame]:
    """Read TSV file and validate required columns."""
    logger.info(f"Reading TSV file: {file_path}")

    df = pd.read_csv(file_path, sep='\t', index_col=[0,1,2])
    
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
    """Create correlation plot for a single file with statistics."""
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

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="Analyze correlation between PBL and PBR from multiple TSV files")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="Input TSV files (space-separated)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output PDF file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script."""
    
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
    try:
        config = PBL_PBR_CorrelationAnalysisConfig(
            input_files=args.input,
            output_path=args.output
        )

        logger.info("=== PBL-PBR Correlation Analysis ===")
        logger.info(f"Processing {len(config.input_files)} input files...")
        
        # Sort files by name
        sorted_files = sorted(config.input_files, key=lambda x: x.name)
        logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")
        
        # Read and process files
        data_dict = {}
        
        for file_path in sorted_files:
            filename = file_path.name
            logger.info(f"Reading {filename}...")
            
            df = read_tsv_file(file_path)
            if df is not None:
                data_dict[filename] = df
        
        if not data_dict:
            logger.error("Error: No valid data found in any input file!")
            sys.exit(1)
        
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
            sys.exit(1)
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
