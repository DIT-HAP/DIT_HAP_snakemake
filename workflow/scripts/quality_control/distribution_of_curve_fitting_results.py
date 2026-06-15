"""
The python script template for the DIT-HAP project.

This script reads a CSV file with statistical data and generates histogram plots
for all numeric columns, arranged in a 4-column subplot layout and saved as PDF.

Typical Usage:
    python distribution_of_curve_fitting_results.py --input input.csv --output output.pdf

Input: CSV file with statistical data
Output: PDF file with histogram plots
Other information: Uses matplotlib for plotting with DIT_HAP styling
"""

# =============================== Imports ===============================
# Standard library imports
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple
import time

# Third-party imports
from loguru import logger
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd

# Plotting imports
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# =============================== Constants ===============================
# The following is for plotting
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_file: Path = Field(..., description="Path to the input file")
    output_file: Path = Field(..., description="Path to the output file")
    bins: int = Field(30, description="Number of bins for histograms")

    @field_validator('input_file')
    def validate_input_file(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True) # Create dir if it doesn't exist
        return v
    
    class Config:
        frozen = True # Makes the model immutable after creation

class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    total_rows: int = Field(..., ge=0, description="Total number of data rows")
    total_columns: int = Field(..., ge=0, description="Total number of columns")
    numeric_columns: int = Field(..., ge=0, description="Number of numeric columns analyzed")
    plots_generated: int = Field(..., ge=0, description="Number of histogram plots created")
    missing_values: int = Field(..., ge=0, description="Total missing values in dataset")
    execution_time: float = Field(..., ge=0.0, description="Script execution time in seconds")


# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove() # Remove default logger
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# =============================== Core Functions ===============================
@logger.catch
def load_and_analyze_data(input_file: Path) -> Tuple[pd.DataFrame, List[str], Dict[str, any]]:
    """Load CSV data and identify numeric columns for analysis."""
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
def create_histogram_plots(df: pd.DataFrame, numeric_columns: List[str], 
                         output_file: Path, bins: int = 30) -> Dict[str, any]:
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
    """Core analysis function for processing curve fitting results and generating distribution plots."""
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


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(
        description="Generate histogram plots for numeric columns in curve fitting results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input CSV file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output PDF file")
    parser.add_argument("--bins", type=int, default=30, help="Number of bins for histograms")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point for generating distribution plots of curve fitting results."""

    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
    try:
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output,
            bins=args.bins
        )

        logger.info(f"Starting processing of {config.input_file}")

        # Run the core analysis/logic
        results = analyze_curve_fitting_results(config)
        
        # Log completion
        logger.success(f"Analysis complete. Results saved to {config.output_file}")
        logger.info(f"Generated {results.plots_generated} plots from {results.numeric_columns} numeric columns")
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()