"""
Filter insertion reads by hard filtering based on read count thresholds at initial timepoints.

Filters insertion reads based on minimum read count thresholds at initial timepoints,
supporting reproducible scientific analysis with comprehensive logging and validation.
"""

import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple
from pydantic import BaseModel, Field, field_validator
import pandas as pd


# =============================== Configuration & Models ===============================

class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths and filtering parameters."""
    input_file: Path = Field(..., description="Path to the input TSV file with insertion reads")
    output_file: Path = Field(..., description="Path to the output TSV file for filtered reads")
    initial_timepoint: str = Field(..., description="Initial timepoint column name for filtering")
    cutoff_threshold: int = Field(..., ge=0, description="Minimum read count threshold")

    @field_validator('input_file')
    def validate_input_file(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        if not v.suffix.lower() in ['.tsv', '.txt']:
            logger.warning(f"Input file may not be TSV format: {v.suffix}")
        return v
    
    @field_validator('output_file')
    def validate_output_path(cls, v: Path) -> Path:
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    @field_validator('initial_timepoint')
    def validate_timepoint(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("Initial timepoint cannot be empty")
        return v.strip()

    class Config:
        frozen = True


class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the filtering analysis."""
    total_insertions: int = Field(..., ge=0, description="Total insertions before filtering")
    retained_insertions: int = Field(..., ge=0, description="Insertions retained after filtering")
    removed_insertions: int = Field(..., ge=0, description="Insertions removed by filtering")
    retention_rate: float = Field(..., ge=0, le=100, description="Percentage of insertions retained")
    samples_processed: int = Field(..., ge=0, description="Number of samples processed")


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
def load_insertion_data(input_file: Path) -> pd.DataFrame:
    """Load insertion reads from TSV file with multi-index structure."""
    logger.info(f"Loading insertion data from: {input_file}")
    
    try:
        df = pd.read_csv(
            input_file, 
            sep="\t", 
            index_col=[0, 1, 2, 3], 
            header=[0, 1]
        )
        logger.success(f"Loaded {df.shape[0]:,} insertions with {df.shape[1]} columns")
        
        if not isinstance(df.columns, pd.MultiIndex):
            raise ValueError("Expected MultiIndex columns with (Sample, Timepoint) structure")
        
        if len(df.columns.levels) != 2:
            raise ValueError("Expected columns to have 2 levels: Sample and Timepoint")
            
        return df
        
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        raise


@logger.catch
def validate_timepoint_exists(df: pd.DataFrame, timepoint: str) -> None:
    """Validate that the specified timepoint exists in the data."""
    available_timepoints = df.columns.get_level_values(1).unique()
    if timepoint not in available_timepoints:
        logger.error(f"Timepoint '{timepoint}' not found in data")
        logger.error(f"Available timepoints: {list(available_timepoints)}")
        raise ValueError(f"Invalid timepoint: {timepoint}")
    
    logger.debug(f"Validated timepoint '{timepoint}' exists in data")


@logger.catch
def apply_hard_filtering(df: pd.DataFrame, config: InputOutputConfig) -> Tuple[pd.DataFrame, AnalysisResult]:
    """Apply hard filtering across all samples based on read count threshold."""
    logger.info("Starting hard filtering process...")
    
    # Validate timepoint exists
    validate_timepoint_exists(df, config.initial_timepoint)
    
    # Display initial data info
    logger.info("=" * 60)
    logger.info("INITIAL DATA SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total insertions: {df.shape[0]:,}")
    logger.info(f"Total samples: {len(df.columns.get_level_values(0).unique())}")
    logger.info(f"Timepoints: {list(df.columns.get_level_values(1).unique())}")
    logger.info(f"Initial timepoint: {config.initial_timepoint}")
    logger.info(f"Cutoff threshold: {config.cutoff_threshold}")
    
    # Process each sample
    filtered_samples = {}
    sample_results = []
    total_insertions = 0
    retained_insertions = 0
    
    for sample_name, sample_data in df.groupby(level="Sample", axis=1):
        logger.debug(f"Processing sample: {sample_name}")
        
        sample_total = len(sample_data)
        mask = sample_data[(sample_name, config.initial_timepoint)] >= config.cutoff_threshold
        sample_retained = mask.sum()
        
        retention_rate = (sample_retained / sample_total * 100 
                         if sample_total > 0 else 0)
        
        total_insertions += sample_total
        retained_insertions += sample_retained
        
        logger.debug(
            f"Sample {sample_name}: {sample_retained:,}/{sample_total:,} "
            f"insertions retained ({retention_rate:.2f}%)"
        )
        
        if sample_retained > 0:
            filtered_samples[sample_name] = sample_data[mask]
        
        sample_results.append({
            'sample_name': sample_name,
            'total_insertions': sample_total,
            'retained_insertions': sample_retained,
            'retention_rate': retention_rate
        })
    
    # Combine filtered samples
    if not filtered_samples:
        logger.warning("No samples retained any insertions after filtering")
        filtered_df = pd.DataFrame(columns=df.columns)
    else:
        filtered_df = pd.concat(filtered_samples.values(), axis=1)
    
    # Calculate overall statistics
    removed_insertions = total_insertions - retained_insertions
    retention_rate = (retained_insertions / total_insertions * 100 
                     if total_insertions > 0 else 0)
    
    result = AnalysisResult(
        total_insertions=total_insertions,
        retained_insertions=retained_insertions,
        removed_insertions=removed_insertions,
        retention_rate=retention_rate,
        samples_processed=len(sample_results)
    )
    
    return filtered_df, result


# =============================== Main Function ===============================

def parse_arguments():
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter insertion reads by hard filtering based on read count thresholds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reads_hard_filtering.py -i raw_reads.tsv -o filtered_reads.tsv -itp 0h -c 5
  python reads_hard_filtering.py -i counts.tsv -o filtered.tsv --init-timepoint YES0 --cutoff 10
        """
    )
    
    parser.add_argument(
        "-i", "--input", 
        type=Path, 
        required=True, 
        help="Input TSV file with insertion reads"
    )
    parser.add_argument(
        "-o", "--output", 
        type=Path, 
        required=True, 
        help="Output TSV file for filtered reads"
    )
    parser.add_argument(
        "-itp", "--init-timepoint", 
        type=str, 
        required=True, 
        help="Initial timepoint column name for filtering"
    )
    parser.add_argument(
        "-c", "--cutoff", 
        type=int, 
        required=True, 
        help="Minimum read count threshold at initial timepoint"
    )
    parser.add_argument(
        "-v", "--verbose", 
        action="store_true", 
        help="Enable verbose logging"
    )
    
    return parser.parse_args()


@logger.catch
def main():
    """Main entry point of the script."""
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    try:
        # Validate input and output paths using the Pydantic model
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output,
            initial_timepoint=args.init_timepoint,
            cutoff_threshold=args.cutoff
        )

        logger.info(f"Starting processing of {config.input_file}")
        
        # Load data
        df = load_insertion_data(config.input_file)
        
        # Apply filtering
        filtered_df, results = apply_hard_filtering(df, config)
        
        # Save results
        logger.info("Saving filtered results...")
        filtered_df.to_csv(config.output_file, sep="\t", header=True, index=True)
        
        # Display summary
        logger.info("=" * 70)
        logger.info("FILTERING SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Total insertions: {results.total_insertions:,}")
        logger.info(f"Retained insertions: {results.retained_insertions:,}")
        logger.info(f"Removed insertions: {results.removed_insertions:,}")
        logger.success(f"Retention rate: {results.retention_rate:.2f}%")
        logger.info(f"Samples processed: {results.samples_processed}")
        
        logger.success(f"Analysis complete. Results saved to {config.output_file}")

    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()