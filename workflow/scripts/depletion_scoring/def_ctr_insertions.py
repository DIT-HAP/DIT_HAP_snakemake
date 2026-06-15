"""
Control insertion selection for transposon depletion analysis.

This module provides functionality to select control insertions from transposon
insertion data for use in depletion analysis. Control insertions are selected
based on stringent criteria to ensure they represent neutral genomic regions
unaffected by selection pressure.

The selection criteria for control insertions include:
1. Location in intergenic regions (non-genic, non-regulatory regions)
2. Minimum distance of 500 bp from any gene boundaries
3. Exclusion of regions with potential regulatory function

Typical Usage:
    python def_ctr_insertions.py --input insertion_counts.tsv --annotation genomic_annotations.tsv --output control_insertions.tsv

Input: Tab-separated insertion count table and genomic annotation table
Output: Tab-separated file containing selected control insertions
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Tuple
from loguru import logger
from pydantic import BaseModel, Field, field_validator
import pandas as pd

# =============================== Constants ===============================
CONTROL_DISTANCE_THRESHOLD = 500  # Minimum distance (bp) from gene boundaries for control insertions

# =============================== Configuration & Models ===============================


class InputOutputConfig(BaseModel):
    """Validate input and output file paths."""
    input_file: Path = Field(..., description="Path to the insertion table")
    annotation_file: Path = Field(..., description="Path to the annotation table")
    output_file: Path = Field(..., description="Path to the output file")

    @field_validator('input_file', 'annotation_file')
    def validate_input_file(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True


class ControlSelectionResult(BaseModel):
    """Hold and validate the results of the analysis."""
    total_insertions_processed: int = Field(..., ge=0, description="Total number of insertions processed")
    control_insertions_selected: int = Field(..., ge=0, description="Number of control insertions selected")
    success_rate: float = Field(..., ge=0.0, le=100.0, description="Percentage of successful operations")

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
def load_and_preprocess_data(counts_file: Path, annotations_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load and preprocess insertion count and genomic annotation tables."""
    counts_df = pd.read_csv(
        counts_file, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t"
    )
    
    # Remove rows with any NA value
    counts_df = counts_df.dropna(axis=0, how="any").copy()
    
    # Load and process annotations
    insertion_annotations = pd.read_csv(
        annotations_file, index_col=[0, 1, 2, 3], sep="\t"
    )
    
    return counts_df, insertion_annotations


@logger.catch
def get_control_insertions(counts_df: pd.DataFrame, insertion_annotations: pd.DataFrame) -> pd.DataFrame:
    """Select control insertions based on stringent genomic criteria."""
    ctr_insertions = insertion_annotations.query(
        f"Type == 'Intergenic region' and Distance_to_region_start > {CONTROL_DISTANCE_THRESHOLD} and Distance_to_region_end > {CONTROL_DISTANCE_THRESHOLD}"
    )
    
    ctr_insertions = ctr_insertions[ctr_insertions.index.isin(counts_df.index)].drop_duplicates(keep="first")
    
    return ctr_insertions


@logger.catch
def save_results(control_insertions: pd.DataFrame, output_file: Path) -> None:
    """Save selected control insertions to a tab-separated file."""
    control_insertions.to_csv(output_file, sep="\t", index=True)
    logger.info(f"Saved {len(control_insertions)} control insertions to {output_file}")


@logger.catch
def process_control_insertions(config: InputOutputConfig) -> ControlSelectionResult:
    """Execute the complete control insertion selection pipeline."""
    logger.info(f"Starting processing of {config.input_file}")
    
    # Load data
    counts_df, insertion_annotations = load_and_preprocess_data(
        config.input_file, config.annotation_file
    )
    
    total_insertions = len(counts_df)
    logger.info(f"Loaded {total_insertions} insertions for processing")
    
    # Select control insertions
    control_insertions = get_control_insertions(counts_df, insertion_annotations)
    control_count = len(control_insertions)
    
    # Save results
    save_results(control_insertions, config.output_file)
    
    # Calculate success metrics
    success_rate = (control_count / total_insertions * 100) if total_insertions > 0 else 0.0
    
    return ControlSelectionResult(
        total_insertions_processed=total_insertions,
        control_insertions_selected=control_count,
        success_rate=success_rate
    )

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Select control insertions for transposon depletion analysis"
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to tab-separated insertion count table")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Path to tab-separated genomic annotation table")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output path for selected control insertions")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
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
            annotation_file=args.annotation,
            output_file=args.output
        )
        
        # Run the core analysis
        results = process_control_insertions(config)
        
        logger.success(f"Analysis complete. Results: {results.model_dump_json()}")
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

# =============================== Final Reminders ===============================