"""
The python script for merging PBL and PBR insertion counts by coordinate and strand.

This script merges insertion data from PBL (left primer) and PBR (right primer) reads,
organizing by chromosome, coordinate, and strand for downstream analysis.

Typical Usage:
    python merge_strand_insertions.py --inputPBL pbl_insertions.tsv --inputPBR pbr_insertions.tsv --output merged_insertions.tsv

Input: PBL and PBR insertion TSV files
Output: Merged insertion TSV file
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple
from pydantic import BaseModel, Field, field_validator
import pandas as pd


# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_pbl: Path = Field(..., description="Path to the PBL input file")
    input_pbr: Path = Field(..., description="Path to the PBR input file")
    output_file: Path = Field(..., description="Path to the output file")

    @field_validator('input_pbl', 'input_pbr')
    def validate_input_files(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True) # Create dir if it doesn't exist
        return v
    
    class Config:
        frozen = True # Makes the model immutable after creation

class MergeResult(BaseModel):
    """Pydantic model to hold and validate the results of the merge operation."""
    total_sites_processed: int = Field(..., ge=0, description="Total number of insertion sites processed")
    total_reads_merged: int = Field(..., ge=0, description="Total number of reads merged")
    coordinate_strand_pairs: int = Field(..., ge=0, description="Number of unique coordinate-strand pairs")


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
def merge_insertion_data(pbl_df: pd.DataFrame, pbr_df: pd.DataFrame) -> pd.DataFrame:
    """Merge PBL and PBR insertion data."""
    logger.info("Merging PBL and PBR insertion data")
    
    merged_df = pd.merge(
        pbl_df,
        pbr_df,
        how="outer",
        on=["Chr", "Coordinate"],
        suffixes=("_PBL", "_PBR"),
    ).fillna(0)
    
    # Create plus strand data
    plus_df = merged_df[["Chr", "Coordinate", "-_PBL", "+_PBR"]].copy()
    plus_df["Strand"] = "+"
    plus_df.rename(columns={"-_PBL": "PBL", "+_PBR": "PBR"}, inplace=True)
    
    # Create minus strand data
    minus_df = merged_df[["Chr", "Coordinate", "+_PBL", "-_PBR"]].copy()
    minus_df["Strand"] = "-"
    minus_df.rename(columns={"+_PBL": "PBL", "-_PBR": "PBR"}, inplace=True)
    
    # Combine and finalize
    final_df = pd.concat([plus_df, minus_df], axis=0)
    final_df = final_df.set_index(["Chr", "Coordinate", "Strand"])
    final_df = final_df.astype(int).sort_index()
    final_df["Reads"] = final_df["PBL"] + final_df["PBR"]
    
    logger.success(f"Merged data: {len(final_df):,} coordinate-strand pairs")
    return final_df


@logger.catch
def save_merged_data(merged_df: pd.DataFrame, output_path: Path) -> None:
    """Save merged data to TSV file."""
    logger.info(f"Writing merged data to: {output_path}")
    merged_df.to_csv(output_path, sep='\t', index=True, header=True)
    
    total_reads = merged_df["Reads"].sum()
    pbl_reads = merged_df["PBL"].sum()
    pbr_reads = merged_df["PBR"].sum()
    
    logger.info(f"Total reads: {total_reads:,}")
    logger.info(f"PBL reads: {pbl_reads:,} ({pbl_reads/total_reads*100:.1f}%)")
    logger.info(f"PBR reads: {pbr_reads:,} ({pbr_reads/total_reads*100:.1f}%)")
    logger.success(f"Output saved to: {output_path}")


@logger.catch 
def main_processing_function(config: InputOutputConfig) -> MergeResult:
    """
    Main processing function to merge strand-specific insertions.
    """
    # Load data
    pbl_data = pd.read_csv(config.input_pbl, sep="\t", header=0)
    pbr_data = pd.read_csv(config.input_pbr, sep="\t", header=0)
    
    # Merge data
    merged_df = merge_insertion_data(pbl_data, pbr_data)
    
    # Save results
    save_merged_data(merged_df, config.output_file)
    
    # Calculate results summary
    total_sites_processed = len(merged_df)
    total_reads_merged = merged_df["Reads"].sum()
    coordinate_strand_pairs = len(merged_df)
    
    # Create result object
    result = MergeResult(
        total_sites_processed=total_sites_processed,
        total_reads_merged=total_reads_merged,
        coordinate_strand_pairs=coordinate_strand_pairs
    )
    
    return result


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="Merge PBL and PBR insertion counts by coordinate and strand")
    parser.add_argument("-i", "--inputPBL", type=Path, required=True, help="Path to the PBL input file")
    parser.add_argument("-j", "--inputPBR", type=Path, required=True, help="Path to the PBR input file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script. Replace this docstring with a relevant one."""

    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
    try:
        config = InputOutputConfig(
            input_pbl=args.inputPBL,
            input_pbr=args.inputPBR,
            output_file=args.output
        )

        logger.info(f"Starting processing of PBL: {config.input_pbl} and PBR: {config.input_pbr}")

        # Run the core analysis/logic. REPLACE THIS WITH YOUR LOGIC.
        results = main_processing_function(config)
        
        # Log results summary
        logger.info("Processing completed:")
        logger.info(f"  - Total sites processed: {results.total_sites_processed:,}")
        logger.info(f"  - Total reads merged: {results.total_reads_merged:,}")
        logger.info(f"  - Coordinate-strand pairs: {results.coordinate_strand_pairs:,}")
        
        logger.success(f"Analysis complete. Results saved to {config.output_file}")
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()