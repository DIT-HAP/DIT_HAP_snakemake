"""
Extract transposon insertion sites from aligned read TSV files.

Processes TSV output from BAM parsing to identify and count transposon insertion
sites based on read alignment coordinates and strand orientation.

Typical Usage:
    python extract_insertion_sites.py --input_tsv input.tsv --output_tsv output.tsv

Input: TSV file with aligned read information
Output: TSV file with insertion site counts
"""

# =============================== Imports ===============================
import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_file: Path = Field(..., description="Path to the input TSV file")
    output_file: Path = Field(..., description="Path to the output TSV file")
    chunk_size: int = Field(500000, ge=10000, le=5000000, description="Rows per chunk")

    @field_validator('input_file')
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


class ExtractionStats(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    total_rows: int = Field(..., ge=0, description="Total number of rows processed")
    valid_rows: int = Field(..., ge=0, description="Number of valid rows processed")
    invalid_rows: int = Field(..., ge=0, description="Number of invalid/skipped rows")
    unique_sites: int = Field(..., ge=0, description="Number of unique insertion sites")
    total_plus_insertions: int = Field(..., ge=0, description="Total + strand insertions")
    total_minus_insertions: int = Field(..., ge=0, description="Total - strand insertions")

    @property
    def total_insertions(self) -> int:
        """Total number of insertions."""
        return self.total_plus_insertions + self.total_minus_insertions
    
    @property
    def validity_rate(self) -> float:
        """Percentage of valid rows."""
        if self.total_rows == 0:
            return 0.0
        return (self.valid_rows / self.total_rows) * 100


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
def calculate_insertion_coordinate(row: pd.Series) -> Optional[int]:
    """Calculate insertion coordinate based on strand orientation."""
    try:
        strand = row['R1_Strand']
        
        if strand == '+':
            # For + strand: TTAA[Genome] - use position after TTAA (ref_start + 4)
            return int(row['R1_Ref_Start']) + 4
        elif strand == '-':
            # For - strand: [Genome]TTAA - use position at end (ref_end)
            return int(row['R1_Ref_End'])
        else:
            return None
            
    except (ValueError, TypeError, KeyError):
        return None


@logger.catch
def create_validation_mask(df: pd.DataFrame) -> pd.Series:
    """Create a validation mask for filtering valid rows."""
    return (
        df['R1_Strand'].notna() & 
        df['R1_Chrom'].notna() & 
        df['R1_Ref_Start'].notna() & 
        df['R1_Ref_End'].notna() &
        df['R1_Strand'].isin(['+', '-'])
    )


@logger.catch
def count_insertions_vectorized(valid_df: pd.DataFrame) -> Dict[Tuple[str, int], Dict[str, int]]:
    """Count insertions using vectorized operations."""
    # Calculate coordinates for all valid rows at once
    coordinates = valid_df.apply(calculate_insertion_coordinate, axis=1)
    valid_coords_df = valid_df[coordinates.notna()].copy()
    valid_coords_df['Insertion_Coordinate'] = coordinates.dropna()
    
    # Group and count using pandas operations
    grouped = valid_coords_df.groupby(['R1_Chrom', 'Insertion_Coordinate', 'R1_Strand']).size()
    
    # Convert to our dictionary format
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    for (chrom, coord, strand), count in grouped.items():
        insertion_counts[(chrom, int(coord))][strand] = count
    
    return dict(insertion_counts)


@logger.catch
def extract_insertion_sites(chunk: pd.DataFrame, chunk_num: int) -> Tuple[Dict[Tuple[str, int], Dict[str, int]], int, int]:
    """Process a single chunk of data to extract insertion sites."""
    chunk_rows = len(chunk)
    
    # Filter valid rows
    valid_mask = create_validation_mask(chunk)
    valid_chunk = chunk[valid_mask].copy()
    valid_rows = len(valid_chunk)
    invalid_rows = chunk_rows - valid_rows
    
    if chunk_num == 1 or chunk_num % 10 == 0:
        retention_rate = (valid_rows / chunk_rows * 100) if chunk_rows > 0 else 0
        logger.info(f"Chunk {chunk_num}: {valid_rows:,}/{chunk_rows:,} valid rows ({retention_rate:.1f}%)")
    
    # Count insertions using vectorized operations
    insertion_counts = count_insertions_vectorized(valid_chunk) if valid_rows > 0 else {}
    
    return insertion_counts, valid_rows, invalid_rows


@logger.catch
def create_output_dataframe(insertion_counts: Dict[Tuple[str, int], Dict[str, int]]) -> Tuple[pd.DataFrame, int, int]:
    """Create output DataFrame from insertion counts."""
    output_data = []
    total_plus = 0
    total_minus = 0
    
    for (chrom, coord), strand_counts in insertion_counts.items():
        plus_count = strand_counts['+']
        minus_count = strand_counts['-']
        
        output_data.append({
            'Chr': chrom,
            'Coordinate': coord,
            '+': plus_count,
            '-': minus_count
        })
        
        total_plus += plus_count
        total_minus += minus_count
    
    output_df = pd.DataFrame(output_data)
    output_df = output_df.sort_values(['Chr', 'Coordinate'])
    
    return output_df, total_plus, total_minus


@logger.catch
def process_chunks(config: InputOutputConfig) -> Tuple[Dict[Tuple[str, int], Dict[str, int]], int, int, int, int]:
    """Process all chunks and aggregate results."""
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    total_rows = 0
    total_valid_rows = 0
    total_invalid_rows = 0
    chunk_count = 0
    
    logger.info("Starting chunked processing...")
    
    try:
        chunk_iterator = pd.read_csv(
            config.input_file,
            sep='\t',
            chunksize=config.chunk_size,
            na_values=['N/A', 'NA', '']
        )
        
        for chunk_df in chunk_iterator:
            chunk_count += 1
            total_rows += len(chunk_df)
            
            if chunk_count % 10 == 0:
                logger.info(f"Processing chunk {chunk_count}, total rows: {total_rows:,}")
            
            chunk_counts, valid_rows, invalid_rows = extract_insertion_sites(chunk_df, chunk_count)
            
            # Aggregate counts
            for key, strand_counts in chunk_counts.items():
                insertion_counts[key]['+'] += strand_counts['+']
                insertion_counts[key]['-'] += strand_counts['-']
            
            total_valid_rows += valid_rows
            total_invalid_rows += invalid_rows
        
        logger.success(f"Completed processing {chunk_count} chunks")
        
    except Exception as e:
        logger.error(f"Error during chunked processing: {e}")
        raise
    
    return dict(insertion_counts), total_rows, total_valid_rows, total_invalid_rows, chunk_count


@logger.catch
def write_empty_output(config: InputOutputConfig, total_rows: int, total_valid_rows: int, total_invalid_rows: int) -> ExtractionStats:
    """Write empty output file when no insertion sites found."""
    logger.warning("No insertion sites found!")
    empty_df = pd.DataFrame(columns=['Chr', 'Coordinate', '+', '-'])
    empty_df.to_csv(config.output_file, sep='\t', index=False)
    
    return ExtractionStats(
        total_rows=total_rows,
        valid_rows=total_valid_rows,
        invalid_rows=total_invalid_rows,
        unique_sites=0,
        total_plus_insertions=0,
        total_minus_insertions=0
    )


@logger.catch
def count_insertion_sites(config: InputOutputConfig) -> ExtractionStats:
    """Main function to extract insertion sites from aligned reads."""
    logger.info(f"Processing TSV file: {config.input_file}")
    logger.info(f"Chunk size: {config.chunk_size:,} rows")
    
    # Process all chunks
    insertion_counts, total_rows, total_valid_rows, total_invalid_rows, chunk_count = process_chunks(config)
    
    # Convert to output format
    logger.info("Preparing output table...")
    
    if not insertion_counts:
        return write_empty_output(config, total_rows, total_valid_rows, total_invalid_rows)
    
    # Create output DataFrame
    output_df, total_plus, total_minus = create_output_dataframe(insertion_counts)
    
    # Write output
    logger.info(f"Writing {len(output_df):,} insertion sites to {config.output_file}")
    output_df.to_csv(config.output_file, sep='\t', index=False)
    
    stats = ExtractionStats(
        total_rows=total_rows,
        valid_rows=total_valid_rows,
        invalid_rows=total_invalid_rows,
        unique_sites=len(output_df),
        total_plus_insertions=total_plus,
        total_minus_insertions=total_minus
    )
    
    # Display simplified summary
    logger.success(f"Processing complete: {stats.unique_sites:,} unique sites, {stats.total_insertions:,} total insertions")
    # Log detailed statistics
    logger.info("Statistics summary:")
    logger.info(f"  Total rows processed: {stats.total_rows:,}")
    logger.info(f"  Valid rows: {stats.valid_rows:,}")
    logger.info(f"  Invalid rows: {stats.invalid_rows:,}")
    logger.info(f"  Unique insertion sites: {stats.unique_sites:,}")
    logger.info(f"  Plus strand insertions: {stats.total_plus_insertions:,}")
    logger.info(f"  Minus strand insertions: {stats.total_minus_insertions:,}")
    logger.info(f"  Total insertions: {stats.total_insertions:,}")
    
    return stats


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(
        description="Extract insertion sites from aligned read TSV files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input TSV file from BAM parsing")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV file with insertion counts")
    parser.add_argument("-c", "--chunk_size", type=int, default=500000, help="Number of rows to process per chunk")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


@logger.catch
def main():
    """Main entry point of the script."""
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    logger.info(f"Pandas version: {pd.__version__}")

    # Validate input and output paths using the Pydantic model
    try:
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output,
            chunk_size=args.chunk_size
        )

        logger.info(f"Starting processing of {config.input_file}")

        # Run the core analysis/logic
        count_insertion_sites(config)
        
        logger.success("Processing completed successfully!")
        return 0
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())