"""
Filter aligned read pairs using YAML configuration.

This script filters aligned read pairs from BAM-derived TSV files using configuration
parameters loaded from a YAML file. It supports independent filtering for R1 and R2
reads with configurable quality thresholds, alignment criteria, and chunked processing
for memory efficiency. All filtering parameters are loaded from the YAML configuration
file to ensure consistent and reproducible filtering across analyses.

Key Features:
- YAML-based configuration for all filtering parameters
- Independent filtering for R1 and R2 reads with separate thresholds
- Support for MAPQ, NCIGAR, and NM value filtering
- Filtering based on supplementary (SA) and secondary (XA) alignments
- Proper pair validation options
- Memory-efficient chunked processing for large files
- Comprehensive validation using Pydantic models
- Structured logging with Loguru

Typical Usage:
    python filter_aligned_reads.py --input-file input.tsv --output-file filtered.tsv --config-file config.yaml
    python filter_aligned_reads.py -i input.tsv -o filtered.tsv --config-file config.yaml -c 100000

Input: TSV file with read pair data from BAM parsing (columns include R1_MAPQ, R2_MAPQ, R1_NCIGAR, R2_NCIGAR, R1_NM, R2_NM, R1_SA, R2_SA, R1_XA, R2_XA, Is_Proper_Pair)
Output: Filtered TSV file with high-quality read pairs that pass all specified criteria
Config: YAML file with aligned_read_filtering section containing read_1_filtering and read_2_filtering parameters
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
import pandas as pd
import yaml
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# =============================== Constants ===============================
NA_VALUES = ['N/A', 'NA', '']
DEFAULT_CHUNK_SIZE = 50000


# =============================== Configuration & Models ===============================

class FilterThresholds(BaseModel):
    """Validation model for read filtering thresholds."""
    
    mapq_threshold: Optional[float] = Field(None, ge=0, le=255, description="Minimum MAPQ score")
    ncigar_value: Optional[int] = Field(None, ge=0, description="Required NCIGAR value")
    nm_threshold: Optional[int] = Field(None, ge=0, description="Maximum mismatches allowed")
    no_sa: bool = Field(False, description="Require no supplementary alignments")
    no_xa: bool = Field(False, description="Require no secondary alignments")
    
    class Config:
        frozen = True


class InputOutputConfig(BaseModel):
    """Complete filtering configuration with validation."""
    
    input_file: Path = Field(..., description="Input TSV file path")
    output_file: Path = Field(..., description="Output TSV file path")
    chunk_size: int = Field(DEFAULT_CHUNK_SIZE, ge=1000, le=10000000, description="Rows per chunk")
    config_data: Dict[str, Any] = Field(..., description="Configuration data")
    
    @field_validator('input_file')
    def validate_input_exists(cls, v):
        if not v.exists():
            raise ValueError(f"File not found: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_dir(cls, v):
        output_dir = v.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True


class AnalysisResult(BaseModel):
    """Statistics from filtering operation."""
    
    total_rows: int = Field(..., ge=0)
    filtered_rows: int = Field(..., ge=0)
    removed_rows: int = Field(..., ge=0)
    retention_rate: float = Field(..., ge=0, le=100)
    chunks_processed: int = Field(..., ge=0)
    
    class Config:
        frozen = True


# =============================== Setup Logging ===============================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for read filtering."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# =============================== Core Functions ===============================
@logger.catch
def load_config_from_yaml(config_file: Path) -> Dict[str, Any]:
    """Load filtering configuration from YAML file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Extract aligned_read_filtering configuration
        if 'aligned_read_filtering' not in config:
            raise ValueError("'aligned_read_filtering' section not found in config file")
        
        filtering_config = config['aligned_read_filtering']
        
        # Convert YAML config to internal format
        internal_config = {}
        for read in ['read_1_filtering', 'read_2_filtering']:
            internal_config[read] = FilterThresholds(
                mapq_threshold=filtering_config.get(read, {}).get("mapq_threshold"),
                ncigar_value=filtering_config.get(read, {}).get("ncigar_value"),
                nm_threshold=filtering_config.get(read, {}).get("nm_threshold"),
                no_sa=filtering_config.get(read, {}).get("no_sa"),
                no_xa=filtering_config.get(read, {}).get("no_xa")
            )
        
        internal_config['require_proper_pair'] = filtering_config.get("require_proper_pair")
        
        logger.info(f"Loaded configuration from: {config_file}")
        logger.debug(f"R1 filters: MAPQ={internal_config['read_1_filtering'].mapq_threshold}, NCIGAR={internal_config['read_1_filtering'].ncigar_value}, NM={internal_config['read_1_filtering'].nm_threshold}")
        logger.debug(f"R2 filters: MAPQ={internal_config['read_2_filtering'].mapq_threshold}, NCIGAR={internal_config['read_2_filtering'].ncigar_value}, NM={internal_config['read_2_filtering'].nm_threshold}")
        logger.debug(f"Pair filters: no_sa={internal_config['read_1_filtering'].no_sa}, no_xa={internal_config['read_1_filtering'].no_xa}, proper_pair={internal_config['require_proper_pair']}")
        
        return internal_config
        
    except Exception as e:
        logger.error(f"Error loading config file {config_file}: {e}")
        raise

@logger.catch
def build_filter_mask(
    chunk: pd.DataFrame,
    r1_filters: FilterThresholds,
    r2_filters: FilterThresholds,
    require_proper_pair: bool
) -> pd.Series:
    """Build boolean mask for filtering read pairs."""
    # Initialize mask with all True
    filter_mask = pd.Series([True] * len(chunk), index=chunk.index)
    
    # Apply R1 filters
    if r1_filters.mapq_threshold is not None:
        filter_mask &= (chunk['R1_MAPQ'] >= r1_filters.mapq_threshold)
        
    if r1_filters.ncigar_value is not None:
        filter_mask &= (chunk['R1_NCIGAR'] <= r1_filters.ncigar_value)
        
    if r1_filters.nm_threshold is not None:
        filter_mask &= (chunk['R1_NM'] <= r1_filters.nm_threshold)
        
    if r1_filters.no_sa:
        filter_mask &= (chunk['R1_SA'].isna() | (chunk['R1_SA'] == 'N/A'))
        
    if r1_filters.no_xa:
        filter_mask &= (chunk['R1_XA'].isna() | (chunk['R1_XA'] == 'N/A'))
    
    # Apply R2 filters
    if r2_filters.mapq_threshold is not None:
        filter_mask &= (chunk['R2_MAPQ'] >= r2_filters.mapq_threshold)
        
    if r2_filters.ncigar_value is not None:
        filter_mask &= (chunk['R2_NCIGAR'] <= r2_filters.ncigar_value)
        
    if r2_filters.nm_threshold is not None:
        filter_mask &= (chunk['R2_NM'] <= r2_filters.nm_threshold)
        
    if r2_filters.no_sa:
        filter_mask &= (chunk['R2_SA'].isna() | (chunk['R2_SA'] == 'N/A'))
        
    if r2_filters.no_xa:
        filter_mask &= (chunk['R2_XA'].isna() | (chunk['R2_XA'] == 'N/A'))
    
    # Apply proper pair filter
    if require_proper_pair:
        filter_mask &= (chunk['Is_Proper_Pair'].str.capitalize() == 'Yes')
    
    return filter_mask




@logger.catch
def process_chunk(
    chunk: pd.DataFrame,
    chunk_num: int,
    config: Dict[str, Any],
    first_chunk: bool
) -> Tuple[pd.DataFrame, bool]:
    """Process a single chunk of data."""
    chunk_rows_before = len(chunk)
    
    # Display info for first chunk
    if first_chunk:
        logger.info("=" * 60)
        logger.info("ORIGINAL DATA INFORMATION")
        logger.info("=" * 60)
        logger.info(f"Columns: {len(chunk.columns)}")
        logger.info(f"First chunk size: {chunk_rows_before:,} rows")
        
        logger.debug("Column Data Types:")
        for col, dtype in chunk.dtypes.items():
            logger.debug(f"  {col}: {dtype}")
    
    # Build and apply filter mask
    filter_mask = build_filter_mask(
        chunk,
        config['read_1_filtering'],
        config['read_2_filtering'],
        config['require_proper_pair']
    )
    
    filtered_chunk = chunk[filter_mask]
    chunk_filtered_rows = len(filtered_chunk)
    
    # Log progress
    if chunk_num == 1 or chunk_num % 10 == 0:
        retention_rate = (chunk_filtered_rows / chunk_rows_before * 100 
                         if chunk_rows_before > 0 else 0)
        logger.info(
            f"Chunk {chunk_num}: {chunk_filtered_rows:,}/{chunk_rows_before:,} "
            f"rows retained ({retention_rate:.1f}%)"
        )
    
    return filtered_chunk, False


@logger.catch
def filter_read_pairs(config: InputOutputConfig) -> AnalysisResult:
    """Filter read pairs using chunked processing."""
    logger.info(f"Loading data from: {config.input_file}")
    
    # Initialize counters
    total_rows = 0
    filtered_rows = 0
    chunk_count = 0
    first_chunk = True
    
    logger.info(f"Processing file in chunks of {config.chunk_size:,} rows...")
    
    try:
        # Create chunk iterator
        chunk_iterator = pd.read_csv(
            config.input_file,
            sep='\t',
            na_values=NA_VALUES,
            chunksize=config.chunk_size
        )
        
        # Process each chunk
        for chunk_df in chunk_iterator:
            chunk_count += 1
            total_rows += len(chunk_df)
            
            if chunk_count % 10 == 0:
                logger.info(f"Processing chunk {chunk_count}, total rows: {total_rows:,}")
            
            # Process chunk
            filtered_chunk, first_chunk = process_chunk(
                chunk_df, chunk_count, config.config_data, first_chunk
            )
            filtered_rows += len(filtered_chunk)
            
            # Write filtered chunk
            if chunk_count == 1:
                filtered_chunk.to_csv(
                    config.output_file, sep='\t', index=False, mode='w'
                )
                logger.info(f"Created output file: {config.output_file}")
            else:
                filtered_chunk.to_csv(
                    config.output_file, sep='\t', index=False, 
                    mode='a', header=False
                )
        
        logger.info(f"Completed processing {chunk_count} chunks")
        
        # Calculate statistics
        removed_rows = total_rows - filtered_rows
        retention_rate = filtered_rows / total_rows * 100 if total_rows > 0 else 0
        
        stats = AnalysisResult(
            total_rows=total_rows,
            filtered_rows=filtered_rows,
            removed_rows=removed_rows,
            retention_rate=retention_rate,
            chunks_processed=chunk_count
        )
        
        # Display summary
        logger.info("=" * 60)
        logger.info("FILTERING SUMMARY")
        logger.info("=" * 60)
        logger.success(f"Total chunks processed: {stats.chunks_processed}")
        logger.info(f"Original read pairs: {stats.total_rows:,}")
        logger.info(f"Filtered read pairs: {stats.filtered_rows:,}")
        logger.info(f"Removed read pairs: {stats.removed_rows:,}")
        logger.success(f"Overall retention rate: {stats.retention_rate:.2f}%")
        logger.info(f"Output written to: {config.output_file}")
        
        # Display sample of filtered data
        if filtered_rows > 0:
            try:
                sample_df = pd.read_csv(config.output_file, sep='\t', nrows=5)
                logger.debug("Sample of filtered data (first 5 rows):")
                logger.debug(f"Shape: {sample_df.shape}")
            except Exception as e:
                logger.warning(f"Could not read sample of filtered data: {e}")
        
        return stats
        
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        raise


# =============================== Main Function ===============================

def parse_arguments():
    """Set command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter aligned read pairs using configuration from YAML file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Input TSV file with read pair data"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Output TSV file for filtered data"
    )
    parser.add_argument(
        "--config", 
        required=True, 
        help="YAML configuration file with filtering parameters"
    )
    
    # Chunking configuration
    parser.add_argument(
        "-c", "--chunk-size", 
        type=int, 
        default=50000, 
        help="Number of rows to process per chunk"
    )
    
    return parser.parse_args()

@logger.catch
def main():
    """Main execution function for read filtering."""
    args = parse_arguments()
    setup_logging()
    
    try:
        # Load configuration from YAML file
        config_data = load_config_from_yaml(Path(args.config))
        
        config = InputOutputConfig(
            input_file=Path(args.input),
            output_file=Path(args.output),
            chunk_size=args.chunk_size,
            config_data=config_data
        )
        
        # Process the file
        filter_read_pairs(config)
        logger.success("Filtering completed successfully")
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()