"""
Extract filtering statistics from read pair filtering log files.

This script processes log files generated during read pair filtering to extract
comprehensive statistics for PBL (left) and PBR (right) read pairs. It calculates
retention rates, processes counts, and generates summary statistics across samples.

The script parses log files containing filtering summaries with the following format:

.. code-block:: text

    ============================================================
    FILTERING SUMMARY
    ============================================================
    Total chunks processed: 12345
    Original read pairs: 1,234,567
    Filtered read pairs: 1,123,456
    Removed read pairs: 111,111
    Overall retention rate: 91.05%
    Output written to: sample_name.PBL.filtered.tsv

Typical Usage:
    python extract_mapping_filtering_statistics.py \
        -i sample1.log sample2.log \
        -o filtering_statistics.tsv \
        -v

Input:
    Log files from read pair filtering process containing filtering summaries.

Output:
    TSV file with filtering statistics including:
    - Original and filtered read pair counts
    - Retention rates per sample and read type
    - Total statistics across PBL and PBR pairs
"""

# =============================== Imports ===============================
import sys
import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional
from loguru import logger
from pydantic import BaseModel, Field, field_validator
import pandas as pd


# =============================== Constants ===============================
SUMMARY_PATTERN = re.compile(
    r".*\| ============================================================\s*\n"
    r".*\| FILTERING SUMMARY\s*\n"
    r".*\| ============================================================\s*\n"
    r".*\| Total chunks processed: (\d+)\s*\n"
    r".*\| Original read pairs: ([\d,]+)\s*\n"
    r".*\| Filtered read pairs: ([\d,]+)\s*\n"
    r".*\| Removed read pairs: ([\d,]+)\s*\n"
    r".*\| Overall retention rate: ([\d.]+)%\s*\n"
    r".*\| Output written to: (.+?\.(?:PBL|PBR)\.filtered\.tsv)"
)


# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_files: List[Path] = Field(..., description="List of log file paths")
    output_file: Path = Field(..., description="Path to save the output statistics TSV")

    @field_validator('input_files')
    def validate_input_files(cls, v):
        """Validate that all input files exist and are readable."""
        for file_path in v:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        return v

    @field_validator('output_file')
    def validate_output_file(cls, v):
        """Validate output file path and create parent directories."""
        v.parent.mkdir(parents=True, exist_ok=True)
        return v

    class Config:
        frozen = True


class FilteringStatistics(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    chunks_processed_pbl: Optional[int] = None
    original_read_pairs_pbl: Optional[int] = None
    filtered_read_pairs_pbl: Optional[int] = None
    removed_read_pairs_pbl: Optional[int] = None
    retention_rate_pbl: Optional[float] = None
    output_file_pbl: Optional[str] = None
    chunks_processed_pbr: Optional[int] = None
    original_read_pairs_pbr: Optional[int] = None
    filtered_read_pairs_pbr: Optional[int] = None
    removed_read_pairs_pbr: Optional[int] = None
    retention_rate_pbr: Optional[float] = None
    output_file_pbr: Optional[str] = None
    total_original_pairs: Optional[int] = None
    total_filtered_pairs: Optional[int] = None
    overall_retention_rate: Optional[float] = None


class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    total_samples_processed: int = Field(..., ge=0, description="Total number of samples processed")
    total_log_files: int = Field(..., ge=0, description="Total number of log files processed")
    output_path: Path = Field(..., description="Path to the output file")


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
def parse_log_file(log_file: Path) -> Dict[str, FilteringStatistics]:
    """Parse a single log file and extract filtering statistics."""
    sample_name = log_file.stem
    logger.info(f"Processing: {sample_name}")
    
    try:
        with open(log_file, "r") as f:
            content = f.read()
    except Exception as e:
        logger.error(f"Error reading {log_file}: {str(e)}")
        return {}
    
    matches = SUMMARY_PATTERN.findall(content)
    
    if not matches:
        logger.warning(f"No filtering summary sections found in: {sample_name}")
        return {}
    
    logger.debug(f"Found {len(matches)} filtering summary sections in {sample_name}")
    
    stats_dict = {}
    
    for match in matches:
        chunks_processed = int(match[0])
        original_pairs = int(match[1].replace(',', ''))
        filtered_pairs = int(match[2].replace(',', ''))
        removed_pairs = int(match[3].replace(',', ''))
        retention_rate = float(match[4]) / 100
        output_path = match[5]
        
        # Determine if this is PBL or PBR based on output path
        if ".PBL.filtered.tsv" in output_path:
            suffix = "pbl"
        elif ".PBR.filtered.tsv" in output_path:
            suffix = "pbr"
        else:
            logger.warning(f"Could not determine PBL/PBR from output path: {output_path}")
            continue
        
        # Update statistics
        stats_dict.update({
            f"chunks_processed_{suffix}": chunks_processed,
            f"original_read_pairs_{suffix}": original_pairs,
            f"filtered_read_pairs_{suffix}": filtered_pairs,
            f"removed_read_pairs_{suffix}": removed_pairs,
            f"retention_rate_{suffix}": retention_rate,
            f"output_file_{suffix}": output_path
        })
        
        logger.debug(f"  {suffix.upper()}: {original_pairs:,} -> {filtered_pairs:,} ({retention_rate*100:.2f}% retained)")
    
    return {sample_name: FilteringStatistics(**stats_dict)}


@logger.catch
def extract_summary_data(log_files: List[Path]) -> Dict[str, FilteringStatistics]:
    """Extract filtering statistics from multiple log files."""
    logger.info(f"Found {len(log_files)} log files with filtering statistics")
    
    all_statistics = {}
    
    for log_file in log_files:
        file_stats = parse_log_file(log_file)
        all_statistics.update(file_stats)
    
    return all_statistics


@logger.catch
def create_dataframe(statistics: Dict[str, FilteringStatistics]) -> pd.DataFrame:
    """Create a pandas DataFrame from filtering statistics dictionary."""
    if not statistics:
        logger.error("No statistics extracted from any log files")
        return pd.DataFrame()
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(
        {sample: stats.model_dump(exclude_none=True) for sample, stats in statistics.items()},
        orient='index'
    )
    
    # Sort columns for better readability
    pbl_cols = [col for col in df.columns if col.endswith('_pbl')]
    pbr_cols = [col for col in df.columns if col.endswith('_pbr')]
    all_cols = sorted(pbl_cols) + sorted(pbr_cols)
    
    # Ensure all expected columns exist
    for col in all_cols:
        if col not in df.columns:
            df[col] = None
    
    df = df.reindex(columns=all_cols)
    
    # Calculate totals
    if 'original_read_pairs_pbl' in df.columns and 'original_read_pairs_pbr' in df.columns:
        df['total_original_pairs'] = df[['original_read_pairs_pbl', 'original_read_pairs_pbr']].sum(axis=1)
    
    if 'filtered_read_pairs_pbl' in df.columns and 'filtered_read_pairs_pbr' in df.columns:
        df['total_filtered_pairs'] = df[['filtered_read_pairs_pbl', 'filtered_read_pairs_pbr']].sum(axis=1)
    
    # Calculate overall retention rate
    if 'total_original_pairs' in df.columns and 'total_filtered_pairs' in df.columns:
        df['overall_retention_rate'] = (
            df['total_filtered_pairs'] / df['total_original_pairs']
        ).round(4)
    
    return df


@logger.catch
def save_results(df: pd.DataFrame, output_file: Path) -> None:
    """Save the filtering statistics DataFrame to a TSV file."""
    if df.empty:
        logger.error("Cannot save empty DataFrame")
        return
    
    # Sort by sample name
    df = df.rename_axis("Sample", axis=0).sort_index()
    
    # Save to file
    df.to_csv(output_file, sep="\t", index=True, float_format="%.2f")
    logger.success(f"Statistics saved to: {output_file}")
    
    # Display summary
    summary_cols = ['total_original_pairs', 'total_filtered_pairs', 'overall_retention_rate']
    available_cols = [col for col in summary_cols if col in df.columns]
    
    if available_cols:
        logger.info("Summary statistics:")
        logger.info(f"\n{df[available_cols].describe()}")


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(
        description="Extract filtering statistics from read pair filtering log files"
    )
    parser.add_argument(
        "-i", "--input", 
        type=Path, 
        nargs="+", 
        required=True, 
        help="Path to the log files"
    )
    parser.add_argument(
        "-o", "--output", 
        type=Path, 
        required=True, 
        help="Path to save the output statistics TSV"
    )
    parser.add_argument(
        "-v", "--verbose", 
        action="store_true", 
        help="Enable verbose logging"
    )
    return parser.parse_args()


@logger.catch
def main():
    """Main entry point of the script. Replace this docstring with a relevant one."""
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)
    
    try:
        # Validate inputs
        config = InputOutputConfig(
            input_files=args.input,
            output_file=args.output
        )
        
        logger.info("Starting filtering statistics extraction")
        
        # Extract statistics
        statistics = extract_summary_data(config.input_files)
        
        if not statistics:
            logger.error("No statistics extracted from any log files")
            sys.exit(1)
        
        # Create DataFrame
        df = create_dataframe(statistics)
        
        # Save results
        save_results(df, config.output_file)
        
        # Create analysis result
        result = AnalysisResult(
            total_samples_processed=len(statistics),
            total_log_files=len(config.input_files),
            output_path=config.output_file
        )
        
        logger.success("Extraction completed successfully!")
        logger.info(f"Processed {result.total_samples_processed} samples from {result.total_log_files} log files")
        logger.info(f"Output saved to: {result.output_path}")
        
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()