#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "loguru",
#     "pyyaml",
#     "pyarrow",
# ]
# ///

"""
Filter Aligned Read Pairs Using YAML Configuration
==================================================

Filter aligned read pairs from BAM-derived Parquet files using thresholds
loaded from a YAML configuration file. R1 and R2 reads are filtered
independently with separate MAPQ, NCIGAR, and NM thresholds, plus optional
rejection of reads that carry supplementary (SA) or secondary (XA) alignments
and an optional proper-pair requirement.

All filtering criteria are read from the ``aligned_read_filtering`` section of
the config so runs stay consistent and reproducible. The input Parquet file is
streamed in row-group batches and each batch is filtered then written to the
output as its own row group, keeping memory use bounded for large files.

Input
-----
- Parquet file with read-pair data from BAM parsing. Expected columns include
  ``R1_MAPQ``, ``R2_MAPQ``, ``R1_NCIGAR``, ``R2_NCIGAR``, ``R1_NM``, ``R2_NM``,
  ``R1_SA``, ``R2_SA``, ``R1_XA``, ``R2_XA``, and ``Is_Proper_Pair``.
- YAML config with an ``aligned_read_filtering`` section containing
  ``read_1_filtering``, ``read_2_filtering``, and ``require_proper_pair``.

Output
------
- Filtered Parquet file containing only the read pairs that pass every
  configured criterion.

Usage
-----
    python filter_aligned_reads.py -i input.parquet -o filtered.parquet --config config.yaml
    python filter_aligned_reads.py -i input.parquet -o filtered.parquet --config config.yaml -c 3000000

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""
# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# 2. Data Processing Imports
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from read_processing.filtering import (  # noqa: E402
    FilterThresholds,
    AnalysisResult,
    load_config_from_yaml,
    strip_read_prefix,
    coerce_column_dtypes,
    build_filter_mask,
    process_chunk,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
DEFAULT_CHUNK_SIZE = 50000

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Complete filtering configuration: paths, chunk size, and loaded thresholds."""
    input_file: Path
    output_file: Path
    chunk_size: int
    config_data: dict[str, Any]

    def __post_init__(self) -> None:
        """Validate the input path and ensure the output directory exists."""
        if not self.input_file.exists():
            raise ValueError(f"File not found: {self.input_file}")
        output_dir = self.output_file.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Statistics from a filtering run."""
    total_rows: int
    filtered_rows: int
    removed_rows: int
    retention_rate: float
    chunks_processed: int


setup_logger()
# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def filter_read_pairs(config: InputOutputConfig) -> AnalysisResult:
    """Stream the input Parquet file in batches, filter each, write to output, and return run statistics."""
    logger.info(f"Loading data from: {config.input_file}")

    # Initialize counters
    total_rows = 0
    filtered_rows = 0
    chunk_count = 0
    first_chunk = True
    writer: pq.ParquetWriter | None = None
    output_schema: pa.Schema | None = None

    logger.info(f"Processing file in chunks of {config.chunk_size:,} rows...")

    parquet_file = pq.ParquetFile(config.input_file)

    try:
        # Process each chunk
        for record_batch in parquet_file.iter_batches(batch_size=config.chunk_size):
            chunk_df = record_batch.to_pandas()
            chunk_count += 1
            total_rows += len(chunk_df)

            if chunk_count % 10 == 0:
                logger.info(f"Processing chunk {chunk_count}, total rows: {total_rows:,}")

            # Process chunk
            filtered_chunk, first_chunk = process_chunk(
                chunk_df,
                chunk_count,
                config.config_data["read_1_filtering"],
                config.config_data["read_2_filtering"],
                config.config_data["require_proper_pair"],
                first_chunk,
            )
            filtered_rows += len(filtered_chunk)

            # Write filtered chunk as its own row group
            if writer is None:
                output_schema = pa.Schema.from_pandas(filtered_chunk, preserve_index=False)
                writer = pq.ParquetWriter(config.output_file, output_schema)
                logger.info(f"Created output file: {config.output_file}")

            table = pa.Table.from_pandas(filtered_chunk, schema=output_schema, preserve_index=False)
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()

    logger.info(f"Completed processing {chunk_count} chunks")

    # Calculate statistics
    removed_rows = total_rows - filtered_rows
    retention_rate = filtered_rows / total_rows * 100 if total_rows > 0 else 0

    stats = AnalysisResult(
        total_rows=total_rows,
        filtered_rows=filtered_rows,
        removed_rows=removed_rows,
        retention_rate=retention_rate,
        chunks_processed=chunk_count,
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
            sample_df = pq.ParquetFile(config.output_file).read_row_group(0).to_pandas().head(5)
            logger.debug("Sample of filtered data (first 5 rows):")
            logger.debug(f"Shape: {sample_df.shape}")
        except Exception as e:
            logger.warning(f"Could not read sample of filtered data: {e}")

    return stats

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Filter aligned read pairs using configuration from YAML file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input Parquet file with read pair data",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output Parquet file for filtered data",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="YAML configuration file with filtering parameters",
    )

    # Chunking configuration
    parser.add_argument(
        "-c", "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Number of rows to process per chunk",
    )

    # Logging configuration
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable DEBUG level logging",
    )

    return parser.parse_args()


def main() -> int:
    """Main orchestrator: load YAML config, then filter read pairs with chunked processing."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        # Load configuration from YAML file
        config_data = load_config_from_yaml(Path(args.config))

        config = InputOutputConfig(
            input_file=Path(args.input),
            output_file=Path(args.output),
            chunk_size=args.chunk_size,
            config_data=config_data,
        )

        # Process the file
        filter_read_pairs(config)
        logger.success("Filtering completed successfully")

    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())



