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
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# 3. Third-party Imports
from loguru import logger
import yaml

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
NA_VALUES = ["N/A", "NA", ""]
DEFAULT_CHUNK_SIZE = 50000

# SAM optional-field type codes: these tags are numeric ('i'/'f') per the SAM
# spec (MC/MD/SA/XA are strings and need no dtype coercion below).
NUMERIC_TAGS = frozenset({"AS", "MQ", "NM", "XS"})

# Base (non-tag) fields that are always populated by parse_bam_to_tsv.py and
# never carry "N/A" (mapq/n_cigar/flag default to 0, never "N/A" -- see
# ReadInfo).
ALWAYS_INT_FIELDS = frozenset({"MAPQ", "NCIGAR", "Flag", "FLAG"})

# Base fields plus numeric tags that CAN be "N/A" -- coerced to nullable float
# to match the dtype pd.read_csv would infer for a column containing NaN.
NULLABLE_NUMERIC_FIELDS = frozenset({"LEN", "Pos", "Ref_Start", "Ref_End"}) | NUMERIC_TAGS

# Columns checked via .isna() in build_filter_mask -- need the "N/A"/"NA"/""
# string sentinels turned into real NaN. Every other string column is left
# untouched: pd.to_numeric(errors="coerce") already treats those sentinels as
# unparseable on its own, and no other string column is null-checked.
ISNA_CHECKED_COLUMNS = ("R1_SA", "R1_XA", "R2_SA", "R2_XA")

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FilterThresholds:
    """Per-read filtering thresholds for MAPQ, NCIGAR, NM, and SA/XA alignments."""
    mapq_threshold: float | None = None
    ncigar_value: int | None = None
    nm_threshold: int | None = None
    no_sa: bool = False
    no_xa: bool = False


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

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for read filtering."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

setup_logger()
# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
def load_config_from_yaml(config_file: Path) -> dict[str, Any]:
    """Load filtering thresholds from the ``aligned_read_filtering`` section of a YAML file."""
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Extract aligned_read_filtering configuration
    if "aligned_read_filtering" not in config:
        raise ValueError("'aligned_read_filtering' section not found in config file")

    filtering_config = config["aligned_read_filtering"]

    # Convert YAML config to internal format
    internal_config: dict[str, Any] = {}
    for read in ["read_1_filtering", "read_2_filtering"]:
        internal_config[read] = FilterThresholds(
            mapq_threshold=filtering_config.get(read, {}).get("mapq_threshold"),
            ncigar_value=filtering_config.get(read, {}).get("ncigar_value"),
            nm_threshold=filtering_config.get(read, {}).get("nm_threshold"),
            no_sa=filtering_config.get(read, {}).get("no_sa"),
            no_xa=filtering_config.get(read, {}).get("no_xa"),
        )

    internal_config["require_proper_pair"] = filtering_config.get("require_proper_pair")

    logger.info(f"Loaded configuration from: {config_file}")
    logger.debug(f"R1 filters: MAPQ={internal_config['read_1_filtering'].mapq_threshold}, NCIGAR={internal_config['read_1_filtering'].ncigar_value}, NM={internal_config['read_1_filtering'].nm_threshold}")
    logger.debug(f"R2 filters: MAPQ={internal_config['read_2_filtering'].mapq_threshold}, NCIGAR={internal_config['read_2_filtering'].ncigar_value}, NM={internal_config['read_2_filtering'].nm_threshold}")
    logger.debug(f"Pair filters: no_sa={internal_config['read_1_filtering'].no_sa}, no_xa={internal_config['read_1_filtering'].no_xa}, proper_pair={internal_config['require_proper_pair']}")

    return internal_config


def strip_read_prefix(column: str) -> str:
    """Strip a leading R1_/R2_ prefix from a column name, if present."""
    if column.startswith("R1_") or column.startswith("R2_"):
        return column[3:]
    return column


def coerce_column_dtypes(chunk: pd.DataFrame) -> pd.DataFrame:
    """Coerce an all-string Parquet chunk to the dtypes pd.read_csv would infer from the equivalent TSV."""
    for column in ISNA_CHECKED_COLUMNS:
        if column in chunk.columns:
            chunk[column] = chunk[column].replace(NA_VALUES, np.nan)
    for column in chunk.columns:
        field = strip_read_prefix(column)
        if field in ALWAYS_INT_FIELDS:
            chunk[column] = chunk[column].astype("int64")
        elif field in NULLABLE_NUMERIC_FIELDS:
            chunk[column] = pd.to_numeric(chunk[column], errors="coerce").astype("float64")
    return chunk


@logger.catch
def build_filter_mask(
    chunk: pd.DataFrame,
    r1_filters: FilterThresholds,
    r2_filters: FilterThresholds,
    require_proper_pair: bool,
) -> pd.Series:
    """Build a boolean mask selecting read pairs that pass every configured criterion."""
    # Initialize mask with all True
    filter_mask = pd.Series([True] * len(chunk), index=chunk.index)

    # Apply R1 filters
    if r1_filters.mapq_threshold is not None:
        filter_mask &= (chunk["R1_MAPQ"] >= r1_filters.mapq_threshold)

    if r1_filters.ncigar_value is not None:
        filter_mask &= (chunk["R1_NCIGAR"] <= r1_filters.ncigar_value)

    if r1_filters.nm_threshold is not None:
        filter_mask &= (chunk["R1_NM"] <= r1_filters.nm_threshold)

    if r1_filters.no_sa:
        filter_mask &= (chunk["R1_SA"].isna() | (chunk["R1_SA"] == "N/A"))

    if r1_filters.no_xa:
        filter_mask &= (chunk["R1_XA"].isna() | (chunk["R1_XA"] == "N/A"))

    # Apply R2 filters
    if r2_filters.mapq_threshold is not None:
        filter_mask &= (chunk["R2_MAPQ"] >= r2_filters.mapq_threshold)

    if r2_filters.ncigar_value is not None:
        filter_mask &= (chunk["R2_NCIGAR"] <= r2_filters.ncigar_value)

    if r2_filters.nm_threshold is not None:
        filter_mask &= (chunk["R2_NM"] <= r2_filters.nm_threshold)

    if r2_filters.no_sa:
        filter_mask &= (chunk["R2_SA"].isna() | (chunk["R2_SA"] == "N/A"))

    if r2_filters.no_xa:
        filter_mask &= (chunk["R2_XA"].isna() | (chunk["R2_XA"] == "N/A"))

    # Apply proper pair filter
    if require_proper_pair:
        filter_mask &= (chunk["Is_Proper_Pair"].str.capitalize() == "Yes")

    return filter_mask
@logger.catch
def process_chunk(
    chunk: pd.DataFrame,
    chunk_num: int,
    config: dict[str, Any],
    first_chunk: bool,
) -> tuple[pd.DataFrame, bool]:
    """Filter a single chunk and log first-chunk diagnostics and periodic progress."""
    chunk = coerce_column_dtypes(chunk)
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
        config["read_1_filtering"],
        config["read_2_filtering"],
        config["require_proper_pair"],
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
                chunk_df, chunk_count, config.config_data, first_chunk
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



