#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "loguru",
#     "pyarrow",
# ]
# ///

"""
Extract Transposon Insertion Sites
===================================

Identify and count transposon insertion sites from the columnar read
alignment table produced by BAM parsing. Each aligned read carries a genomic
chromosome, reference start/end coordinates, and a strand orientation; the
insertion coordinate is derived from the strand-specific position of the TTAA
target motif.

The core algorithm streams the input in row-group batches (bounded memory),
keeps rows with a valid ``+``/``-`` strand and complete coordinates, and
computes the insertion coordinate per row: for the ``+`` strand the site is
``R1_Ref_Start + 4`` (the position immediately after ``TTAA``), and for the
``-`` strand it is ``R1_Ref_End``. Rows are then grouped by
``(chromosome, coordinate, strand)`` and counted, aggregating per-strand
tallies across all chunks into a single table.

Input
-----
- A Parquet file of aligned reads with columns ``R1_Strand``, ``R1_Chrom``,
  ``R1_Ref_Start``, and ``R1_Ref_End``.

Output
------
- A tab-separated TSV of insertion sites with columns ``Chr``, ``Coordinate``,
  ``+``, ``-`` (per-strand insertion counts), sorted by ``Chr`` then
  ``Coordinate`` (``sep="\\t"``, ``index=False``).

Usage
-----
    python extract_insertion_sites.py -i aligned.tsv -o insertions.tsv -c 500000
    python extract_insertion_sites.py --input aligned.tsv --output insertions.tsv --verbose

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
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd
import pyarrow.parquet as pq

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from read_processing.insertions import (  # noqa: E402
    InsertionCounts,
    ExtractionStats,
    extract_insertion_sites,
    create_output_dataframe,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
MIN_CHUNK_SIZE = 10000
MAX_CHUNK_SIZE = 5000000

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Configuration for input/output paths and chunked-processing parameters."""
    input_file: Path
    output_file: Path
    chunk_size: int = 500000

    def __post_init__(self) -> None:
        """Validate input path and chunk size, then ensure the output directory exists."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        if not (MIN_CHUNK_SIZE <= self.chunk_size <= MAX_CHUNK_SIZE):
            raise ValueError(
                f"Chunk size must be between {MIN_CHUNK_SIZE:,} and {MAX_CHUNK_SIZE:,}: {self.chunk_size}"
            )
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


setup_logger()

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def process_chunks(config: InputOutputConfig) -> tuple[InsertionCounts, int, int, int, int]:
    """Process all chunks and aggregate results."""
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    total_rows = 0
    total_valid_rows = 0
    total_invalid_rows = 0
    chunk_count = 0

    logger.info("Starting chunked processing...")

    parquet_file = pq.ParquetFile(config.input_file)

    for record_batch in parquet_file.iter_batches(batch_size=config.chunk_size):
        chunk_df = record_batch.to_pandas()
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

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Extract insertion sites from aligned read TSV files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input TSV file from BAM parsing")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV file with insertion counts")
    parser.add_argument("-c", "--chunk_size", type=int, default=500000, help="Number of rows to process per chunk")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate config, extract insertion sites, and write output."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    logger.info(f"Pandas version: {pd.__version__}")

    try:
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output,
            chunk_size=args.chunk_size,
        )

        logger.info(f"Starting processing of {config.input_file}")

        # Run the core analysis/logic
        count_insertion_sites(config)

        logger.success("Processing completed successfully!")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
