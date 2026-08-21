#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
#     "pysam",
#     "pyarrow",
# ]
# ///

"""
BAM to Parquet Read-Pair Parser
================================

Extract per-read-pair alignment summaries from a QNAME-sorted BAM/SAM file and
write them as a columnar table. For each query template the parser records
read 1 and read 2 mapping quality, alignment length, CIGAR string, strand,
number of CIGAR operations, reference name, position, reference start/end,
SAM flag, a configurable set of alignment tags, and the proper-pair status.

The core algorithm is a single streaming pass over the alignment file: it relies
on the input being QNAME-sorted so that all alignments sharing a query name are
adjacent. Unmapped, secondary, and supplementary alignments are skipped; the
first primary read 1 and read 2 seen for a query name are retained. When a new
query name is encountered the previous pair is formatted and buffered, keeping
memory usage bounded regardless of file size (rows are flushed to disk every
``PARQUET_BATCH_SIZE`` pairs). Multi-threaded BAM decompression is supported
via pysam.

Input
-----
- A QNAME-sorted BAM (``.bam``, read as ``rb``) or SAM (any other suffix, read
  as ``r``) file.

Output
------
- A Parquet file with one row per read pair. The 14 base numeric fields
  (``MAPQ``, ``LEN``, ``NCIGAR``, ``Pos``, ``Ref_Start``, ``Ref_End``, ``Flag``
  for both reads) are stored as ``int64`` (nullable: a missing value is a real
  null, not the ``"N/A"`` sentinel). Every other column -- ``QueryName``, CIGAR,
  strand, chrom, the SAM tag columns, and ``Is_Proper_Pair`` -- is stored as a
  string, since tags can be non-numeric and are consumed as strings downstream.
  Typing the base numeric columns natively (rather than as strings) is markedly
  cheaper to encode into Parquet and is verified to leave the downstream
  ``filter_aligned_reads`` output byte- and dtype-identical. Tags are emitted in
  sorted order.

Usage
-----
    python parse_bam_to_tsv.py -i input.bam -o output.parquet -t 8
    python parse_bam_to_tsv.py --input input.bam --output output.parquet --threads 8 --verbose

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
from dataclasses import dataclass, field
from pathlib import Path

# 2. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
from read_processing.bam import (  # noqa: E402
    ReadInfo,
    ReadPairInfo,
    extract_read_info,
    determine_proper_pair_status,
    process_read_pair,
    format_output_line,
    build_header,
    build_schema,
    flush_batch,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# Default tags to extract from BAM/SAM files
DEFAULT_TAGS = ["AS", "MC", "MD", "MQ", "NM", "SA", "XA", "XS"]

# Progress reporting intervals
READ_PROGRESS_INTERVAL = 2000000
PAIR_PROGRESS_INTERVAL = 500000

# Read pairs buffered in memory before each Parquet row-group flush
PARQUET_BATCH_SIZE = 500000

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Input/output paths and processing parameters for BAM to Parquet conversion."""
    input_bam: Path
    output_file: Path
    threads: int = 4
    tag_list: list[str] = field(default_factory=lambda: list(DEFAULT_TAGS))

    def __post_init__(self) -> None:
        """Validate the input file and thread count, then ensure the output directory exists."""
        if not self.input_bam.exists():
            raise ValueError(f"Input file not found: {self.input_bam}")
        if not 1 <= self.threads <= 32:
            raise ValueError(f"threads must be between 1 and 32, got {self.threads}")
        output_dir = self.output_file.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)


setup_logger()

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def process_bam_file(config: InputOutputConfig) -> None:
    """Process a QNAME-sorted BAM/SAM file and write the read-pair Parquet output."""
    logger.info(f"Starting BAM processing with {config.threads} threads")
    logger.info(f"Input: {config.input_bam}")
    logger.info(f"Output: {config.output_file}")
    logger.info(f"Extracting tags: {', '.join(config.tag_list)}")

    # Ensure tags are sorted for consistent output
    sorted_tags = sorted(config.tag_list)

    # Build header and schema
    header_fields = build_header(sorted_tags)
    schema = build_schema(header_fields)

    # Initialize tracking variables
    current_qname = None
    current_r1 = None
    current_r2 = None
    processed_qname_count = 0
    read_count = 0
    buffered_rows: list[list] = []

    # Open BAM/SAM file
    mode = "rb" if str(config.input_bam).endswith(".bam") else "r"
    samfile = pysam.AlignmentFile(
        str(config.input_bam), mode, threads=config.threads
    )

    logger.info(
        "Processing reads in streaming mode (requires qname-sorted input)"
    )

    # Statistics/dictionary encoding buy nothing here: the output is read once,
    # sequentially, in full by filter_aligned_reads.py (no predicate pushdown),
    # and is temp()-deleted immediately after. Disabling both cuts writer
    # overhead with zero effect on the decoded values.
    with pq.ParquetWriter(
        config.output_file, schema, write_statistics=False, use_dictionary=False
    ) as writer:
        for read in samfile:
            read_count += 1

            if read_count % READ_PROGRESS_INTERVAL == 0:
                logger.info(f"Processed {read_count // 1000000}M alignments")

            # Skip unmapped, secondary, and supplementary alignments
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name

            # Process completed pair when encountering new qname
            if qname != current_qname:
                if current_qname is not None:
                    pair_info = process_read_pair(
                        current_qname, current_r1, current_r2, sorted_tags
                    )
                    buffered_rows.append(format_output_line(pair_info, sorted_tags))

                    processed_qname_count += 1
                    if processed_qname_count % PAIR_PROGRESS_INTERVAL == 0:
                        logger.info(f"Written {processed_qname_count} read pairs")
                    if len(buffered_rows) >= PARQUET_BATCH_SIZE:
                        flush_batch(writer, buffered_rows, schema)
                        buffered_rows = []

                # Reset for new qname
                current_qname = qname
                current_r1 = None
                current_r2 = None

            # Store reads
            if read.is_read1:
                if current_r1 is None:
                    current_r1 = read
            elif read.is_read2:
                if current_r2 is None:
                    current_r2 = read

        # Process final pair
        if current_qname is not None:
            pair_info = process_read_pair(
                current_qname, current_r1, current_r2, sorted_tags
            )
            buffered_rows.append(format_output_line(pair_info, sorted_tags))
            processed_qname_count += 1

        if buffered_rows:
            flush_batch(writer, buffered_rows, schema)

        samfile.close()

    logger.success("Processing complete!")
    logger.info(f"Total alignments scanned: {read_count:,}")
    logger.info(f"Total read pairs written: {processed_qname_count:,}")
    logger.info(f"Output saved to: {config.output_file}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Extract comprehensive summary for read pairs from BAM/SAM files"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Path to the input BAM/SAM file (must be qname-sorted)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Path to the output Parquet file",
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads for BAM decompression (default: 4)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable DEBUG level logging",
    )

    return parser.parse_args()


def main() -> int:
    """Main orchestrator for BAM to Parquet conversion."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = InputOutputConfig(
            input_bam=args.input,
            output_file=args.output,
            threads=args.threads,
        )

        # Execute BAM processing pipeline
        process_bam_file(config)

    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
