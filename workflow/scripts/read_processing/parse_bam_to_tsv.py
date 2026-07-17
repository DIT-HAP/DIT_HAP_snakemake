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
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

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

# Base (non-tag) numeric columns emitted by format_output_line as native ints.
# These are always integers straight from pysam (mapq/n_cigar/flag default to 0;
# length/pos/ref_start/ref_end are nullable -> stored as int64 with real nulls).
# SAM tag columns are NOT here: a tag can be non-numeric, so tags stay strings.
# Names must match build_header exactly (note R1_Flag vs R2_FLAG).
INT64_COLUMNS = frozenset(
    {
        "R1_MAPQ", "R1_LEN", "R1_NCIGAR", "R1_Pos", "R1_Ref_Start", "R1_Ref_End", "R1_Flag",
        "R2_MAPQ", "R2_LEN", "R2_NCIGAR", "R2_Pos", "R2_Ref_Start", "R2_Ref_End", "R2_FLAG",
    }
)

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


@dataclass(kw_only=True, slots=True, frozen=True)
class ReadInfo:
    """Alignment information extracted from a single pysam AlignedSegment."""
    mapq: int = 0
    length: int | None = None
    cigar: str = "N/A"
    strand: str = "N/A"
    n_cigar: int = 0
    chrom: str = "N/A"
    pos: int | None = None
    ref_start: int | None = None
    ref_end: int | None = None
    flag: int = 0
    tags: dict[str, str] = field(default_factory=dict)


@dataclass(kw_only=True, slots=True, frozen=True)
class ReadPairInfo:
    """Paired-end read information container for a single query template."""
    qname: str
    read1: ReadInfo | None = None
    read2: ReadInfo | None = None
    is_proper_pair: str = "N/A"

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for BAM to Parquet conversion."""
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
# NOTE: the per-read-pair functions below are intentionally NOT decorated with
# @logger.catch. It is called ~4x per read pair (~24M times on a full lane), and
# its wrapper overhead is measurable; more importantly, on a computation function
# @logger.catch swallows exceptions and returns None instead of propagating,
# which silently corrupts downstream output (see the batch-3 annotation fix in
# docs/optimization_plan.md). @logger.catch stays only on the I/O-boundary
# process_bam_file, where a real error should be logged and abort the run.
def extract_read_info(read: pysam.AlignedSegment | None, tag_list: list[str]) -> ReadInfo:
    """Extract alignment information from a pysam AlignedSegment (or a placeholder if None)."""
    if read is None:
        return ReadInfo(tags={tag: "N/A" for tag in tag_list})

    # Extract basic information
    strand = "-" if read.is_reverse else "+"

    # Calculate NCIGAR
    n_cigar = 0
    if read.cigartuples:
        n_cigar = len(read.cigartuples)
    elif read.cigarstring:
        n_cigar = 0

    # Extract tags
    read_tags = dict(read.get_tags())
    formatted_tags: dict[str, str] = {}

    for tag_name in tag_list:
        value = read_tags.get(tag_name, "N/A")
        if value is True:
            formatted_tags[tag_name] = "True"
        elif value is False:
            formatted_tags[tag_name] = "False"
        elif value is None:
            formatted_tags[tag_name] = "N/A"
        else:
            formatted_tags[tag_name] = str(value)

    return ReadInfo(
        mapq=read.mapping_quality if read.mapping_quality is not None else 0,
        length=read.query_alignment_length,
        cigar=read.cigarstring if read.cigarstring else "N/A",
        strand=strand,
        n_cigar=n_cigar,
        chrom=read.reference_name if read.reference_name is not None else "N/A",
        pos=read.reference_start
        if read.reference_start is not None and read.reference_start != -1
        else None,
        ref_start=read.reference_start if read.reference_start is not None else None,
        ref_end=read.reference_end if read.reference_end is not None else None,
        flag=read.flag if read.flag is not None else 0,
        tags=formatted_tags,
    )


def determine_proper_pair_status(
    read1: pysam.AlignedSegment | None, read2: pysam.AlignedSegment | None
) -> str:
    """Determine whether the reads form a proper pair."""
    if read1 and read1.is_paired:
        return "Yes" if read1.is_proper_pair else "No"
    elif read2 and read2.is_paired:
        return "Yes" if read2.is_proper_pair else "No"
    elif (read1 and not read1.is_paired) or (read2 and not read2.is_paired):
        return "Single_End_Or_Flag_Issue"
    return "N/A"


def process_read_pair(
    qname: str,
    read1: pysam.AlignedSegment | None,
    read2: pysam.AlignedSegment | None,
    tag_list: list[str],
) -> ReadPairInfo:
    """Process a read pair and return its consolidated information."""
    is_proper_pair = determine_proper_pair_status(read1, read2)

    r1_info = extract_read_info(read1, tag_list)
    r2_info = extract_read_info(read2, tag_list)

    return ReadPairInfo(
        qname=qname, read1=r1_info, read2=r2_info, is_proper_pair=is_proper_pair
    )


def format_output_line(pair_info: ReadPairInfo, tag_list: list[str]) -> list:
    """Format ReadPairInfo into an ordered list of output fields (native int for numeric base fields, str/None otherwise)."""
    output_line = [pair_info.qname]

    # Read 1 information. Numeric base fields are emitted as native int (or None
    # for the nullable ones) so the Parquet writer can encode them as int64
    # directly; strings/tags stay as-is. build_schema declares the matching type.
    r1 = pair_info.read1 or ReadInfo(tags={tag: "N/A" for tag in tag_list})
    output_line.extend(
        [
            r1.mapq,
            r1.length,
            r1.cigar,
            r1.strand,
            r1.n_cigar,
            r1.chrom,
            r1.pos,
            r1.ref_start,
            r1.ref_end,
            r1.flag,
        ]
    )
    for tag in tag_list:
        output_line.append(r1.tags.get(tag, "N/A"))

    # Read 2 information
    r2 = pair_info.read2 or ReadInfo(tags={tag: "N/A" for tag in tag_list})
    output_line.extend(
        [
            r2.mapq,
            r2.length,
            r2.cigar,
            r2.strand,
            r2.n_cigar,
            r2.chrom,
            r2.pos,
            r2.ref_start,
            r2.ref_end,
            r2.flag,
        ]
    )
    for tag in tag_list:
        output_line.append(r2.tags.get(tag, "N/A"))

    output_line.append(pair_info.is_proper_pair)
    return output_line


def build_header(tag_list: list[str]) -> list[str]:
    """Build the ordered list of output column names."""
    header_fields = [
        "QueryName",
        "R1_MAPQ",
        "R1_LEN",
        "R1_CIGAR",
        "R1_Strand",
        "R1_NCIGAR",
        "R1_Chrom",
        "R1_Pos",
        "R1_Ref_Start",
        "R1_Ref_End",
        "R1_Flag",
    ]

    for tag_name in tag_list:
        header_fields.append(f"R1_{tag_name}")

    header_fields.extend(
        [
            "R2_MAPQ",
            "R2_LEN",
            "R2_CIGAR",
            "R2_Strand",
            "R2_NCIGAR",
            "R2_Chrom",
            "R2_Pos",
            "R2_Ref_Start",
            "R2_Ref_End",
            "R2_FLAG",
        ]
    )

    for tag_name in tag_list:
        header_fields.append(f"R2_{tag_name}")

    header_fields.append("Is_Proper_Pair")
    return header_fields


def build_schema(header_fields: list[str]) -> pa.Schema:
    """Build a Parquet schema: base numeric fields as int64, all other columns as string."""
    return pa.schema(
        [
            (name, pa.int64() if name in INT64_COLUMNS else pa.string())
            for name in header_fields
        ]
    )


def flush_batch(
    writer: pq.ParquetWriter, buffered_rows: list[list], schema: pa.Schema
) -> None:
    """Write a batch of buffered rows to the Parquet file as one row group, typing each column per the schema."""
    columns = list(zip(*buffered_rows))
    arrays = [
        pa.array(column, type=schema.field(index).type)
        for index, column in enumerate(columns)
    ]
    table = pa.Table.from_arrays(arrays, schema=schema)
    writer.write_table(table)


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
