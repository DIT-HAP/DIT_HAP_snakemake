"""Per-read-pair extraction and formatting for BAM-to-Parquet conversion.

Extracts per-read-pair alignment summaries from pysam ``AlignedSegment``
objects and formats them into the flat row shape written to the output
Parquet table: read 1 and read 2 mapping quality, alignment length, CIGAR
string, strand, number of CIGAR operations, reference name, position,
reference start/end, SAM flag, a configurable set of alignment tags, and the
proper-pair status.

``parse_bam_to_tsv.py``'s ``process_bam_file`` is the streaming loop, reader,
and writer over the QNAME-sorted BAM/SAM file and stays in that script per
Phase 4 convention 3.

Input
-----
- A pysam ``AlignedSegment`` (or ``None``) for read 1 and read 2 of a query
  template, and the ordered list of SAM tags to extract.

Output
------
- ``ReadInfo`` and ``ReadPairInfo`` dataclasses, a formatted output row, the
  output column header, and the corresponding Parquet schema.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass, field

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
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
# DATACLASSES
# =============================================================================
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
