"""Filtered-alignment BAM export helpers.

Reconstruct a filtered BAM from an existing coordinate-sorted BAM and the
filtered Parquet produced by ``filter_aligned_reads``. The Parquet supplies the
authoritative query-name whitelist (one row per read pair that passed
filtering); the BAM supplies the original records, so SEQ, QUAL, SAM tags, and
flags are preserved exactly. Unmapped, secondary, and supplementary records
(FLAG 0x4 / 0x100 / 0x800) are skipped because ``parse_bam_to_tsv`` excluded
them from the Parquet, and each retained record is tagged with the read group
of its source (timepoint + PBL/PBR fragment).

``export_filtered_bam.py``'s ``main`` is the CLI/I/O shell; the functions here
are the shared implementation.

Input
-----
- Filtered Parquet with a ``QueryName`` column.
- Coordinate-sorted BAM containing the original alignments.

Output
------
- ``ExportStats`` summary and one unsorted BAM written with ``@SQ``/``@RG``.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-10
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path

import pyarrow.parquet as pq
import pysam

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
QUERY_NAME_COLUMN = "QueryName"
QUERY_NAME_BATCH_SIZE = 1_000_000

# FLAG bits that exclude a record from the exported BAM:
#   0x4   unmapped
#   0x100 secondary alignment
#   0x800 supplementary (chimeric/split) alignment
EXCLUDED_FLAG_MASK = 0x4 | 0x100 | 0x800


# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ExportStats:
    """Record counts from one BAM export run."""

    pbl_pairs: int
    pbr_pairs: int
    pbl_records: int
    pbr_records: int

    @property
    def total_records(self) -> int:
        """Total records written to the output BAM."""
        return self.pbl_records + self.pbr_records

    @property
    def expected_records(self) -> int:
        """Expected records when every filtered pair has both mates in the BAM."""
        return 2 * (self.pbl_pairs + self.pbr_pairs)


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
def count_parquet_pairs(parquet_path: Path) -> int:
    """Return the number of read-pair rows in a filtered Parquet file."""
    with pq.ParquetFile(parquet_path) as parquet_file:
        return parquet_file.metadata.num_rows


def load_query_names(parquet_path: Path) -> set[str]:
    """Load the QueryName whitelist from a filtered Parquet file."""
    query_names: set[str] = set()

    with pq.ParquetFile(parquet_path) as parquet_file:
        for batch in parquet_file.iter_batches(
            batch_size=QUERY_NAME_BATCH_SIZE,
            columns=[QUERY_NAME_COLUMN],
        ):
            query_names.update(batch.column(0).to_pylist())

    return query_names


def build_export_header(
    source_bam: Path,
    read_groups: list[dict[str, str]],
) -> pysam.AlignmentHeader:
    """Copy the BAM header and append the export read groups."""
    with pysam.AlignmentFile(str(source_bam), "rb") as source:
        header_dict = source.header.to_dict()

    header_dict.setdefault("HD", {})["SO"] = "unsorted"
    header_dict.setdefault("RG", []).extend(read_groups)
    return pysam.AlignmentHeader.from_dict(header_dict)


def write_filtered_records(
    output_bam: Path,
    header: pysam.AlignmentHeader,
    sources: list[tuple[Path, set[str], str]],
    threads: int,
) -> list[int]:
    """Write whitelisted primary records from each source BAM, tagging the read group."""
    counts: list[int] = []

    with pysam.AlignmentFile(str(output_bam), "wb", header=header) as output:
        for source_bam, query_names, read_group in sources:
            written = 0
            with pysam.AlignmentFile(str(source_bam), "rb", threads=threads) as source:
                for read in source:
                    if read.flag & EXCLUDED_FLAG_MASK:
                        continue
                    if read.query_name not in query_names:
                        continue

                    read.set_tag("RG", read_group)
                    output.write(read)
                    written += 1

            counts.append(written)

    return counts
