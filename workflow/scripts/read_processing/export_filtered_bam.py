"""
Export Filtered Alignments as BAM
=================================

Reconstruct a BAM for the read pairs that passed ``filter_aligned_reads``.
The filtered Parquet supplies the authoritative ``QueryName`` whitelist; the
coordinate-sorted BAM (``results/4_sorted``) supplies the original alignment
records, so SEQ, QUAL, SAM tags, and flags are preserved. PBL and PBR records
are combined into one BAM, each tagged with its own ``@RG``. The output is
written unsorted; the Snakemake rule sorts and indexes it with samtools.

Input
-----
- PBL/PBR filtered Parquet files with a ``QueryName`` column.
- PBL/PBR coordinate-sorted BAM files with the original alignments.

Output
------
- One unsorted BAM containing the whitelisted primary records of both
  fragments, with read groups ``{sample}_{timepoint}_{condition}_{PBL,PBR}``.

Usage
-----
    python export_filtered_bam.py \
        --pbl-parquet PBL.filtered.parquet --pbr-parquet PBR.filtered.parquet \
        --pbl-bam PBL.sorted.bam --pbr-bam PBR.sorted.bam \
        --sample S --timepoint T --condition C \
        -o output.unsorted.bam --threads 4

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-11-10
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = (SCRIPT_DIR / "../../src").resolve()
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from logging_setup import setup_logger
from read_processing.bam_export import (
    ExportStats,
    build_export_header,
    count_parquet_pairs,
    load_query_names,
    write_filtered_records,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
DEFAULT_THREADS = 4


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class ExportConfig:
    """Paths and identifiers for one sample/timepoint/condition BAM export."""

    pbl_parquet: Path
    pbr_parquet: Path
    pbl_bam: Path
    pbr_bam: Path
    sample: str
    timepoint: str
    condition: str
    output_bam: Path
    threads: int = DEFAULT_THREADS

    def __post_init__(self) -> None:
        """Validate inputs and ensure the output directory exists."""
        for path in (self.pbl_parquet, self.pbr_parquet, self.pbl_bam, self.pbr_bam):
            if not path.exists():
                raise ValueError(f"Input file not found: {path}")
        if not 1 <= self.threads <= 32:
            raise ValueError(f"threads must be between 1 and 32, got {self.threads}")
        self.output_bam.parent.mkdir(parents=True, exist_ok=True)

    @property
    def read_group_base(self) -> str:
        """Read-group prefix shared by both fragments."""
        return f"{self.sample}_{self.timepoint}_{self.condition}"

    @property
    def read_groups(self) -> list[dict[str, str]]:
        """Read-group header records for the PBL and PBR fragments."""
        return [
            {"ID": f"{self.read_group_base}_PBL", "SM": self.sample},
            {"ID": f"{self.read_group_base}_PBR", "SM": self.sample},
        ]


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
def run_export(config: ExportConfig) -> ExportStats:
    """Load both whitelists and write the combined unsorted BAM."""
    pbl_names = load_query_names(config.pbl_parquet)
    pbr_names = load_query_names(config.pbr_parquet)
    logger.info(f"PBL whitelist: {len(pbl_names):,} query names")
    logger.info(f"PBR whitelist: {len(pbr_names):,} query names")

    header = build_export_header(config.pbl_bam, config.read_groups)
    sources = [
        (config.pbl_bam, pbl_names, f"{config.read_group_base}_PBL"),
        (config.pbr_bam, pbr_names, f"{config.read_group_base}_PBR"),
    ]

    counts = write_filtered_records(config.output_bam, header, sources, config.threads)
    return ExportStats(
        pbl_pairs=count_parquet_pairs(config.pbl_parquet),
        pbr_pairs=count_parquet_pairs(config.pbr_parquet),
        pbl_records=counts[0],
        pbr_records=counts[1],
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Export filtered alignments from coordinate-sorted BAMs as one unsorted BAM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--pbl-parquet",
        type=Path,
        required=True,
        help="PBL filtered Parquet file.",
    )
    parser.add_argument(
        "--pbr-parquet",
        type=Path,
        required=True,
        help="PBR filtered Parquet file.",
    )
    parser.add_argument(
        "--pbl-bam",
        type=Path,
        required=True,
        help="PBL coordinate-sorted BAM file.",
    )
    parser.add_argument(
        "--pbr-bam",
        type=Path,
        required=True,
        help="PBR coordinate-sorted BAM file.",
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="Sample name for the read group SM field.",
    )
    parser.add_argument(
        "--timepoint",
        type=str,
        required=True,
        help="Timepoint name for the read group ID.",
    )
    parser.add_argument(
        "--condition",
        type=str,
        required=True,
        help="Condition name for the read group ID.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output unsorted BAM file.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="BAM decompression threads.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable DEBUG level logging.",
    )
    return parser.parse_args()


def main() -> int:
    """Create the export configuration, run the export, and report record counts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    try:
        config = ExportConfig(
            pbl_parquet=args.pbl_parquet,
            pbr_parquet=args.pbr_parquet,
            pbl_bam=args.pbl_bam,
            pbr_bam=args.pbr_bam,
            sample=args.sample,
            timepoint=args.timepoint,
            condition=args.condition,
            output_bam=args.output,
            threads=args.threads,
        )
        stats = run_export(config)

        logger.info("=" * 60)
        logger.info("BAM EXPORT SUMMARY")
        logger.info("=" * 60)
        logger.info(f"PBL: {stats.pbl_pairs:,} pairs -> {stats.pbl_records:,} records")
        logger.info(f"PBR: {stats.pbr_pairs:,} pairs -> {stats.pbr_records:,} records")
        logger.success(f"Total records written: {stats.total_records:,}")

        if stats.total_records != stats.expected_records:
            logger.warning(
                f"Record count {stats.total_records:,} differs from 2 x pair count "
                f"{stats.expected_records:,}; some filtered pairs may lack a mate in the BAM."
            )

        logger.success(f"Saved unsorted BAM to {config.output_bam}")
    except Exception as e:  # noqa: BLE001 - CLI entry point logs and exits non-zero
        logger.exception(f"BAM export failed: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
