#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "biopython",
#     "loguru",
#     "pandas",
# ]
# ///

"""
Concatenate Timepoint Insertion Counts with Target-Sequence Annotation
======================================================================

Concatenate per-timepoint insertion count files for a single sample/condition
into a wide, multi-indexed matrix, then annotate every insertion site with the
4 bp genomic target sequence taken from the reference genome.

Each timepoint file is matched to its timepoint label by the ``_<timepoint>_``
token embedded in the file name, loaded with a 3-level row MultiIndex
(``Chr``, coordinate, strand), and combined with an outer join so that
insertions absent from some timepoints are preserved. The concatenated matrix
carries a 2-level column MultiIndex ``(Timepoint, ReadType)``; the target
sequence extracted from the reference is appended as an extra ``Target`` row
index level. The matrix is finally split by read type (PBL, PBR, Reads) and
written to three separate tab-separated files with missing counts filled as 0.

Input
-----
- Multiple per-timepoint insertion count files (TSV) with a 3-level row
  MultiIndex on columns 0-2 and read-type columns (PBL / PBR / Reads).
- A reference genome FASTA file used to extract the 4 bp target sequence.

Output
------
- Three tab-separated matrices (PBL, PBR, Reads), each written with the full
  row index and header preserved and missing counts filled as integer 0.

Usage
-----
    python concatenate_timepoint_data.py -s sample_cond -i t1.tsv t2.tsv -tp T1 T2 -g genome.fasta -ol out.PBL.tsv -or out.PBR.tsv -o out.Reads.tsv
    python concatenate_timepoint_data.py -s sample_cond -i t1.tsv t2.tsv -tp T1 T2 -g genome.fasta -ol out.PBL.tsv -or out.PBR.tsv -o out.Reads.tsv --verbose

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

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from read_processing.sequence_extraction import (  # noqa: E402
    AnalysisResult,
    load_reference_data,
    extract_target_sequence,
    process_concatenation_data,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Configuration and validated paths for timepoint concatenation."""
    input_files: list[Path]
    genome_file: Path
    output_pbl: Path
    output_pbr: Path
    output_reads: Path
    sample_name: str
    timepoints: list[str]

    def __post_init__(self) -> None:
        """Validate input/genome paths and timepoint count, then create output dirs."""
        missing = [f for f in self.input_files if not f.exists()]
        if missing:
            raise ValueError(f"Input files not found: {missing}")
        if not self.genome_file.exists():
            raise ValueError(f"Genome file not found: {self.genome_file}")
        if len(self.timepoints) != len(self.input_files):
            raise ValueError(
                f"Number of timepoints ({len(self.timepoints)}) must match number of input files "
                f"({len(self.input_files)})"
            )
        for output_path in (self.output_pbl, self.output_pbr, self.output_reads):
            output_path.parent.mkdir(parents=True, exist_ok=True)


setup_logger()

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def save_processed_data(concatenated: pd.DataFrame, config: InputOutputConfig) -> None:
    """Save concatenated timepoint data to separate output files for each read type."""
    logger.info("Saving concatenated data...")

    # Save PBL data
    if "PBL" in concatenated.columns.get_level_values(1):
        pbl_data = concatenated.xs("PBL", level=1, axis=1)
        pbl_data.fillna(0).astype(int).to_csv(
            config.output_pbl, index=True, sep="\t"
        )
        logger.success(f"Saved PBL data to {config.output_pbl}")
    else:
        logger.warning("No PBL data found in concatenated results")

    # Save PBR data
    if "PBR" in concatenated.columns.get_level_values(1):
        pbr_data = concatenated.xs("PBR", level=1, axis=1)
        pbr_data.fillna(0).astype(int).to_csv(
            config.output_pbr, index=True, sep="\t"
        )
        logger.success(f"Saved PBR data to {config.output_pbr}")
    else:
        logger.warning("No PBR data found in concatenated results")

    # Save Reads data
    if "Reads" in concatenated.columns.get_level_values(1):
        reads_data = concatenated.xs("Reads", level=1, axis=1)
        reads_data.fillna(0).astype(int).to_csv(
            config.output_reads, index=True, sep="\t"
        )
        logger.success(f"Saved Reads data to {config.output_reads}")
    else:
        logger.warning("No Reads data found in concatenated results")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Concatenate insertion data across timepoints with target sequence annotation"
    )
    parser.add_argument("-s", "--sample", type=str, required=True, help="Sample name")
    parser.add_argument("-i", "--input", type=Path, nargs="+", required=True, help="Path to the input insertion count files")
    parser.add_argument("-tp", "--timepoints", type=str, nargs="+", required=True, help="Timepoint names")
    parser.add_argument("-g", "--genome", type=Path, required=True, help="Reference genome FASTA file")
    parser.add_argument("-ol", "--output_pbl", type=Path, required=True, help="Output PBL file path")
    parser.add_argument("-or", "--output_pbr", type=Path, required=True, help="Output PBR file path")
    parser.add_argument("-o", "--output_reads", type=Path, required=True, help="Output reads file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: concatenate timepoint insertion data with target-sequence annotation."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        # Validate input and output paths using the config dataclass
        config = InputOutputConfig(
            input_files=args.input,
            timepoints=args.timepoints,
            genome_file=args.genome,
            output_pbl=args.output_pbl,
            output_pbr=args.output_pbr,
            output_reads=args.output_reads,
            sample_name=args.sample,
        )

        logger.info(f"Starting processing of {config.sample_name}")

        # Load reference genome
        ref_dict = load_reference_data(config.genome_file)

        # Run the core analysis/logic
        concatenated, results = process_concatenation_data(
            config.input_files, config.timepoints, ref_dict
        )

        # Save results
        save_processed_data(concatenated, config)

        logger.success(
            f"Analysis complete. Results saved to {config.output_pbl}, "
            f"{config.output_pbr}, and {config.output_reads}"
        )

        # Display summary statistics
        logger.info("=" * 60)
        logger.info("CONCATENATION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Sample: {config.sample_name}")
        logger.info(f"Timepoints: {', '.join(config.timepoints)}")
        logger.success(f"Unique insertion sites: {results.num_insertions:,}")
        logger.info(f"Chromosomes: {results.num_chromosomes}")
        if results.total_reads > 0:
            logger.info("\nRead counts:")
            logger.info(f"  PBL reads: {results.total_pbl_reads:,}")
            logger.info(f"  PBR reads: {results.total_pbr_reads:,}")
            logger.info(f"  Total reads: {results.total_reads:,}")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

