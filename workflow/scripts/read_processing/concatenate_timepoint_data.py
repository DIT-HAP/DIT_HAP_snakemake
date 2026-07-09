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
from Bio import SeqIO
from loguru import logger

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


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Summary statistics of the concatenation analysis."""
    num_timepoints: int
    num_insertions: int
    num_chromosomes: int
    total_pbl_reads: int
    total_pbr_reads: int
    total_reads: int

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
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
@logger.catch
def load_reference_data(genome_path: Path) -> dict:
    """Load reference genome sequences from FASTA file for target sequence extraction."""
    logger.info(f"Loading reference genome: {genome_path}")

    ref_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    logger.success(f"Loaded {len(ref_dict)} sequences from genome")

    # Log chromosome names
    chroms = list(ref_dict.keys())[:5]
    if len(ref_dict) > 5:
        logger.debug(f"Chromosomes: {chroms} ... and {len(ref_dict) - 5} more")
    else:
        logger.debug(f"Chromosomes: {chroms}")

    return ref_dict


@logger.catch
def extract_target_sequence(chrom: str, coordinate: int, ref_dict: dict) -> str:
    """Extract 4bp target sequence from reference genome at insertion coordinate."""
    try:
        if chrom not in ref_dict:
            logger.warning(f"Chromosome {chrom} not found in reference")
            return "NNNN"

        # Extract 4bp target sequence
        # Coordinate is 1-based, convert to 0-based
        start = coordinate - 4
        end = coordinate

        if start < 0:
            logger.warning(f"Coordinate {coordinate} too close to chromosome start")
            return "NNNN"

        seq = str(ref_dict[chrom].seq[start:end])

        if len(seq) != 4:
            logger.warning(f"Could not extract 4bp at {chrom}:{coordinate}")
            return "NNNN"

        return seq.upper()

    except Exception as e:
        logger.warning(f"Error extracting target at {chrom}:{coordinate}: {e}")
        return "NNNN"


@logger.catch
def process_concatenation_data(
    config: InputOutputConfig,
    ref_dict: dict,
) -> tuple[pd.DataFrame, AnalysisResult]:
    """Concatenate insertion data across multiple timepoints with target sequence annotation."""
    logger.info(f"Concatenating {len(config.timepoints)} timepoints")

    tp_files = {}
    for tp in config.timepoints:
        for file in config.input_files:
            if f"_{tp}_" in file.name:
                tp_files[tp] = file

    # Load all timepoint files
    dfs = {}
    for tp, file in tp_files.items():
        logger.debug(f"Loading timepoint {tp} from {file}")
        df = pd.read_csv(file, header=0, index_col=[0, 1, 2], sep="\t")
        logger.debug(f"  Loaded {len(df)} insertions for {tp}")
        dfs[tp] = df

    # Concatenate all timepoints
    logger.info("Concatenating dataframes...")
    concatenated = pd.concat(dfs, axis=1, join="outer")

    # Sort by timepoint names and coordinates
    concatenated = concatenated.sort_index(
        level=0, axis=1, key=lambda x: x.str.lower()
    ).sort_index(axis=0)

    logger.success(f"Concatenated {len(concatenated)} unique insertion sites")

    # Add target sequence information
    logger.info("Adding target sequences...")
    target_sequences = []

    for idx in concatenated.index:
        chrom = idx[0]
        coordinate = idx[1]
        target = extract_target_sequence(chrom, coordinate, ref_dict)
        target_sequences.append(target)

    # Add target as new index level
    concatenated = concatenated.set_index(
        pd.Series(target_sequences, name="Target", index=concatenated.index),
        append=True,
    )

    # Count unique targets
    unique_targets = concatenated.index.get_level_values("Target").unique()
    logger.info(f"Found {len(unique_targets)} unique target sequences")

    # Log target distribution if interesting
    target_counts = concatenated.index.get_level_values("Target").value_counts()
    if "TTAA" in target_counts.index:
        ttaa_fraction = target_counts["TTAA"] / len(concatenated) * 100
        logger.info(f"TTAA targets: {target_counts['TTAA']} ({ttaa_fraction:.1f}%)")

    # Calculate read totals before creating frozen stats object
    total_pbl_reads = 0
    total_pbr_reads = 0
    total_reads = 0

    if "PBL" in concatenated.columns.get_level_values(1):
        pbl_data = concatenated.xs("PBL", level=1, axis=1)
        total_pbl_reads = int(pbl_data.sum().sum())

    if "PBR" in concatenated.columns.get_level_values(1):
        pbr_data = concatenated.xs("PBR", level=1, axis=1)
        total_pbr_reads = int(pbr_data.sum().sum())

    if "Reads" in concatenated.columns.get_level_values(1):
        reads_data = concatenated.xs("Reads", level=1, axis=1)
        total_reads = int(reads_data.sum().sum())

    # Create statistics object with calculated values (can't modify after creation due to frozen=True)
    result = AnalysisResult(
        num_timepoints=len(config.timepoints),
        num_insertions=len(concatenated),
        num_chromosomes=concatenated.index.get_level_values("Chr").nunique(),
        total_pbl_reads=total_pbl_reads,
        total_pbr_reads=total_pbr_reads,
        total_reads=total_reads,
    )

    return concatenated, result


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
        concatenated, results = process_concatenation_data(config, ref_dict)

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

