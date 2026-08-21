"""Target-sequence extraction and timepoint concatenation for insertion counts.

Concatenates per-timepoint insertion count files for a single sample/condition
into a wide, multi-indexed matrix, then annotates every insertion site with the
4 bp genomic target sequence taken from the reference genome.

Each timepoint file is loaded with a 3-level row MultiIndex (``Chr``,
coordinate, strand) and combined with an outer join so that insertions absent
from some timepoints are preserved. The concatenated matrix carries a 2-level
column MultiIndex ``(Timepoint, ReadType)``; the target sequence extracted
from the reference is appended as an extra ``Target`` row index level.

``concatenate_timepoint_data.py``'s ``save_processed_data`` is pure
orchestration plus file I/O with no separable domain logic to extract (per
Phase 4 convention 2): splitting the concatenated matrix by read type and
writing each to its own TSV. It dissolves entirely into that script's
``main()`` and contributes nothing here.

Input
-----
- Per-timepoint insertion count DataFrames and a reference-genome dict, as
  prepared by the invoking script.

Output
------
- ``AnalysisResult`` dataclass and the concatenated, target-annotated
  DataFrame.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from loguru import logger

# =============================================================================
# DATACLASSES
# =============================================================================
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
    input_files: list[Path],
    timepoints: list[str],
    ref_dict: dict,
) -> tuple[pd.DataFrame, AnalysisResult]:
    """Concatenate insertion data across multiple timepoints with target sequence annotation."""
    logger.info(f"Concatenating {len(timepoints)} timepoints")

    tp_files = {}
    for tp in timepoints:
        for file in input_files:
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
        num_timepoints=len(timepoints),
        num_insertions=len(concatenated),
        num_chromosomes=concatenated.index.get_level_values("Chr").nunique(),
        total_pbl_reads=total_pbl_reads,
        total_pbr_reads=total_pbr_reads,
        total_reads=total_reads,
    )

    return concatenated, result
