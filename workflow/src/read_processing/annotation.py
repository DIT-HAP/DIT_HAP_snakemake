"""Genomic feature annotation of transposon insertion sites.

Annotates transposon insertion sites with genomic features including genes,
intergenic regions, and coding sequences. Insertion coordinates are intersected
against a genome-region BED annotation using pybedtools, then distances to the
start/stop codon, the affected amino-acid residue, the reading frame, and the
insertion direction relative to each gene are computed.

Boundary duplicates (insertions falling exactly on a region edge, or matching
both a coding and an intergenic region) are resolved so that each insertion
keeps a single, most-specific annotation.

Input
-----
- Insertion file (TSV or CSV) with columns: Chr, Coordinate, Strand, Target.
- Genome-region BED file (TSV) with region intervals and their metadata.

Output
------
- Tab-separated annotation table with per-insertion genomic features, codon
  distances, affected residues, and insertion direction.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger
from pybedtools import BedTool

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Aggregate statistics describing the annotation outcome."""
    total_insertions: int
    annotated_insertions: int
    coding_insertions: int
    intergenic_insertions: int
    unique_genes: int
    forward_insertions: int
    reverse_insertions: int

    @property
    def coding_percentage(self) -> float:
        """Percentage of insertions falling in coding regions."""
        if self.total_insertions == 0:
            return 0.0
        return (self.coding_insertions / self.total_insertions) * 100.0

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_insertion_data(input_path: Path) -> pd.DataFrame:
    """Load insertion sites from a TSV or CSV file into a BED-ready DataFrame."""
    logger.info(f"Loading insertion data from {input_path}")

    file_ext = input_path.suffix.lower()

    if file_ext == ".tsv":
        df = pd.read_csv(input_path, header=0, usecols=[0, 1, 2, 3], sep="\t")
    elif file_ext == ".csv":
        df = pd.read_csv(input_path, header=0, usecols=[0, 1, 2, 3])
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    # Rename coordinate column for BED processing
    df = df.rename(columns={"Coordinate": "End"})
    df.insert(1, "Start", df["End"])

    logger.success(f"Loaded {len(df)} insertion sites")
    logger.debug("For BED processing, Start and End are the same (point insertions)")

    return df


@logger.catch
def load_genome_regions(region_path: Path) -> pd.DataFrame:
    """Load genome region annotations from a BED file into a DataFrame."""
    logger.info(f"Loading genome regions from {region_path}")

    df = pd.read_csv(region_path, sep="\t", header=0)
    df = df.rename(columns={"#Chr": "Chr"})

    logger.success(f"Loaded {len(df)} genome regions")

    # Log region types if available
    if "Type" in df.columns:
        type_counts = df["Type"].value_counts()
        logger.debug("Region types:")
        for region_type, count in type_counts.items():
            logger.debug(f"  {region_type}: {count}")

    return df


def calculate_codon_distances_vectorized(annotated_df: pd.DataFrame) -> pd.DataFrame:
    """Compute strand-aware distances/fractions to start and stop codons using vectorized operations."""
    # Initialize output columns
    codon_df = pd.DataFrame(index=annotated_df.index)

    # Mask for non-intergenic
    non_intergenic = annotated_df["Type"] != "Intergenic region"
    plus_strand = annotated_df["Strand_Interval"] == "+"

    # For + strand: start=region_start, stop=region_end; for - strand: reverse
    codon_df["Distance_to_start_codon"] = np.where(
        non_intergenic & plus_strand,
        annotated_df["Distance_to_region_start"],
        np.where(
            non_intergenic & ~plus_strand,
            annotated_df["Distance_to_region_end"],
            np.nan
        )
    )
    codon_df["Distance_to_stop_codon"] = np.where(
        non_intergenic & plus_strand,
        annotated_df["Distance_to_region_end"],
        np.where(
            non_intergenic & ~plus_strand,
            annotated_df["Distance_to_region_start"],
            np.nan
        )
    )
    codon_df["Fraction_to_start_codon"] = np.where(
        non_intergenic & plus_strand,
        annotated_df["Fraction_to_region_start"],
        np.where(
            non_intergenic & ~plus_strand,
            annotated_df["Fraction_to_region_end"],
            np.nan
        )
    )
    codon_df["Fraction_to_stop_codon"] = np.where(
        non_intergenic & plus_strand,
        annotated_df["Fraction_to_region_end"],
        np.where(
            non_intergenic & ~plus_strand,
            annotated_df["Fraction_to_region_start"],
            np.nan
        )
    )

    return codon_df


def calculate_affected_residue_vectorized(annotated_df: pd.DataFrame) -> pd.DataFrame:
    """Compute the affected amino-acid residue index and reading frame using vectorized operations."""
    residue_df = pd.DataFrame(index=annotated_df.index)

    # Mask for rows where calculation applies
    non_coding_mask = (annotated_df["Type"] == "Intergenic region") | (annotated_df["Type"] == "Non-coding gene")

    # Get accumulated CDS base (default 0 if missing), convert to float
    cds_base = pd.to_numeric(
        annotated_df.get("Accumulated_CDS_bases", pd.Series(0.0, index=annotated_df.index)),
        errors='coerce'
    ).fillna(0.0)

    # Add offset for CDS features
    is_cds = annotated_df["Feature"] == "CDS"
    plus_strand = annotated_df["Strand_Interval"] == "+"

    cds_offset = np.where(
        is_cds & plus_strand,
        annotated_df["Coordinate"].astype(int) - annotated_df["Start_Interval"].astype(int),
        np.where(
            is_cds & ~plus_strand,
            annotated_df["End_Interval"].astype(int) - annotated_df["Coordinate"].astype(int),
            0
        )
    )

    cds_base_total = cds_base + cds_offset

    # Calculate residue and frame
    residue_df["Residue_affected"] = np.where(
        non_coding_mask,
        np.nan,
        (cds_base_total // 3 + 1)
    )
    residue_df["Residue_frame"] = np.where(
        non_coding_mask,
        np.nan,
        (cds_base_total % 3)
    )

    return residue_df


def assign_insertion_direction_vectorized(annotated_df: pd.DataFrame) -> pd.Series:
    """Determine the insertion direction (Forward/Reverse/NaN) relative to the gene using vectorized operations."""
    intergenic = annotated_df["Type"] == "Intergenic region"
    same_strand = annotated_df["Strand"] == annotated_df["Strand_Interval"]

    # Use None instead of np.nan to avoid dtype promotion error with str in numpy 2.x
    return pd.Series(
        np.where(
            intergenic,
            None,
            np.where(same_strand, "Forward", "Reverse")
        ),
        index=annotated_df.index,
        dtype=object
    )


@logger.catch
def drop_boundary_duplicates(sub_df: pd.DataFrame) -> pd.DataFrame:
    """Remove boundary-duplicated insertions, preferring coding over intergenic."""
    # Remove insertions exactly at boundaries
    boundary_mask = (
        (sub_df["Distance_to_start_codon"] == 0) |
        (sub_df["Distance_to_stop_codon"] == 0)
    )
    sub_df = sub_df[~boundary_mask]

    # If only one type remains, return it
    if sub_df["Type"].nunique() == 1:
        return sub_df

    # Otherwise, prefer coding regions over intergenic
    coding_mask = sub_df["Type"] != "Intergenic region"
    if coding_mask.any():
        return sub_df[coding_mask]

    return sub_df


def annotate_insertions(
    insertions_df: pd.DataFrame,
    regions_df: pd.DataFrame,
) -> tuple[pd.DataFrame, AnalysisResult]:
    """Intersect insertions with genome regions and derive per-insertion feature annotations."""
    logger.info("Annotating insertions with genomic features...")

    # Convert to BedTool objects
    insertions_bed = BedTool.from_dataframe(insertions_df)
    regions_bed = BedTool.from_dataframe(regions_df)

    # Get column names for merging
    insertion_cols = insertions_df.columns.tolist()
    region_cols = [
        f"{col}_Interval" if col in insertion_cols else col
        for col in regions_df.columns
    ]

    logger.debug(f"Insertion columns: {len(insertion_cols)}")
    logger.debug(f"Region columns: {len(region_cols)}")

    # Intersect insertions with regions
    intersected = insertions_bed.intersect(regions_bed, wa=True, wb=True)

    # Convert back to DataFrame
    annotated_df = intersected.to_dataframe(names=insertion_cols + region_cols)

    # Clean up DataFrame
    annotated_df.drop(columns=["Start"], inplace=True)
    annotated_df.rename(columns={"End": "Coordinate"}, inplace=True)

    # Replace "." with NaN
    annotated_df.replace(r"^\.$", np.nan, inplace=True, regex=True)

    logger.info("Calculating distances and annotations...")

    # Calculate distances to region boundaries
    annotated_df["Distance_to_region_start"] = (
        annotated_df["Coordinate"] - annotated_df["ParentalRegion_start"]
    )
    annotated_df["Distance_to_region_end"] = (
        annotated_df["ParentalRegion_end"] - annotated_df["Coordinate"]
    )
    annotated_df["Fraction_to_region_start"] = (
        annotated_df["Distance_to_region_start"] /
        annotated_df["ParentalRegion_length"]
    )
    annotated_df["Fraction_to_region_end"] = (
        annotated_df["Distance_to_region_end"] /
        annotated_df["ParentalRegion_length"]
    )

    # Calculate codon distances (vectorized)
    codon_distances = calculate_codon_distances_vectorized(annotated_df)
    annotated_df = pd.concat([annotated_df, codon_distances], axis=1)

    # Calculate affected residues (vectorized)
    residue_info = calculate_affected_residue_vectorized(annotated_df)
    annotated_df = pd.concat([annotated_df, residue_info], axis=1)

    # Assign insertion direction (vectorized)
    annotated_df["Insertion_direction"] = assign_insertion_direction_vectorized(annotated_df)

    logger.info("Removing boundary duplicates...")

    # Drop duplicated insertions at boundaries. Iterate manually rather than
    # groupby().apply() so the Chr/Coordinate/Strand grouping columns are not
    # dropped from the result (pandas >=2.2 excludes them from the group frame
    # passed to apply(), and drop_boundary_duplicates() never re-adds them).
    deduped_groups = []
    for _, group in annotated_df.groupby(["Chr", "Coordinate", "Strand"]):
        deduped_groups.append(drop_boundary_duplicates(group))
    annotated_df = pd.concat(deduped_groups).reset_index(drop=True)

    # Calculate statistics
    stats = AnalysisResult(
        total_insertions=len(insertions_df),
        annotated_insertions=len(annotated_df),
        coding_insertions=len(annotated_df[annotated_df["Type"] != "Intergenic region"]),
        intergenic_insertions=len(annotated_df[annotated_df["Type"] == "Intergenic region"]),
        unique_genes=annotated_df["GeneName"].nunique() if "GeneName" in annotated_df.columns else 0,
        forward_insertions=len(annotated_df[annotated_df["Insertion_direction"] == "Forward"]),
        reverse_insertions=len(annotated_df[annotated_df["Insertion_direction"] == "Reverse"]),
    )

    return annotated_df, stats
