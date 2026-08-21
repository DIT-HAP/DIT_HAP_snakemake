"""Tests for genomic feature annotation of transposon insertion sites."""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pybedtools", reason="requires an env with pybedtools installed")

from read_processing.annotation import (
    assign_insertion_direction_vectorized,
    calculate_affected_residue_vectorized,
    calculate_codon_distances_vectorized,
    drop_boundary_duplicates,
)


# =============================================================================
# TESTS
# =============================================================================
def test_assign_insertion_direction_same_strand_is_forward():
    """A coding insertion on the same strand as its gene is Forward."""
    df = pd.DataFrame({
        "Type": ["Gene"],
        "Strand": ["+"],
        "Strand_Interval": ["+"],
    })
    direction = assign_insertion_direction_vectorized(df)
    assert direction.tolist() == ["Forward"]


def test_assign_insertion_direction_opposite_strand_is_reverse():
    """A coding insertion on the opposite strand from its gene is Reverse."""
    df = pd.DataFrame({
        "Type": ["Gene"],
        "Strand": ["+"],
        "Strand_Interval": ["-"],
    })
    direction = assign_insertion_direction_vectorized(df)
    assert direction.tolist() == ["Reverse"]


def test_assign_insertion_direction_intergenic_is_none():
    """An intergenic insertion has no direction."""
    df = pd.DataFrame({
        "Type": ["Intergenic region"],
        "Strand": ["+"],
        "Strand_Interval": ["+"],
    })
    direction = assign_insertion_direction_vectorized(df)
    assert direction.tolist() == [None]


def test_drop_boundary_duplicates_removes_boundary_row():
    """A row sitting exactly on the start or stop codon boundary is dropped."""
    sub_df = pd.DataFrame({
        "Type": ["Gene", "Gene"],
        "Distance_to_start_codon": [0, 5],
        "Distance_to_stop_codon": [10, 8],
    })
    result = drop_boundary_duplicates(sub_df)
    assert result["Distance_to_start_codon"].tolist() == [5]


def test_drop_boundary_duplicates_prefers_coding_over_intergenic():
    """A coding-vs-intergenic conflict at the same site is resolved by keeping coding."""
    sub_df = pd.DataFrame({
        "Type": ["Gene", "Intergenic region"],
        "Distance_to_start_codon": [5, np.nan],
        "Distance_to_stop_codon": [8, np.nan],
    })
    result = drop_boundary_duplicates(sub_df)
    assert result["Type"].tolist() == ["Gene"]
    assert len(result) == 1


def test_calculate_codon_distances_plus_strand():
    """+ strand codon distances map directly to region-start/end distances."""
    df = pd.DataFrame({
        "Type": ["Gene"],
        "Strand_Interval": ["+"],
        "Distance_to_region_start": [5],
        "Distance_to_region_end": [15],
        "Fraction_to_region_start": [0.25],
        "Fraction_to_region_end": [0.75],
    })
    codon_df = calculate_codon_distances_vectorized(df)
    assert codon_df["Distance_to_start_codon"].tolist() == [5]
    assert codon_df["Distance_to_stop_codon"].tolist() == [15]


def test_calculate_codon_distances_minus_strand():
    """- strand codon distances are reversed relative to region-start/end distances."""
    df = pd.DataFrame({
        "Type": ["Gene"],
        "Strand_Interval": ["-"],
        "Distance_to_region_start": [5],
        "Distance_to_region_end": [15],
        "Fraction_to_region_start": [0.25],
        "Fraction_to_region_end": [0.75],
    })
    codon_df = calculate_codon_distances_vectorized(df)
    assert codon_df["Distance_to_start_codon"].tolist() == [15]
    assert codon_df["Distance_to_stop_codon"].tolist() == [5]


def test_calculate_codon_distances_intergenic_is_nan():
    """Intergenic rows get NaN codon distances regardless of strand."""
    df = pd.DataFrame({
        "Type": ["Intergenic region"],
        "Strand_Interval": ["+"],
        "Distance_to_region_start": [5],
        "Distance_to_region_end": [15],
        "Fraction_to_region_start": [0.25],
        "Fraction_to_region_end": [0.75],
    })
    codon_df = calculate_codon_distances_vectorized(df)
    assert np.isnan(codon_df["Distance_to_start_codon"].iloc[0])
    assert np.isnan(codon_df["Distance_to_stop_codon"].iloc[0])


def test_calculate_affected_residue_cds_plus_strand():
    """A CDS insertion on the + strand computes residue index from its offset into the interval."""
    df = pd.DataFrame({
        "Type": ["Gene"],
        "Feature": ["CDS"],
        "Strand_Interval": ["+"],
        "Coordinate": [106],
        "Start_Interval": [100],
        "End_Interval": [200],
        "Accumulated_CDS_bases": [0.0],
    })
    residue_df = calculate_affected_residue_vectorized(df)
    # offset = 106 - 100 = 6; residue = 6 // 3 + 1 = 3; frame = 6 % 3 = 0
    assert residue_df["Residue_affected"].tolist() == [3.0]
    assert residue_df["Residue_frame"].tolist() == [0.0]


def test_calculate_affected_residue_intergenic_is_nan():
    """An intergenic row has no residue or frame."""
    df = pd.DataFrame({
        "Type": ["Intergenic region"],
        "Feature": ["intergenic"],
        "Strand_Interval": ["+"],
        "Coordinate": [106],
        "Start_Interval": [100],
        "End_Interval": [200],
        "Accumulated_CDS_bases": [0.0],
    })
    residue_df = calculate_affected_residue_vectorized(df)
    assert np.isnan(residue_df["Residue_affected"].iloc[0])
    assert np.isnan(residue_df["Residue_frame"].iloc[0])
