"""Tests for insertion-site extraction and counting."""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd

from read_processing.insertions import (
    calculate_insertion_coordinates_vectorized,
    count_insertions_vectorized,
    create_validation_mask,
)


# =============================================================================
# TESTS
# =============================================================================
def test_create_validation_mask_passes_valid_row():
    """A row with complete coordinates and a valid strand passes."""
    df = pd.DataFrame({
        "R1_Strand": ["+"],
        "R1_Chrom": ["chr1"],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [150],
    })
    mask = create_validation_mask(df)
    assert mask.tolist() == [True]


def test_create_validation_mask_fails_missing_strand():
    """A row with a missing strand is invalid."""
    df = pd.DataFrame({
        "R1_Strand": [None],
        "R1_Chrom": ["chr1"],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [150],
    })
    mask = create_validation_mask(df)
    assert mask.tolist() == [False]


def test_create_validation_mask_fails_missing_chrom():
    """A row with a missing chromosome is invalid."""
    df = pd.DataFrame({
        "R1_Strand": ["+"],
        "R1_Chrom": [None],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [150],
    })
    mask = create_validation_mask(df)
    assert mask.tolist() == [False]


def test_create_validation_mask_fails_missing_start():
    """A row with a missing reference start is invalid."""
    df = pd.DataFrame({
        "R1_Strand": ["+"],
        "R1_Chrom": ["chr1"],
        "R1_Ref_Start": [None],
        "R1_Ref_End": [150],
    })
    mask = create_validation_mask(df)
    assert mask.tolist() == [False]


def test_create_validation_mask_fails_missing_end():
    """A row with a missing reference end is invalid."""
    df = pd.DataFrame({
        "R1_Strand": ["+"],
        "R1_Chrom": ["chr1"],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [None],
    })
    mask = create_validation_mask(df)
    assert mask.tolist() == [False]


def test_create_validation_mask_fails_invalid_strand_value():
    """A row with a strand value other than + or - is invalid."""
    df = pd.DataFrame({
        "R1_Strand": ["0"],
        "R1_Chrom": ["chr1"],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [150],
    })
    mask = create_validation_mask(df)
    assert mask.tolist() == [False]


def test_calculate_insertion_coordinates_plus_strand():
    """+ strand coordinate is R1_Ref_Start + 4."""
    df = pd.DataFrame({
        "R1_Strand": ["+"],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [150],
    })
    coords = calculate_insertion_coordinates_vectorized(df)
    assert coords.tolist() == [104]


def test_calculate_insertion_coordinates_minus_strand():
    """- strand coordinate is R1_Ref_End."""
    df = pd.DataFrame({
        "R1_Strand": ["-"],
        "R1_Ref_Start": [100],
        "R1_Ref_End": [150],
    })
    coords = calculate_insertion_coordinates_vectorized(df)
    assert coords.tolist() == [150]


def test_count_insertions_aggregates_same_site():
    """Two reads landing on the same site aggregate into one count."""
    df = pd.DataFrame({
        "R1_Strand": ["+", "+"],
        "R1_Chrom": ["chr1", "chr1"],
        "R1_Ref_Start": [100, 100],
        "R1_Ref_End": [150, 160],
    })
    counts = count_insertions_vectorized(df)
    assert counts == {("chr1", 104): {"+": 2, "-": 0}}


def test_count_insertions_keeps_distinct_sites_separate():
    """Two distinct sites are counted independently."""
    df = pd.DataFrame({
        "R1_Strand": ["+", "-"],
        "R1_Chrom": ["chr1", "chr2"],
        "R1_Ref_Start": [100, 200],
        "R1_Ref_End": [150, 250],
    })
    counts = count_insertions_vectorized(df)
    assert counts == {
        ("chr1", 104): {"+": 1, "-": 0},
        ("chr2", 250): {"+": 0, "-": 1},
    }
