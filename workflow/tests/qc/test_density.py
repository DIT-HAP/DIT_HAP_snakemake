"""Tests for insertion density, gap, and strand statistics."""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
import pytest

from qc.density import (
    calculate_gap_statistics,
    calculate_gini_coefficient,
    calculate_strand_statistics,
)


# =============================================================================
# TESTS: calculate_gini_coefficient
# =============================================================================
def test_gini_coefficient_perfectly_even_array_is_zero():
    """Equal values give zero inequality."""
    values = np.array([5, 5, 5, 5])
    assert calculate_gini_coefficient(values) == pytest.approx(0.0, abs=1e-12)


def test_gini_coefficient_maximally_skewed_array_approaches_one():
    """One nonzero value among many zeros drives Gini close to 1.

    For n values where only the last (largest) is nonzero, the closed-form
    result is (n-1)/n, e.g. n=100 -> 0.99.
    """
    values = np.array([0] * 99 + [100])
    result = calculate_gini_coefficient(values)
    assert result == pytest.approx(0.99, abs=1e-12)
    assert result > 0.95


def test_gini_coefficient_single_element_is_zero():
    """A single-element array always returns 0.0.

    With n=1 the cumsum short-circuit only fires for a zero value; for a
    nonzero single value, cumsum[-1] == the value itself, and the formula
    reduces to (2*1*v)/(1*v) - (1+1)/1 = 2 - 2 = 0. So both [0] and [5]
    (or any single value) yield exactly 0.0 -- there is only one data point,
    so there is no inequality to measure between points.
    """
    assert calculate_gini_coefficient(np.array([5])) == pytest.approx(0.0, abs=1e-12)
    assert calculate_gini_coefficient(np.array([0])) == pytest.approx(0.0, abs=1e-12)


def test_gini_coefficient_empty_array_is_zero():
    """An empty array returns 0.0 via the explicit len==0 guard."""
    assert calculate_gini_coefficient(np.array([])) == 0.0


def test_gini_coefficient_all_zero_short_circuits_to_zero():
    """All-zero values (e.g. a fully-depleted gene) return 0.0, not NaN from 0/0."""
    assert calculate_gini_coefficient(np.array([0, 0, 0])) == 0.0


# =============================================================================
# TESTS: calculate_gap_statistics
# =============================================================================
def _make_gene_insertions(coordinates: list[int], start: int, end: int, length: int) -> pd.DataFrame:
    """Build a minimal gene_insertions frame indexed by (Chr, Coordinate, Strand, Target)."""
    index = pd.MultiIndex.from_tuples(
        [("chr1", coord, "+", "geneA") for coord in coordinates],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    return pd.DataFrame(
        {
            "ParentalRegion_start": [start] * len(coordinates),
            "ParentalRegion_end": [end] * len(coordinates),
            "ParentalRegion_length": [length] * len(coordinates),
        },
        index=index,
    )


def test_gap_statistics_two_insertions_known_gap():
    """Two insertions at coordinates 120 and 150, gene start=end=coords, yield one gap.

    Boundary set = sorted({120, 150}) (start/end coincide with the insertion
    coordinates themselves). Gap = 150 - 120 - 1 = 29 (the count of bases
    strictly between the two insertion sites). With only one gap, min == max
    == mean == median == 29 and sd == 0. Fraction = 29 / length = 29 / 200 = 0.145.
    """
    gene_insertions = _make_gene_insertions([120, 150], start=120, end=150, length=200)
    result = calculate_gap_statistics(gene_insertions)

    assert result["num_gaps"] == 1
    assert result["largest_gap"] == 29
    assert result["smallest_gap"] == 29
    assert result["mean_gap_length"] == pytest.approx(29.0)
    assert result["median_gap_length"] == pytest.approx(29.0)
    assert result["gap_length_sd"] == pytest.approx(0.0)
    assert result["largest_gap_fraction"] == pytest.approx(0.145)
    assert result["all_gap_lengths"] == "29"
    # A single gap value carries no spread, so its own Gini is 0.
    assert result["gini_coefficient_of_location"] == pytest.approx(0.0)


def test_gap_statistics_single_insertion_no_gap_possible():
    """A single insertion sitting exactly at start==end==coordinate has zero gaps.

    Reading the code: coordinates_with_start_and_end collapses to a single
    unique value {120}, so the gaps list comprehension iterates zero times.
    ``gaps`` is empty, taking the early-return branch, which reports all
    zero-valued numeric fields, empty-string joined fields, and NaN for the
    location Gini coefficient (documented, not guessed).
    """
    gene_insertions = _make_gene_insertions([120], start=120, end=120, length=50)
    result = calculate_gap_statistics(gene_insertions)

    assert result["num_gaps"] == 0
    assert result["largest_gap"] == 0
    assert result["smallest_gap"] == 0
    assert result["mean_gap_length"] == 0
    assert result["all_gap_lengths"] == ""
    assert result["all_gap_lengths_fraction"] == ""
    assert np.isnan(result["gini_coefficient_of_location"])


# =============================================================================
# TESTS: calculate_strand_statistics
# =============================================================================
def test_strand_statistics_three_forward_one_reverse():
    """3 Forward / 1 Reverse insertions give a 0.75/0.25 split and 0.25 bias.

    forward_count=3, reverse_count=1, total=4 ->
    forward_preference = 3/4 = 0.75, reverse_preference = 1/4 = 0.25,
    strand_bias = |0.75 - 0.5| = 0.25. All four coordinates are distinct
    (100, 110, 120, 130), so no coordinate has both a Forward and a Reverse
    insertion -> paired_sites = 0, paired_sites_fraction = 0/4 = 0.0.
    """
    index = pd.MultiIndex.from_tuples(
        [
            ("chr1", 100, "+", "geneA"),
            ("chr1", 110, "+", "geneA"),
            ("chr1", 120, "+", "geneA"),
            ("chr1", 130, "+", "geneA"),
        ],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    gene_insertions = pd.DataFrame(
        {"Insertion_direction": ["Forward", "Forward", "Forward", "Reverse"]},
        index=index,
    )

    result = calculate_strand_statistics(gene_insertions)

    assert result["forward_insertions"] == 3
    assert result["reverse_insertions"] == 1
    assert result["forward_preference"] == pytest.approx(0.75)
    assert result["reverse_preference"] == pytest.approx(0.25)
    assert result["strand_bias"] == pytest.approx(0.25)
    assert result["paired_sites"] == 0
    assert result["paired_sites_fraction"] == pytest.approx(0.0)


def test_strand_statistics_paired_site_at_same_coordinate():
    """A coordinate with both Forward and Reverse insertions counts as one paired site."""
    index = pd.MultiIndex.from_tuples(
        [
            ("chr1", 100, "+", "geneA"),
            ("chr1", 100, "-", "geneA"),
            ("chr1", 110, "+", "geneA"),
        ],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )
    gene_insertions = pd.DataFrame(
        {"Insertion_direction": ["Forward", "Reverse", "Forward"]},
        index=index,
    )

    result = calculate_strand_statistics(gene_insertions)

    # total_sites counts unique Coordinate values (100, 110) -> 2 sites.
    assert result["paired_sites"] == 1
    assert result["paired_sites_fraction"] == pytest.approx(0.5)
