"""Tests for read-count cutoff retention and binned distribution computation."""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
import pytest

from qc.read_counts import calculate_cutoff_statistics, compute_binned_distribution


# =============================================================================
# TESTS: calculate_cutoff_statistics
# =============================================================================
def test_cutoff_statistics_straddling_cutoff_uses_greater_equal():
    """Rows with values [5, 10, 10, 15] and cutoff=10 keep the >= 10 rows.

    The comparison in the code is ``df[initial_time_point] >= cutoff``, so a
    value exactly equal to the cutoff (10) is retained, not dropped:
    kept = {10, 10, 15} -> rows_kept = 3 of 4 -> pct_rows_kept = 75.0.
    original_counts = 5+10+10+15 = 40; counts_kept = 10+10+15 = 35 ->
    pct_counts_kept = 35/40*100 = 87.5.
    """
    df = pd.DataFrame({"T0": [5, 10, 10, 15]})
    result = calculate_cutoff_statistics(df, "T0", 10)

    assert result["original_rows"] == 4
    assert result["rows_kept"] == 3
    assert result["pct_rows_kept"] == pytest.approx(75.0)
    assert result["original_counts"] == pytest.approx(40.0)
    assert result["counts_kept"] == pytest.approx(35.0)
    assert result["pct_counts_kept"] == pytest.approx(87.5)


def test_cutoff_statistics_value_just_below_cutoff_is_dropped():
    """A value one unit below the cutoff is excluded, confirming the >= boundary."""
    df = pd.DataFrame({"T0": [9, 10]})
    result = calculate_cutoff_statistics(df, "T0", 10)

    assert result["rows_kept"] == 1
    assert result["pct_rows_kept"] == pytest.approx(50.0)
    assert result["counts_kept"] == pytest.approx(10.0)


def test_cutoff_statistics_missing_column_returns_none():
    """A missing initial timepoint column raises ValueError internally.

    ``calculate_cutoff_statistics`` is decorated with ``@logger.catch``, which
    matches the established convention elsewhere in this codebase (see
    ``workflow/tests/figure_render/test_density.py``): the ValueError is
    logged and swallowed, and the function returns None rather than
    propagating the exception to the caller.
    """
    df = pd.DataFrame({"T1": [1, 2]})
    result = calculate_cutoff_statistics(df, "T0", 10)
    assert result is None


def test_cutoff_statistics_non_numeric_column_returns_none():
    """A non-numeric initial timepoint column also raises internally and is swallowed to None."""
    df = pd.DataFrame({"T0": ["a", "b"]})
    result = calculate_cutoff_statistics(df, "T0", 10)
    assert result is None


# =============================================================================
# TESTS: compute_binned_distribution
# =============================================================================
def test_binned_distribution_known_values_bins_4():
    """log10([1, 10, 100, 1000]) = [0, 1, 2, 3]; 4 equal-width bins span [0, 3].

    numpy.histogram splits [0, 3] into 4 bins of width 0.75:
    [0, 0.75), [0.75, 1.5), [1.5, 2.25), [2.25, 3.0] (last bin closed on the
    right). Each of the four log values (0, 1, 2, 3) falls in a distinct bin,
    so every bin has exactly one count.
    """
    df = pd.DataFrame({"T0": [1, 10, 100, 1000]})
    result = compute_binned_distribution(df, "sampleA", bins=4)

    assert list(result.columns) == ["sample", "timepoint", "bin_left", "bin_right", "count"]
    assert len(result) == 4
    assert result["sample"].tolist() == ["sampleA"] * 4
    assert result["timepoint"].tolist() == ["T0"] * 4
    assert result["count"].tolist() == [1, 1, 1, 1]
    np.testing.assert_allclose(result["bin_left"].tolist(), [0.0, 0.75, 1.5, 2.25])
    np.testing.assert_allclose(result["bin_right"].tolist(), [0.75, 1.5, 2.25, 3.0])


def test_binned_distribution_ignores_non_positive_values():
    """Zero and negative values are dropped before the log10/histogram step.

    Only the two positive values (1, 100) participate; log10 -> [0, 2],
    binned into bins=2 of width 1.0: [0, 1) and [1, 2], one count each.
    """
    df = pd.DataFrame({"T0": [1, 0, -5, 100]})
    result = compute_binned_distribution(df, "sampleA", bins=2)

    assert result["count"].tolist() == [1, 1]
    np.testing.assert_allclose(result["bin_left"].tolist(), [0.0, 1.0])
    np.testing.assert_allclose(result["bin_right"].tolist(), [1.0, 2.0])


def test_binned_distribution_no_positive_values_emits_empty_marker_row():
    """A column with no positive values emits a single NaN-bin, zero-count marker row.

    Reading the code: when ``positive`` is empty, the function short-circuits
    before calling np.histogram and appends a one-row frame with bin_left and
    bin_right set to NaN and count 0, so the sample/timepoint group is not
    silently dropped from the output.
    """
    df = pd.DataFrame({"T0": [0, -1, -2]})
    result = compute_binned_distribution(df, "sampleA", bins=4)

    assert len(result) == 1
    assert result["count"].tolist() == [0]
    assert np.isnan(result["bin_left"].iloc[0])
    assert np.isnan(result["bin_right"].iloc[0])


def test_binned_distribution_no_numeric_columns_returns_none():
    """A frame with no numeric columns raises ValueError internally, swallowed by @logger.catch to None."""
    df = pd.DataFrame({"T0": ["a", "b", "c"]})
    result = compute_binned_distribution(df, "sampleA", bins=4)
    assert result is None
