"""Tests for read-count cutoff retention statistics.

The binned-distribution tests died with the pre-binned pipeline: the read count
figure now histograms raw values directly in the rendering layer.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
import pytest

from qc.read_counts import calculate_cutoff_statistics


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
