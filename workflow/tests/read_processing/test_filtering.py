"""Tests for threshold-based filtering of aligned read pairs."""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import pandas as pd

from read_processing.filtering import (
    FilterThresholds,
    build_filter_mask,
    coerce_column_dtypes,
    strip_read_prefix,
)


# =============================================================================
# TESTS
# =============================================================================
def test_strip_read_prefix_strips_r1():
    """An R1_ prefix is stripped."""
    assert strip_read_prefix("R1_MAPQ") == "MAPQ"


def test_strip_read_prefix_strips_r2():
    """An R2_ prefix is stripped."""
    assert strip_read_prefix("R2_NM") == "NM"


def test_strip_read_prefix_passes_through_unprefixed():
    """A column without an R1_/R2_ prefix is returned unchanged."""
    assert strip_read_prefix("Chr") == "Chr"


def test_coerce_column_dtypes_on_synthetic_string_chunk():
    """An all-string chunk is coerced to the dtypes pd.read_csv would infer."""
    chunk = pd.DataFrame({
        "R1_MAPQ": ["30"],
        "R2_MAPQ": ["25"],
        "R1_NM": ["2"],
        "R2_NM": ["N/A"],
        "R1_SA": ["N/A"],
        "R2_SA": [""],
        "R1_XA": ["chr2,100,+,50M,0;"],
        "R2_XA": ["NA"],
        "Other_Col": ["hello"],
    })
    coerced = coerce_column_dtypes(chunk)

    assert coerced["R1_MAPQ"].dtype == np.dtype("int64")
    assert coerced["R1_MAPQ"].iloc[0] == 30
    assert coerced["R2_MAPQ"].iloc[0] == 25

    assert coerced["R1_NM"].dtype == np.dtype("float64")
    assert coerced["R1_NM"].iloc[0] == 2.0
    assert coerced["R2_NM"].dtype == np.dtype("float64")
    assert pd.isna(coerced["R2_NM"].iloc[0])

    assert pd.isna(coerced["R1_SA"].iloc[0])
    assert pd.isna(coerced["R2_SA"].iloc[0])
    assert coerced["R1_XA"].iloc[0] == "chr2,100,+,50M,0;"
    assert pd.isna(coerced["R2_XA"].iloc[0])

    assert coerced["Other_Col"].iloc[0] == "hello"


def test_build_filter_mask_row_inside_every_threshold_and_each_failure():
    """Row 0 passes every threshold; each other row fails exactly one."""
    chunk = pd.DataFrame({
        "R1_MAPQ": [40, 10, 40, 40, 40, 40, 40],
        "R2_MAPQ": [40, 40, 40, 40, 40, 40, 40],
        "R1_NCIGAR": [1, 1, 5, 1, 1, 1, 1],
        "R2_NCIGAR": [1, 1, 1, 1, 1, 1, 1],
        "R1_NM": [0, 0, 0, 10, 0, 0, 0],
        "R2_NM": [0, 0, 0, 0, 0, 0, 0],
        "R1_SA": [np.nan, np.nan, np.nan, np.nan, "chr1,100,+,50M,60,0;", np.nan, np.nan],
        "R2_SA": [np.nan] * 7,
        "R1_XA": [np.nan, np.nan, np.nan, np.nan, np.nan, "chr2,200,-,50M,2;", np.nan],
        "R2_XA": [np.nan] * 7,
        "Is_Proper_Pair": ["Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "No"],
    })
    r1_filters = FilterThresholds(
        mapq_threshold=30, ncigar_value=2, nm_threshold=3, no_sa=True, no_xa=True,
    )
    r2_filters = FilterThresholds()

    mask = build_filter_mask(chunk, r1_filters, r2_filters, require_proper_pair=True)

    assert mask.tolist() == [True, False, False, False, False, False, False]
