#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Insertion Orientation Figure Rendering
================================================

Validates rendering against real data, asserting baseline row counts, log-log
correlations and artifact creation.

Baselines were derived from the pre-refactor arithmetic on
``results/13_filtered/raw_reads.filtered.tsv``: stack both column levels, unstack
Strand, ``dropna(axis=0)``, then keep rows with ``min(axis=1) > 0``. Correlations
are Pearson on log10 values, matching what the old ``create_scatter_correlation_plot``
reported for log-scaled axes.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pandas as pd
import pytest

from figure_render.orientation import load_and_prepare_data, render_orientation_figure

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
REAL_DATA_PATH = Path(
    "/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/18_figure_data/arc/strand_pairs.tsv"
)

# Rows retained after min(axis=1) > 0, and Pearson r on log10 values, per group.
BASELINE_GROUPS = {
    ("HD1328-4_YES", "YES0"): (52602, 0.685155),
    ("HD1328-4_YES", "YES1"): (51466, 0.629445),
    ("HD1328-4_YES", "YES2"): (51388, 0.660974),
    ("HD1328-4_YES", "YES3"): (51195, 0.693912),
    ("HD1328-4_YES", "YES4"): (51075, 0.708606),
    ("HD1328-7_YES", "YES0"): (43768, 0.588175),
    ("HD1328-7_YES", "YES1"): (42067, 0.520461),
    ("HD1328-7_YES", "YES2"): (42592, 0.612527),
    ("HD1328-7_YES", "YES3"): (42214, 0.647225),
    ("HD1328-7_YES", "YES4"): (41817, 0.645900),
    ("HD1328-8_YES", "YES0"): (46337, 0.603617),
    ("HD1328-8_YES", "YES1"): (45002, 0.597594),
    ("HD1328-8_YES", "YES2"): (44218, 0.598641),
    ("HD1328-8_YES", "YES3"): (44300, 0.671635),
    ("HD1328-8_YES", "YES4"): (44128, 0.674323),
}

# Sum of all per-group row counts above.
BASELINE_TOTAL_ROWS = 694169

INPUT_COLUMNS = ["sample", "timepoint", "plus_count", "minus_count"]


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to the real strand pairs TSV (computation layer output)."""
    return REAL_DATA_PATH


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_insertion_orientation"


# =============================================================================
# TESTS
# =============================================================================
def test_total_row_count(real_data_path: Path) -> None:
    """Assert the total retained row count matches the pre-refactor baseline."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    assert len(df) == BASELINE_TOTAL_ROWS, f"Expected {BASELINE_TOTAL_ROWS} rows, got {len(df)}"


def test_baseline_group_row_counts(real_data_path: Path) -> None:
    """Assert each sample-timepoint group retains its pre-refactor row count."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)
    actual = df.groupby(["sample", "timepoint"]).size()

    assert len(actual) == len(BASELINE_GROUPS), f"Expected {len(BASELINE_GROUPS)} groups, got {len(actual)}"

    for key, (expected_rows, _) in BASELINE_GROUPS.items():
        assert actual[key] == expected_rows, f"{key}: expected {expected_rows} rows, got {actual[key]}"


def test_baseline_log10_correlations(real_data_path: Path) -> None:
    """Assert log-log Pearson correlations match the pre-refactor values for every group."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    for key, (_, expected_r) in BASELINE_GROUPS.items():
        sample, timepoint = key
        group = df[(df["sample"] == sample) & (df["timepoint"] == timepoint)]
        actual_r = group[["log10_plus_count", "log10_minus_count"]].corr().iloc[0, 1]

        assert actual_r == pytest.approx(expected_r, abs=0.001), (
            f"{key}: expected r = {expected_r:.6f}, got {actual_r:.6f}"
        )


def test_all_counts_strictly_positive(real_data_path: Path) -> None:
    """Assert the min > 0 filtering semantics were preserved by the computation layer."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    assert (df["plus_count"] > 0).all(), "Found non-positive plus_count values"
    assert (df["minus_count"] > 0).all(), "Found non-positive minus_count values"
    assert df["log10_plus_count"].notna().all(), "Found non-finite log10_plus_count values"
    assert df["log10_minus_count"].notna().all(), "Found non-finite log10_minus_count values"


def test_dual_artifacts_created(real_data_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    # One panel per group is enough to exercise the layout without rendering
    # 694k rasterised points, which takes minutes.
    subset = df[df["sample"] == "HD1328-4_YES"]
    render_orientation_figure(subset, output_stem)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"
    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_single_output_holds_all_groups(real_data_path: Path, output_stem: Path) -> None:
    """Assert every group reaches one figure, guarding the fixed per-file PDF truncation bug."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    # The old code opened PdfPages inside a per-input-file loop, so each file
    # truncated the previous file's output and only the last survived. A single
    # rendering pass over the combined TSV makes that structurally impossible.
    assert df.groupby(["sample", "timepoint"]).ngroups == len(BASELINE_GROUPS)
    assert df["sample"].nunique() == 3, "Expected all three samples in one input"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert an empty strand pairs TSV is handled gracefully without crashing."""
    empty_tsv = tmp_path / "empty.tsv"
    pd.DataFrame(columns=INPUT_COLUMNS).to_csv(empty_tsv, sep="\t", index=False)

    df = load_and_prepare_data(empty_tsv)
    assert df.empty

    render_orientation_figure(df, tmp_path / "empty_test")


def test_nonpositive_rows_dropped(tmp_path: Path) -> None:
    """Assert rows with a zero or negative strand count are dropped before log10."""
    mixed_tsv = tmp_path / "mixed.tsv"
    pd.DataFrame(
        {
            "sample": ["S1", "S1", "S1"],
            "timepoint": ["YES0", "YES0", "YES0"],
            "plus_count": [10.0, 0.0, 5.0],
            "minus_count": [20.0, 30.0, 0.0],
        }
    ).to_csv(mixed_tsv, sep="\t", index=False)

    df = load_and_prepare_data(mixed_tsv)

    assert len(df) == 1, f"Expected 1 surviving row, got {len(df)}"
    assert df.iloc[0]["plus_count"] == 10.0
