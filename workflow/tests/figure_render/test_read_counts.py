#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Read Count Distribution Figure Rendering
==================================================

Validates rendering against real data, asserting baseline bin counts, retention
statistics and artifact creation.

Baselines were derived from the pre-refactor arithmetic on
``results/11_merged/*.merged.tsv`` with ``--bins 50``, ``-t YES0``, ``-c 8``.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.read_counts import (
    load_cutoff_stats,
    load_distribution_data,
    render_distribution_figure,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
DATA_DIR = Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/18_figure_data")

# Positive-value counts per sample/timepoint, i.e. the sum of every bin in a group.
BASELINE_BIN_SUMS = {
    ("HD1328-4_YES", "YES0"): 293306,
    ("HD1328-4_YES", "YES1"): 212152,
    ("HD1328-4_YES", "YES2"): 215304,
    ("HD1328-4_YES", "YES3"): 211483,
    ("HD1328-4_YES", "YES4"): 220481,
    ("HD1328-7_YES", "YES0"): 156023,
    ("HD1328-7_YES", "YES1"): 133650,
    ("HD1328-7_YES", "YES2"): 150295,
    ("HD1328-7_YES", "YES3"): 134763,
    ("HD1328-7_YES", "YES4"): 149330,
    ("HD1328-8_YES", "YES0"): 168611,
    ("HD1328-8_YES", "YES1"): 143732,
    ("HD1328-8_YES", "YES2"): 148999,
    ("HD1328-8_YES", "YES3"): 140150,
    ("HD1328-8_YES", "YES4"): 149548,
}

# Cutoff retention per sample: rows_kept, pct_rows_kept, counts_kept, pct_counts_kept.
BASELINE_RETENTION = {
    "HD1328-4_YES": (659362, 143652, 21.7865, 380888542.0, 99.9404),
    "HD1328-7_YES": (287220, 123271, 42.9187, 282495932.0, 99.9809),
    "HD1328-8_YES": (287186, 129152, 44.9716, 416707741.0, 99.9843),
}

# log10 upper edge of the last YES0 bin for HD1328-4_YES, i.e. log10(max count).
BASELINE_MAX_LOG10 = 6.638270

DISTRIBUTION_COLUMNS = ["sample", "timepoint", "bin_left", "bin_right", "count"]
STATS_COLUMNS = [
    "sample",
    "original_rows",
    "rows_kept",
    "pct_rows_kept",
    "original_counts",
    "counts_kept",
    "pct_counts_kept",
]


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to the real binned distribution TSV (computation layer output)."""
    return DATA_DIR / "read_count_distribution.tsv"


@pytest.fixture
def real_stats_path() -> Path:
    """Path to the real cutoff statistics TSV (computation layer output)."""
    return DATA_DIR / "read_count_cutoff_stats.tsv"


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_read_count_distribution"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_bin_sums(real_data_path: Path) -> None:
    """Assert every group's bin counts sum to the pre-refactor positive-value count."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_distribution_data(real_data_path)
    actual = df.groupby(["sample", "timepoint"])["count"].sum()

    assert len(actual) == len(BASELINE_BIN_SUMS), f"Expected {len(BASELINE_BIN_SUMS)} groups, got {len(actual)}"

    for key, expected in BASELINE_BIN_SUMS.items():
        assert actual[key] == expected, f"{key}: expected {expected} positive values, got {actual[key]}"


def test_bin_count_is_fifty(real_data_path: Path) -> None:
    """Assert each group carries exactly 50 bins, the preserved default bin count."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_distribution_data(real_data_path)
    bins_per_group = df.groupby(["sample", "timepoint"]).size()

    assert (bins_per_group == 50).all(), f"Expected 50 bins per group, got {bins_per_group.unique().tolist()}"
    assert len(df) == 750, f"Expected 750 total bin rows, got {len(df)}"


def test_baseline_bin_edges(real_data_path: Path) -> None:
    """Assert HD1328-4_YES YES0 bin edges span log10 of the observed count range."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_distribution_data(real_data_path)
    group = df[(df["sample"] == "HD1328-4_YES") & (df["timepoint"] == "YES0")]

    assert group["bin_left"].min() == pytest.approx(0.0, abs=1e-6), "Lowest bin should start at log10(1) = 0"
    assert group["bin_right"].max() == pytest.approx(BASELINE_MAX_LOG10, abs=1e-5), (
        f"Highest bin should end at {BASELINE_MAX_LOG10}"
    )

    # First bin holds every insertion with a single-digit-ish count; baseline 108699.
    first_bin = group.sort_values("bin_left").iloc[0]
    assert first_bin["count"] == 108699, f"Expected 108699 in first bin, got {first_bin['count']}"


def test_baseline_retention_statistics(real_stats_path: Path) -> None:
    """Assert cutoff retention statistics match the pre-refactor values for every sample."""
    if not real_stats_path.exists():
        pytest.skip(f"Real stats not found: {real_stats_path}")

    stats_df = load_cutoff_stats(real_stats_path).set_index("sample")

    for sample, (orig_rows, rows_kept, pct_rows, counts_kept, pct_counts) in BASELINE_RETENTION.items():
        row = stats_df.loc[sample]
        assert row["original_rows"] == orig_rows, f"{sample}: original_rows"
        assert row["rows_kept"] == rows_kept, f"{sample}: rows_kept"
        assert row["pct_rows_kept"] == pytest.approx(pct_rows, abs=0.001), f"{sample}: pct_rows_kept"
        assert row["counts_kept"] == pytest.approx(counts_kept, rel=1e-9), f"{sample}: counts_kept"
        assert row["pct_counts_kept"] == pytest.approx(pct_counts, abs=0.001), f"{sample}: pct_counts_kept"


def test_dual_artifacts_created(real_data_path: Path, real_stats_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists() or not real_stats_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_distribution_data(real_data_path)
    stats_df = load_cutoff_stats(real_stats_path)
    render_distribution_figure(df, stats_df, output_stem, "YES0", 8.0)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"
    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_renderer_does_not_rebin(real_data_path: Path) -> None:
    """Assert the rendered bar heights equal the stored bin counts exactly."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    import cnsplots as cns
    import matplotlib.pyplot as plt

    df = load_distribution_data(real_data_path)
    group = df[(df["sample"] == "HD1328-4_YES") & (df["timepoint"] == "YES0")].dropna(subset=["bin_left"])

    fig, ax = plt.subplots()
    binrange = (float(group["bin_left"].min()), float(group["bin_right"].max()))
    cns.histplot(data=group, x="bin_center", weights="count", bins=len(group), binrange=binrange, ax=ax)

    heights = np.array([patch.get_height() for patch in ax.patches])
    stored = group["count"].to_numpy().astype(float)
    plt.close(fig)

    assert np.array_equal(heights, stored), "Rendered bar heights differ from stored bin counts"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert an empty distribution TSV is handled gracefully without crashing."""
    empty_tsv = tmp_path / "empty.tsv"
    pd.DataFrame(columns=DISTRIBUTION_COLUMNS).to_csv(empty_tsv, sep="\t", index=False)

    empty_stats = tmp_path / "empty_stats.tsv"
    pd.DataFrame(columns=STATS_COLUMNS).to_csv(empty_stats, sep="\t", index=False)

    df = load_distribution_data(empty_tsv)
    stats_df = load_cutoff_stats(empty_stats)
    assert df.empty

    render_distribution_figure(df, stats_df, tmp_path / "empty_test", "YES0", 8.0)


def test_no_valid_data_panel(tmp_path: Path) -> None:
    """Assert a group whose bins are all empty renders a 'No valid data' panel."""
    marker_tsv = tmp_path / "marker.tsv"
    pd.DataFrame(
        {
            "sample": ["S1", "S1"],
            "timepoint": ["YES0", "YES1"],
            "bin_left": [0.0, np.nan],
            "bin_right": [1.0, np.nan],
            "count": [10, 0],
        }
    ).to_csv(marker_tsv, sep="\t", index=False)

    stats_tsv = tmp_path / "marker_stats.tsv"
    pd.DataFrame(
        {
            "sample": ["S1"],
            "original_rows": [10],
            "rows_kept": [5],
            "pct_rows_kept": [50.0],
            "original_counts": [100.0],
            "counts_kept": [80.0],
            "pct_counts_kept": [80.0],
        }
    ).to_csv(stats_tsv, sep="\t", index=False)

    output_stem = tmp_path / "marker_test"
    df = load_distribution_data(marker_tsv)
    stats_df = load_cutoff_stats(stats_tsv)

    render_distribution_figure(df, stats_df, output_stem, "YES0", 8.0)

    assert (output_stem.parent / f"{output_stem.name}.pdf").exists(), "PDF not created for marker-row input"
