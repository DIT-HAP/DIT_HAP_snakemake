#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Insertion Density Figure Rendering
=============================================

Validates rendering against real data, asserting baseline statistics and
artifact creation.

Baselines were derived by running the pre-refactor
``insertion_density_analysis.py`` on
``results/13_filtered/raw_reads.filtered.tsv`` /
``results/12_concatenated/annotations.tsv`` with initial/final timepoints
YES0/YES4, and confirmed byte-identical after the refactor via diff.

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

from figure_render.density import load_density_data, render_density_figure

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# Point at the path rule insertion_density_data actually writes. Asserting
# against the old reports/ copy would let these tests pass even if the
# refactored computation rule broke, since that artifact predates the split.
REAL_DATA_PATH = Path(
    "/data/c/yangyusheng_optimized/DIT_HAP_snakemake"
    "/projects/HD_DIT_HAP/results/18_figure_data/insertion_density_analysis.tsv"
)

# Baseline from the pre-refactor script, YES0/YES4, confirmed byte-identical post-refactor.
BASELINE_TOTAL_GENES = 4936
BASELINE_TOTAL_INSERTIONS = 71469
BASELINE_MEAN_DENSITY_INITIAL = 6.250128646677472
BASELINE_MEAN_DENSITY_FINAL = 5.9771130470016205
BASELINE_MEAN_DENSITY_LOG2FC = -0.0616185170178282
BASELINE_MEAN_GINI_DEPTH_INITIAL = 0.5084866288492706
BASELINE_MEAN_STRAND_BIAS = 0.07037115072933549

REQUIRED_COLUMNS = [
    "insertion_density_per_kb_initial",
    "insertion_density_per_kb_final",
    "insertion_density_log2fc",
    "total_reads_initial",
    "total_reads_final",
    "gini_coefficient_of_depth_initial",
    "gini_coefficient_of_depth_final",
]


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to the real density statistics TSV (computation layer output)."""
    return REAL_DATA_PATH


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_insertion_density"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_gene_count(real_data_path: Path) -> None:
    """Assert the total gene count matches the pre-refactor baseline."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_density_data(real_data_path)
    assert len(df) == BASELINE_TOTAL_GENES, f"Expected {BASELINE_TOTAL_GENES} genes, got {len(df)}"


def test_baseline_total_insertions(real_data_path: Path) -> None:
    """Assert total insertions analyzed matches the pre-refactor baseline."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_density_data(real_data_path)
    total = int(df["total_insertions"].sum())
    assert total == BASELINE_TOTAL_INSERTIONS, f"Expected {BASELINE_TOTAL_INSERTIONS} insertions, got {total}"


def test_baseline_mean_metrics(real_data_path: Path) -> None:
    """Assert mean density, gini, and strand bias statistics match the pre-refactor baseline."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_density_data(real_data_path)

    assert df["insertion_density_per_kb_initial"].mean() == pytest.approx(BASELINE_MEAN_DENSITY_INITIAL, abs=1e-6)
    assert df["insertion_density_per_kb_final"].mean() == pytest.approx(BASELINE_MEAN_DENSITY_FINAL, abs=1e-6)
    assert df["insertion_density_log2fc"].mean() == pytest.approx(BASELINE_MEAN_DENSITY_LOG2FC, abs=1e-6)
    assert df["gini_coefficient_of_depth_initial"].mean() == pytest.approx(BASELINE_MEAN_GINI_DEPTH_INITIAL, abs=1e-6)
    assert df["strand_bias"].mean() == pytest.approx(BASELINE_MEAN_STRAND_BIAS, abs=1e-6)


def test_required_columns_present(real_data_path: Path) -> None:
    """Assert the density statistics table carries every column the figure draws."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_density_data(real_data_path)
    for col in REQUIRED_COLUMNS:
        assert col in df.columns, f"Missing required column: {col}"


def test_dual_artifacts_created(real_data_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_density_data(real_data_path)
    render_density_figure(df, output_stem, "YES0", "YES4")

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"
    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert an empty density statistics TSV is handled gracefully without crashing."""
    empty_tsv = tmp_path / "empty.tsv"
    pd.DataFrame(columns=[*REQUIRED_COLUMNS, "FYPOviability"]).to_csv(empty_tsv, sep="\t", index=False)

    df = load_density_data(empty_tsv)
    assert df.empty

    render_density_figure(df, tmp_path / "empty_test", "YES0", "YES4")

    pdf_path = tmp_path / "empty_test.pdf"
    assert pdf_path.exists(), "PDF artifact not created for empty data"


def test_missing_required_column_returns_none(tmp_path: Path) -> None:
    """Assert a density statistics TSV missing a required column is rejected (logger.catch swallows the raise)."""
    bad_tsv = tmp_path / "bad.tsv"
    pd.DataFrame({"insertion_density_per_kb_initial": [1.0, 2.0]}).to_csv(bad_tsv, sep="\t", index=False)

    # load_density_data is decorated with @logger.catch, matching the established
    # convention in plot_pbl_pbr_correlation.py / plot_read_count_distribution.py:
    # the ValueError is logged and the function returns None rather than propagating.
    result = load_density_data(bad_tsv)
    assert result is None
