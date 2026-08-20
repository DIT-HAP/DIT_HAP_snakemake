#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Gene Coverage Figure Rendering
==========================================

Validates rendering against real data, asserting baseline coverage statistics
and artifact creation.

Baselines were derived by running the pre-refactor ``gene_coverage_analysis.py``
on ``results/14_insertion_level_depletion_analysis/LFC.tsv`` /
``results/12_concatenated/annotations.tsv`` /
``resources/pombase_data/2026-06-01/Gene_metadata/gene_viability.tsv``.

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

from figure_render.coverage import load_coverage_data, render_coverage_figure

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
REAL_DATA_PATH = Path(
    "/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/18_figure_data/arc/gene_coverage_stats.tsv"
)

# Baseline from the pre-refactor script: category -> (covered, not_covered, total, coverage_pct)
BASELINE_STATS = {
    "viable": (3185, 423, 3608, 88.276),
    "inviable": (1090, 117, 1207, 90.307),
    "condition-dependent": (160, 15, 175, 91.429),
    "unknown": (111, 7584, 7695, 1.442),
}

REQUIRED_COLUMNS = ["category", "covered", "not_covered", "total", "coverage_pct"]


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to the real coverage statistics TSV (computation layer output)."""
    return REAL_DATA_PATH


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_gene_coverage"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_category_count(real_data_path: Path) -> None:
    """Assert the number of viability categories matches the pre-refactor baseline."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_coverage_data(real_data_path)
    assert len(df) == len(BASELINE_STATS), f"Expected {len(BASELINE_STATS)} categories, got {len(df)}"


def test_baseline_coverage_stats(real_data_path: Path) -> None:
    """Assert covered/not_covered/total/coverage_pct match the pre-refactor baseline per category."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_coverage_data(real_data_path).set_index("category")

    for category, (covered, not_covered, total, coverage_pct) in BASELINE_STATS.items():
        row = df.loc[category]
        assert int(row["covered"]) == covered, f"{category}: expected covered={covered}, got {row['covered']}"
        assert int(row["not_covered"]) == not_covered, f"{category}: expected not_covered={not_covered}, got {row['not_covered']}"
        assert int(row["total"]) == total, f"{category}: expected total={total}, got {row['total']}"
        assert row["coverage_pct"] == pytest.approx(coverage_pct, abs=0.01), (
            f"{category}: expected coverage_pct={coverage_pct}, got {row['coverage_pct']}"
        )


def test_display_order_preserved(real_data_path: Path) -> None:
    """Assert categories appear in the fixed PomBase display order."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_coverage_data(real_data_path)
    assert df["category"].tolist() == ["viable", "inviable", "condition-dependent", "unknown"]


def test_dual_artifacts_created(real_data_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_coverage_data(real_data_path)
    render_coverage_figure(df, output_stem)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"
    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert an empty coverage statistics TSV is handled gracefully without crashing."""
    empty_tsv = tmp_path / "empty.tsv"
    pd.DataFrame(columns=REQUIRED_COLUMNS).to_csv(empty_tsv, sep="\t", index=False)

    df = load_coverage_data(empty_tsv)
    assert df.empty

    render_coverage_figure(df, tmp_path / "empty_test")

    pdf_path = tmp_path / "empty_test.pdf"
    assert pdf_path.exists(), "PDF artifact not created for empty data"


def test_zero_total_category_handled(tmp_path: Path) -> None:
    """Assert a category with zero covered/not_covered genes renders a placeholder, not a crash."""
    zero_tsv = tmp_path / "zero.tsv"
    pd.DataFrame(
        {
            "category": ["viable", "empty_category"],
            "covered": [10, 0],
            "not_covered": [5, 0],
            "total": [15, 0],
            "coverage_pct": [66.667, 0.0],
        }
    ).to_csv(zero_tsv, sep="\t", index=False)

    df = load_coverage_data(zero_tsv)
    render_coverage_figure(df, tmp_path / "zero_test")

    pdf_path = tmp_path / "zero_test.pdf"
    assert pdf_path.exists(), "PDF artifact not created when a category has zero genes"


def test_missing_required_column_returns_none(tmp_path: Path) -> None:
    """Assert a coverage TSV missing a required column is rejected (logger.catch swallows the raise)."""
    bad_tsv = tmp_path / "bad.tsv"
    pd.DataFrame({"category": ["viable"]}).to_csv(bad_tsv, sep="\t", index=False)

    result = load_coverage_data(bad_tsv)
    assert result is None
