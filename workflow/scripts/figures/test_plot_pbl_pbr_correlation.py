#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for PBL-PBR Correlation Figure Rendering
==============================================

Validates rendering against real data, asserting baseline statistics and artifact creation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-13
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str(SCRIPT_DIR.resolve()))
from plot_pbl_pbr_correlation import load_and_prepare_data, render_correlation_figure  # noqa: E402


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to real PBL-PBR pairs TSV (computed layer output)."""
    # Note: paths relative to main checkout, not worktree
    return Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/18_figure_data/pbl_pbr_pairs.tsv")


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_pbl_pbr_correlation"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_statistics(real_data_path: Path) -> None:
    """Assert baseline row counts and correlation for YES0 group against real data."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    # Filter to YES0 group (HD1328-4_YES0_YES.tsv)
    yes0_group = df[(df['sample'] == 'HD1328-4') & (df['timepoint'] == 'YES0') & (df['condition'] == 'YES')]

    # Assert row count matches baseline
    assert len(yes0_group) == 131057, f"Expected 131057 rows for YES0 group, got {len(yes0_group)}"

    # Compute Pearson correlation on log10 columns
    correlation = yes0_group[['log10_pbl', 'log10_pbr']].corr().iloc[0, 1]

    # Assert correlation matches baseline (tolerance 0.001)
    assert abs(correlation - 0.8521) < 0.001, f"Expected PCC ≈ 0.8521, got {correlation:.4f}"


def test_dual_artifacts_created(real_data_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)
    render_correlation_figure(df, output_stem)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"

    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_total_row_count(real_data_path: Path) -> None:
    """Assert total row count matches sum of both files (131057 + 128039 = 259096)."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_and_prepare_data(real_data_path)

    # Total should be 259096 (131057 + 128039)
    assert len(df) == 259096, f"Expected 259096 total rows, got {len(df)}"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert empty data case is handled gracefully."""
    # Create empty TSV with correct schema
    empty_tsv = tmp_path / "empty.tsv"
    empty_df = pd.DataFrame(columns=['sample', 'timepoint', 'condition', 'pbl', 'pbr'])
    empty_df.to_csv(empty_tsv, sep='\t', index=False)

    output_stem = tmp_path / "empty_test"

    # Load should return empty dataframe
    df = load_and_prepare_data(empty_tsv)
    assert df.empty

    # Render should not crash
    render_correlation_figure(df, output_stem)
