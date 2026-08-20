#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Replicate-Branch MA Plot Figure Rendering
===================================================

Validates rendering against real data, asserting baseline statistics, NaN-padj
colour semantics, and artifact creation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-16
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pandas as pd
import pytest

from figure_render.ma_plot_replicates import (
    NONSIGNIFICANT_COLOR,
    PADJ_THRESHOLD,
    REQUIRED_COLUMNS,
    SIGNIFICANT_COLOR,
    load_ma_data,
    render_ma_figure,
    significance_colors,
)


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to real replicate-branch MA values TSV."""
    return Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/18_figure_data/ma_values_replicates.tsv")


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_ma_plot_replicates"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_statistics(real_data_path: Path) -> None:
    """Assert baseline row counts and schema against real data."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_ma_data(real_data_path)

    # HD_DIT_HAP: 93596 insertions x 4 non-initial timepoints (YES1-YES4)
    assert len(df) == 374_384, f"Expected 374384 rows, got {len(df)}"

    timepoints = sorted(df['timepoint'].unique())
    assert len(timepoints) == 4, f"Expected 4 timepoints, got {timepoints}"

    for col in ("baseMean", "log2FoldChange", "padj"):
        assert df[col].dtype in ['float64', 'float32'], f"{col} should be float, got {df[col].dtype}"

    # baseMean feeds a log x-axis, so it must be strictly positive
    assert (df['baseMean'] > 0).all(), "baseMean must be positive for log-scale x-axis"

    # Independent filtering leaves NaN padj values that must survive loading
    assert df['padj'].isna().any(), "Expected some NaN padj values from independent filtering"


def test_nan_padj_renders_nonsignificant() -> None:
    """Assert NaN padj maps to the non-significant colour, as in pydeseq2."""
    padj = pd.Series([0.001, PADJ_THRESHOLD, 0.9, float('nan')])
    colors = significance_colors(padj)

    assert colors.tolist() == [
        SIGNIFICANT_COLOR,
        NONSIGNIFICANT_COLOR,  # padj == threshold is not < threshold
        NONSIGNIFICANT_COLOR,
        NONSIGNIFICANT_COLOR,  # NaN < threshold is False in pydeseq2 too
    ]
    assert colors.notna().all(), "Colour mapping must not leave NaN entries"


def test_dual_artifacts_created(real_data_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_ma_data(real_data_path)
    render_ma_figure(df, output_stem)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"

    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_missing_column_rejected(tmp_path: Path) -> None:
    """Assert a TSV missing padj is refused, not silently plotted as all-gray."""
    bad_tsv = tmp_path / "bad.tsv"
    bad_df = pd.DataFrame({"timepoint": ["YES1"], "baseMean": [10.0], "log2FoldChange": [0.5]})
    bad_df.to_csv(bad_tsv, sep='\t', index=False)

    # The loader raises ValueError, which @logger.catch logs and converts into a
    # None return, so the caller never receives a frame lacking significance data.
    assert load_ma_data(bad_tsv) is None


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert empty data case is handled gracefully."""
    empty_tsv = tmp_path / "empty.tsv"
    empty_df = pd.DataFrame(columns=list(REQUIRED_COLUMNS))
    empty_df.to_csv(empty_tsv, sep='\t', index=False)

    output_stem = tmp_path / "empty_test"

    df = load_ma_data(empty_tsv)
    assert df.empty

    # Render should not crash
    render_ma_figure(df, output_stem)
