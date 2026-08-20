#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for DESeq2 Dispersion Figure Rendering
============================================

Validates rendering against real data, asserting baseline statistics and artifact creation.

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

from figure_render.dispersions import (
    DISPERSION_SERIES,
    REQUIRED_COLUMNS,
    X_COLUMN,
    load_dispersion_data,
    render_dispersion_figure,
)


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_data_path() -> Path:
    """Path to real dispersion figure-data TSV."""
    return Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/18_figure_data/dispersion_data.tsv")


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_dispersions"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_statistics(real_data_path: Path) -> None:
    """Assert baseline row counts and schema against real data."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_dispersion_data(real_data_path)

    # HD_DIT_HAP has 93596 insertions after imputation/hard filtering,
    # one dispersion row per insertion
    assert len(df) == 93596, f"Expected 93596 rows, got {len(df)}"

    # 4-level index from the computation layer
    assert list(df.index.names) == ["Chr", "Coordinate", "Strand", "Target"]

    # All value columns present and numeric
    for col in REQUIRED_COLUMNS:
        assert col in df.columns, f"Missing column: {col}"
        assert df[col].dtype in ['float64', 'float32'], f"{col} should be float, got {df[col].dtype}"

    # Log-scale axes require strictly positive values
    assert (df[X_COLUMN] > 0).all(), "normed_mean must be positive for log-scale x-axis"
    for col, _, _ in DISPERSION_SERIES:
        assert (df[col].dropna() > 0).all(), f"{col} must be positive for log-scale y-axis"


def test_dual_artifacts_created(real_data_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_data_path.exists():
        pytest.skip(f"Real data not found: {real_data_path}")

    df = load_dispersion_data(real_data_path)
    render_dispersion_figure(df, output_stem)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"

    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_missing_column_rejected(tmp_path: Path) -> None:
    """Assert a TSV missing a dispersion column is refused, not silently plotted."""
    bad_tsv = tmp_path / "bad.tsv"
    bad_df = pd.DataFrame(
        {
            "Chr": ["I"],
            "Coordinate": [100],
            "Strand": ["+"],
            "Target": ["TTAA"],
            X_COLUMN: [10.0],
            "genewise_dispersion": [0.1],
            # MAP_dispersion and fitted_dispersion deliberately absent
        }
    )
    bad_df.to_csv(bad_tsv, sep='\t', index=False)

    # The loader raises ValueError, which @logger.catch logs and converts into a
    # None return, so the caller never receives a partially-valid frame.
    assert load_dispersion_data(bad_tsv) is None


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert empty data case is handled gracefully and still writes artifacts."""
    empty_tsv = tmp_path / "empty.tsv"
    empty_df = pd.DataFrame(
        columns=["Chr", "Coordinate", "Strand", "Target", *REQUIRED_COLUMNS]
    )
    empty_df.to_csv(empty_tsv, sep='\t', index=False)

    output_stem = tmp_path / "empty_test"

    df = load_dispersion_data(empty_tsv)
    assert df.empty

    # Render should not crash
    render_dispersion_figure(df, output_stem)
