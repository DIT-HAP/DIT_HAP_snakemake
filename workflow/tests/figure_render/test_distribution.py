#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Distribution of Curve Fitting Results Figure Rendering
=================================================================

Validates rendering against real data, asserting baseline statistics and artifact creation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import importlib.util
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from figure_render.histogram import render_histogram_grid_figure


def _load_distribution_script() -> ModuleType:
    """Load the distribution CLI script by path; workflow/scripts/figures has no __init__.py."""
    path = Path(__file__).resolve().parents[2] / "scripts" / "figures" / "plot_distribution_of_curve_fitting.py"
    spec = importlib.util.spec_from_file_location("_script_plot_distribution_of_curve_fitting", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCRIPT = _load_distribution_script()
load_fitting_stats = _SCRIPT.load_fitting_stats


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_stats_path() -> Path:
    """Path to real fitting statistics TSV."""
    return Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv")


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_distribution"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_statistics(real_stats_path: Path) -> None:
    """Assert baseline counts and mean statistics against real data."""
    if not real_stats_path.exists():
        pytest.skip(f"Real data not found: {real_stats_path}")

    df, metric_cols = load_fitting_stats(real_stats_path)

    # Total successful fits: 93560
    assert len(df) == 93560, f"Expected 93560 successful fits, got {len(df)}"

    # Should have 15 metric columns
    assert len(metric_cols) == 15, f"Expected 15 metric columns, got {len(metric_cols)}"

    # Check mean R2 matches baseline
    mean_r2 = df['R2'].mean()
    assert abs(mean_r2 - 0.111230) < 0.001, f"Expected mean R2 ≈ 0.111230, got {mean_r2:.6f}"

    # Check mean RMSE matches baseline
    mean_rmse = df['RMSE'].mean()
    assert abs(mean_rmse - 0.163) < 0.01, f"Expected mean RMSE ≈ 0.163, got {mean_rmse:.3f}"


def test_dual_artifacts_created(real_stats_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_stats_path.exists():
        pytest.skip(f"Real data not found: {real_stats_path}")

    df, metric_cols = load_fitting_stats(real_stats_path)
    render_histogram_grid_figure(df, output_stem, value_columns=metric_cols, bins=30)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"

    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert empty data case is handled gracefully."""
    # Create empty TSV with correct schema
    empty_stats = tmp_path / "empty_stats.tsv"
    stats_df = pd.DataFrame(columns=['Chr', 'Coordinate', 'Strand', 'Target', 'Status', 'A', 'DR', 'DL', 'R2', 'RMSE'])
    stats_df.to_csv(empty_stats, sep='\t', index=False)

    output_stem = tmp_path / "empty_test"

    # Load should return empty dataframe
    df, metric_cols = load_fitting_stats(empty_stats)
    assert df.empty

    # Render should not crash
    render_histogram_grid_figure(df, output_stem, value_columns=metric_cols, bins=30)
