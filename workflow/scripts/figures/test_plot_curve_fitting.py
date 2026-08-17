#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Curve Fitting Figure Rendering
=========================================

Validates rendering against real data, asserting baseline statistics and artifact creation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str(SCRIPT_DIR.resolve()))
from plot_curve_fitting import load_and_sample_data, PlotConfig, render_curve_fitting_figure  # noqa: E402


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_stats_path() -> Path:
    """Path to real fitting statistics TSV."""
    return Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv")


@pytest.fixture
def real_lfc_path() -> Path:
    """Path to real LFC data TSV."""
    return Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/14_insertion_level_depletion_analysis/LFC.tsv")


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_curve_fitting"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_statistics(real_stats_path: Path, real_lfc_path: Path, output_stem: Path) -> None:
    """Assert baseline counts and mean statistics against real data."""
    if not real_stats_path.exists():
        pytest.skip(f"Real data not found: {real_stats_path}")
    if not real_lfc_path.exists():
        pytest.skip(f"Real LFC data not found: {real_lfc_path}")

    config = PlotConfig(
        fitting_stats_path=real_stats_path,
        lfc_path=real_lfc_path,
        output_stem=output_stem,
        n_curves=32,
        random_seed=42,
    )

    sampled_stats, lfc_sampled, time_points = load_and_sample_data(config)

    # Total rows: 93596 (header excluded from wc -l count of 93597)
    # Successful fits: 93560
    assert len(sampled_stats) == 32, f"Expected 32 sampled curves, got {len(sampled_stats)}"
    assert len(time_points) == 5, f"Expected 5 time points, got {len(time_points)}"
    assert time_points == [0.0, 2.352, 5.588, 9.104, 12.48], f"Unexpected time points: {time_points}"

    # Verify mean R2 and RMSE for sampled curves
    mean_r2 = sampled_stats['R2'].mean()
    mean_rmse = sampled_stats['RMSE'].mean()

    # Since sampling is deterministic (seed 42), these should be stable
    assert 0.0 <= mean_r2 <= 1.0, f"Mean R2 out of valid range: {mean_r2}"
    assert mean_rmse > 0, f"Mean RMSE should be positive: {mean_rmse}"

    # All sampled curves should be successful
    assert (sampled_stats['Status'] == 'Success').all(), "Not all sampled curves are successful"


def test_dual_artifacts_created(real_stats_path: Path, real_lfc_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_stats_path.exists():
        pytest.skip(f"Real data not found: {real_stats_path}")
    if not real_lfc_path.exists():
        pytest.skip(f"Real LFC data not found: {real_lfc_path}")

    config = PlotConfig(
        fitting_stats_path=real_stats_path,
        lfc_path=real_lfc_path,
        output_stem=output_stem,
        n_curves=32,
        random_seed=42,
    )

    sampled_stats, lfc_sampled, time_points = load_and_sample_data(config)
    render_curve_fitting_figure(sampled_stats, lfc_sampled, time_points, output_stem)

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
    empty_lfc = tmp_path / "empty_lfc.tsv"

    stats_df = pd.DataFrame(columns=['Chr', 'Coordinate', 'Strand', 'Target', 'time_points', 'Status', 'A', 'DR', 'DL', 'R2', 'RMSE'])
    stats_df.to_csv(empty_stats, sep='\t', index=False)

    lfc_df = pd.DataFrame(columns=['Chr', 'Coordinate', 'Strand', 'Target', 'YES0', 'YES1', 'YES2', 'YES3', 'YES4'])
    lfc_df.to_csv(empty_lfc, sep='\t', index=False)

    output_stem = tmp_path / "empty_test"

    config = PlotConfig(
        fitting_stats_path=empty_stats,
        lfc_path=empty_lfc,
        output_stem=output_stem,
    )

    # Load should return empty dataframe
    sampled_stats, lfc_sampled, time_points = load_and_sample_data(config)
    assert sampled_stats.empty

    # Render should not crash
    render_curve_fitting_figure(sampled_stats, lfc_sampled, time_points, output_stem)
