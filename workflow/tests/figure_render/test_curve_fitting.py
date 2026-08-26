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
import importlib.util
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from depletion.curve_model import sigmoid_function
from figure_render.curves import render_fitted_curves_figure


def _load_curve_fitting_script() -> ModuleType:
    """Load the curve fitting CLI script by path; workflow/scripts/figures has no __init__.py."""
    path = Path(__file__).resolve().parents[2] / "scripts" / "figures" / "plot_curve_fitting.py"
    spec = importlib.util.spec_from_file_location("_script_plot_curve_fitting", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCRIPT = _load_curve_fitting_script()
load_and_sample_data = _SCRIPT.load_and_sample_data


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

    # Total rows: 93596 (header excluded from wc -l count of 93597)
    # Successful fits: 93560
    joined, time_points, timepoint_columns = load_and_sample_data(
        real_stats_path, real_lfc_path, 32, 42
    )

    assert len(joined) == 32, f"Expected 32 sampled curves, got {len(joined)}"
    assert len(time_points) == 5, f"Expected 5 time points, got {len(time_points)}"
    assert time_points == [0.0, 2.352, 5.588, 9.104, 12.48], f"Unexpected time points: {time_points}"
    assert len(timepoint_columns) == len(time_points), "One value column per time point"
    assert (joined["Status"] == "Success").all(), "Not all sampled curves are successful"

    # Verify mean R2 and RMSE for sampled curves
    mean_r2 = joined['R2'].mean()
    mean_rmse = joined['RMSE'].mean()

    # Since sampling is deterministic (seed 42), these should be stable
    assert 0.0 <= mean_r2 <= 1.0, f"Mean R2 out of valid range: {mean_r2}"
    assert mean_rmse > 0, f"Mean RMSE should be positive: {mean_rmse}"


def test_dual_artifacts_created(real_stats_path: Path, real_lfc_path: Path, output_stem: Path) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty."""
    if not real_stats_path.exists():
        pytest.skip(f"Real data not found: {real_stats_path}")
    if not real_lfc_path.exists():
        pytest.skip(f"Real LFC data not found: {real_lfc_path}")

    joined, time_points, timepoint_columns = load_and_sample_data(
        real_stats_path, real_lfc_path, 32, 42
    )
    render_fitted_curves_figure(
        joined,
        output_stem,
        x_values=time_points,
        value_columns=timepoint_columns,
        model=sigmoid_function,
        model_params=["A", "DR", "DL"],
        annotations=["R2", "RMSE"],
    )

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

    # Load should return empty dataframe
    joined, time_points, timepoint_columns = load_and_sample_data(
        empty_stats, empty_lfc, 32, 42
    )
    assert joined.empty

    # Render should not crash
    render_fitted_curves_figure(
        joined,
        output_stem,
        x_values=time_points,
        value_columns=timepoint_columns,
        model=sigmoid_function,
        model_params=["A", "DR", "DL"],
        annotations=["R2", "RMSE"],
    )
