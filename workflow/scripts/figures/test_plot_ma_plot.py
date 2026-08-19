#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for MA Plot Figure Rendering
===================================

Validates loading and rendering of the wide-format baseMean/LFC tables,
against both real project data and small synthetic fixtures.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-18
Version:  2.0.0
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
from plot_ma_plot import load_ma_data, render_ma_figure, Orientation, PlotConfig  # noqa: E402


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_basemean_path() -> Path:
    """Path to the real project baseMean TSV."""
    return Path(
        "/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/"
        "results/14_insertion_level_depletion_analysis/baseMean.tsv"
    )


@pytest.fixture
def real_lfc_path() -> Path:
    """Path to the real project LFC TSV."""
    return Path(
        "/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/"
        "results/14_insertion_level_depletion_analysis/LFC.tsv"
    )


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_ma_plot"


@pytest.fixture
def small_index() -> pd.MultiIndex:
    """A tiny 4-level row index shared by synthetic baseMean/LFC fixtures."""
    return pd.MultiIndex.from_tuples(
        [("chr1", 100, "+", "geneA"), ("chr1", 200, "-", "geneB"), ("chr2", 50, "+", "geneC")],
        names=["Chr", "Coordinate", "Strand", "Target"],
    )


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_statistics(real_basemean_path: Path, real_lfc_path: Path) -> None:
    """Assert baseline shape against real project data."""
    if not (real_basemean_path.exists() and real_lfc_path.exists()):
        pytest.skip(f"Real data not found: {real_basemean_path}, {real_lfc_path}")

    basemean_df, lfc_df = load_ma_data(real_basemean_path, real_lfc_path)

    assert len(basemean_df) == len(lfc_df) > 0
    assert list(lfc_df.columns) == list(basemean_df.columns)


def test_load_ma_data_rejects_mismatched_timepoints(small_index: pd.MultiIndex, tmp_path: Path) -> None:
    """Assert a timepoint present in LFC but absent from baseMean is rejected."""
    basemean_df = pd.DataFrame({"YES0": [10.0, 20.0, 30.0]}, index=small_index)
    lfc_df = pd.DataFrame({"YES0": [0.0, 0.0, 0.0], "YES1": [0.5, -0.3, 0.0]}, index=small_index)

    basemean_path = tmp_path / "baseMean.tsv"
    lfc_path = tmp_path / "LFC.tsv"
    basemean_df.to_csv(basemean_path, sep="\t")
    lfc_df.to_csv(lfc_path, sep="\t")

    assert load_ma_data(basemean_path, lfc_path) is None


@pytest.mark.parametrize("orientation", [Orientation.VERTICAL, Orientation.HORIZONTAL])
def test_dual_artifacts_created(small_index: pd.MultiIndex, output_stem: Path, orientation: Orientation) -> None:
    """Assert that both PDF and PNG artifacts are created and non-empty for synthetic data, in either orientation."""
    basemean_df = pd.DataFrame(
        {"YES0": [10.0, 20.0, 30.0], "YES1": [10.0, 20.0, 30.0]}, index=small_index
    )
    lfc_df = pd.DataFrame(
        {"YES0": [0.0, 0.0, 0.0], "YES1": [0.5, -0.3, 0.0]}, index=small_index
    )

    render_ma_figure(basemean_df, lfc_df, output_stem, orientation)

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"

    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


@pytest.mark.parametrize("orientation", [Orientation.VERTICAL, Orientation.HORIZONTAL])
def test_empty_data_handling(tmp_path: Path, orientation: Orientation) -> None:
    """Assert empty data case is handled gracefully in either orientation."""
    empty_basemean = pd.DataFrame(columns=["YES0", "YES1"])
    empty_lfc = pd.DataFrame(columns=["YES0", "YES1"])

    output_stem = tmp_path / "empty_test"

    # Render should not crash
    render_ma_figure(empty_basemean, empty_lfc, output_stem, orientation)


def test_config_rejects_missing_input(tmp_path: Path) -> None:
    """Assert PlotConfig rejects a non-existent input path."""
    with pytest.raises(ValueError, match="does not exist"):
        PlotConfig(
            basemean_path=tmp_path / "nope_baseMean.tsv",
            lfc_path=tmp_path / "nope_LFC.tsv",
            output_stem=tmp_path / "out",
        )
