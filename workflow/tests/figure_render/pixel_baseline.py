#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Pixel baseline harness for figure render refactoring.

Renders every figure from ``HD_DIT_HAP`` data and compares PNGs byte-exactly.
PNG output was verified deterministic (max abs diff 0 across repeat renders),
so comparison needs no tolerance. PDFs embed timestamps and are not compared.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
from PIL import Image

# =============================================================================
# CONSTANTS
# =============================================================================
PROJECT_ROOT = Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake")
ARC_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/18_figure_data/arc"
DEPLETION_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/14_insertion_level_depletion_analysis"
FITTING_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/15_insertion_level_curve_fitting"

# Figures whose pixels are expected to change by design (spec section
# "Expected pixel changes"): correlation adopts orientation's grid layout,
# curve_fitting gets the house-style palette.
EXPECTED_TO_CHANGE = frozenset({"correlation", "curve_fitting"})


# =============================================================================
# CORE LOGIC
# =============================================================================
def compare_png(baseline: Path, candidate: Path) -> tuple[bool, int]:
    """Compare two PNGs pixel-exactly, returning (identical, max_abs_diff)."""
    a = np.asarray(Image.open(baseline).convert("RGB")).astype(int)
    b = np.asarray(Image.open(candidate).convert("RGB")).astype(int)

    if a.shape != b.shape:
        return False, -1

    diff = int(np.abs(a - b).max())
    return diff == 0, diff


def render_all_baseline_figures(out_dir: Path) -> dict[str, Path]:
    """Render all ten figures into out_dir, returning {figure_name: png_path}.

    Imports are function-local so this module stays importable while modules are
    mid-refactor: a figure whose module has been replaced raises ImportError,
    which is recorded rather than aborting the whole sweep.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rendered: dict[str, Path] = {}

    def _run(name: str, fn) -> None:
        stem = out_dir / name
        try:
            fn(stem)
        except ImportError as exc:
            print(f"SKIP {name}: {exc}")
            return
        png = stem.parent / f"{stem.name}.review.png"
        if png.exists():
            rendered[name] = png
        else:
            print(f"MISSING {name}: no PNG produced")

    _run("correlation", _render_correlation)
    _run("orientation", _render_orientation)
    _run("read_counts", _render_read_counts)
    _run("density", _render_density)
    _run("coverage", _render_coverage)
    _run("dispersions", _render_dispersions)
    _run("ma_plot_replicates", _render_ma_replicates)
    _run("ma_plot", _render_ma_plot)
    _run("distribution", _render_distribution)
    _run("curve_fitting", _render_curve_fitting)

    return rendered


def _render_correlation(stem: Path) -> None:
    from figure_render.correlation import load_and_prepare_data, render_correlation_figure
    df = load_and_prepare_data(ARC_DIR / "pbl_pbr_pairs.tsv")
    render_correlation_figure(df, stem)


def _render_orientation(stem: Path) -> None:
    from figure_render.orientation import load_and_prepare_data, render_orientation_figure
    df = load_and_prepare_data(ARC_DIR / "strand_pairs.tsv")
    render_orientation_figure(df, stem)


def _render_read_counts(stem: Path) -> None:
    from figure_render.read_counts import (
        load_cutoff_stats, load_distribution_data, render_distribution_figure,
    )
    df = load_distribution_data(ARC_DIR / "read_count_distribution.tsv")
    stats = load_cutoff_stats(ARC_DIR / "read_count_cutoff_stats.tsv")
    render_distribution_figure(df, stats, stem, "YES0", 8.0)


def _render_density(stem: Path) -> None:
    from figure_render.density import load_density_data, render_density_figure
    df = load_density_data(ARC_DIR / "insertion_density_analysis.tsv")
    render_density_figure(df, stem, "YES0", "YES4")


def _render_coverage(stem: Path) -> None:
    from figure_render.coverage import load_coverage_data, render_coverage_figure
    df = load_coverage_data(ARC_DIR / "gene_coverage_stats.tsv")
    render_coverage_figure(df, stem)


def _render_dispersions(stem: Path) -> None:
    from figure_render.dispersions import load_dispersion_data, render_dispersion_figure
    df = load_dispersion_data(ARC_DIR / "dispersion_data.tsv")
    render_dispersion_figure(df, stem)


def _render_ma_replicates(stem: Path) -> None:
    from figure_render.ma_plot_replicates import load_ma_data, render_ma_figure
    df = load_ma_data(ARC_DIR / "ma_values.tsv")
    render_ma_figure(df, stem)


def _render_ma_plot(stem: Path) -> None:
    from figure_render.ma_plot import Orientation, load_ma_data, render_ma_figure
    basemean, lfc = load_ma_data(DEPLETION_DIR / "baseMean.tsv", DEPLETION_DIR / "LFC.tsv")
    render_ma_figure(basemean, lfc, stem, Orientation.VERTICAL)


def _render_distribution(stem: Path) -> None:
    from figure_render.distribution import load_fitting_stats, render_distribution_figure
    df, metric_cols = load_fitting_stats(FITTING_DIR / "insertion_level_fitting_statistics.tsv")
    render_distribution_figure(df, metric_cols, stem, 30)


def _render_curve_fitting(stem: Path) -> None:
    from figure_render.curve_fitting import load_and_sample_data, render_curve_fitting_figure
    stats, lfc, time_points = load_and_sample_data(
        FITTING_DIR / "insertion_level_fitting_statistics.tsv",
        DEPLETION_DIR / "LFC.tsv",
        32,
        42,
    )
    render_curve_fitting_figure(stats, lfc, time_points, stem)
