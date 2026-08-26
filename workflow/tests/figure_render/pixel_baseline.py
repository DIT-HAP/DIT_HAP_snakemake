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


def _load_script(stem: str):
    """Load a figure CLI script by stem; workflow/scripts/figures has no __init__.py."""
    import importlib.util

    path = PROJECT_ROOT / "workflow" / "scripts" / "figures" / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"_script_{stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _render_correlation(stem: Path) -> None:
    from figure_render.scatter import render_grouped_regression_figure
    script = _load_script("plot_pbl_pbr_correlation")
    df = script.load_and_prepare_data(ARC_DIR / "pbl_pbr_pairs.tsv")
    render_grouped_regression_figure(
        df, stem, x=script.X_COLUMN, y=script.Y_COLUMN,
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL,
        row_key="sample", col_key="timepoint",
    )


def _render_orientation(stem: Path) -> None:
    from figure_render.scatter import render_grouped_regression_figure
    script = _load_script("plot_insertion_orientation")
    df = script.load_and_prepare_data(ARC_DIR / "strand_pairs.tsv")
    render_grouped_regression_figure(
        df, stem, x=script.X_COLUMN, y=script.Y_COLUMN,
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL,
        row_key="sample", col_key="timepoint",
    )


def _render_read_counts(stem: Path) -> None:
    from figure_render.histogram import render_prebinned_histogram_figure
    script = _load_script("plot_read_count_distribution")
    df = script.load_distribution_data(ARC_DIR / "read_count_distribution.tsv")
    stats = script.load_cutoff_stats(ARC_DIR / "read_count_cutoff_stats.tsv")
    initial_time_point = "YES0"
    cutoff = 8.0
    render_prebinned_histogram_figure(
        df, stem,
        row_key="sample", col_key="timepoint",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL,
        marker_value=float(np.log10(cutoff)),
        marker_label=f"Cutoff = {cutoff:.2g}",
        marker_on_col_value=initial_time_point,
        footer_lines=script.build_retention_footer(df, stats),
        footer_header=f"Cutoff applied to '{initial_time_point}' (>= {cutoff:.2g}):",
    )


def _render_density(stem: Path) -> None:
    from figure_render.scatter import render_scatter_grid_figure
    script = _load_script("plot_insertion_density")
    df = script.load_density_data(ARC_DIR / "insertion_density_analysis.tsv")
    render_scatter_grid_figure(
        df, stem,
        panels=script.build_density_panels("YES0", "YES4"),
        hue=script.VIABILITY_COLUMN,
        hue_order=script.VIABILITY_HUE_ORDER,
    )


def _render_coverage(stem: Path) -> None:
    from figure_render.composition import render_composition_figure
    script = _load_script("plot_gene_coverage")
    df = script.load_coverage_data(ARC_DIR / "gene_coverage_stats.tsv")
    render_composition_figure(
        df, stem,
        category_column="category", percentage_column="coverage_pct",
        part_column="covered", whole_column="not_covered", total_column="total",
        part_label=script.COVERED_LABEL, whole_label=script.NOT_COVERED_LABEL,
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL, title=script.TITLE,
        donut_unit=script.DONUT_UNIT,
    )


def _render_dispersions(stem: Path) -> None:
    from figure_render.series import render_series_scatter_figure
    script = _load_script("plot_dispersions")
    df = script.load_dispersion_data(ARC_DIR / "dispersion_data.tsv")
    render_series_scatter_figure(
        df, stem,
        x=script.X_COLUMN, series=script.DISPERSION_SERIES,
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL, title=script.TITLE,
    )


def _render_ma_replicates(stem: Path) -> None:
    from figure_render.ma import Orientation, render_ma_figure
    script = _load_script("plot_ma_plot_replicates")
    df = script.load_ma_data(ARC_DIR / "ma_values.tsv")
    panels, colors = script.build_ma_panels(df)
    render_ma_figure(
        panels, stem,
        abundance_label=script.ABUNDANCE_LABEL, effect_label=script.EFFECT_LABEL,
        title_prefix=script.TITLE_PREFIX,
        orientation=Orientation.HORIZONTAL, stack=True,
        point_colors=colors,
        panel_width=script.PANEL_WIDTH, panel_height=script.PANEL_HEIGHT,
        share_axes=False,
    )


def _render_ma_plot(stem: Path) -> None:
    from figure_render.ma import Orientation, render_ma_figure
    script = _load_script("plot_ma_plot")
    basemean_df, lfc_df = script.load_ma_data(
        DEPLETION_DIR / "baseMean.tsv", DEPLETION_DIR / "LFC.tsv",
    )
    panels = script.build_ma_panels(basemean_df, lfc_df)
    render_ma_figure(
        panels, stem,
        abundance_label=script.ABUNDANCE_LABEL, effect_label=script.EFFECT_LABEL,
        title_prefix=script.TITLE_PREFIX, orientation=Orientation.VERTICAL,
    )


def _render_distribution(stem: Path) -> None:
    from figure_render.histogram import render_histogram_grid_figure
    script = _load_script("plot_distribution_of_curve_fitting")
    df, metric_cols = script.load_fitting_stats(
        FITTING_DIR / "insertion_level_fitting_statistics.tsv",
    )
    render_histogram_grid_figure(df, stem, value_columns=metric_cols, bins=30)


def _render_curve_fitting(stem: Path) -> None:
    from depletion.curve_model import sigmoid_function
    from figure_render.curves import render_fitted_curves_figure
    script = _load_script("plot_curve_fitting")
    joined, time_points, timepoint_columns = script.load_and_sample_data(
        FITTING_DIR / "insertion_level_fitting_statistics.tsv",
        DEPLETION_DIR / "LFC.tsv",
        32,
        42,
    )
    render_fitted_curves_figure(
        joined, stem,
        x_values=time_points,
        value_columns=timepoint_columns,
        model=sigmoid_function,
        model_params=script.MODEL_PARAM_COLUMNS,
        annotations=script.ANNOTATION_COLUMNS,
        xlabel=script.X_LABEL,
        ylabel=script.Y_LABEL,
    )
