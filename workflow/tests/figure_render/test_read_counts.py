#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Read Count Distribution Figure (single-stage)
=======================================================

The pipeline no longer produces pre-binned TSVs; the figure script reads raw
read-count tables and histograms them in memory. Tests therefore drive the
script's assembly function against the real HD_DIT_HAP merged tables (baseline
numbers below were recomputed bit-exactly from the same tables at refactor
time) and the new grouped log-scale renderer on synthetic frames.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-27
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import importlib.util
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from figure_render.histogram import render_grouped_histogram_figure


def _load_read_counts_script() -> ModuleType:
    """Load the read counts CLI script by path; workflow/scripts/figures has no __init__.py."""
    path = Path(__file__).resolve().parents[2] / "scripts" / "figures" / "plot_read_count_distribution.py"
    spec = importlib.util.spec_from_file_location("_script_plot_read_count_distribution", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCRIPT = _load_read_counts_script()
assemble_distribution = _SCRIPT.assemble_distribution
X_LABEL = _SCRIPT.X_LABEL
Y_LABEL = _SCRIPT.Y_LABEL

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
PROJECT_ROOT = Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake")
MERGED_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/11_merged"

# Baseline retention per sample, recomputed bit-exactly from the same raw
# merged tables the single-stage script now consumes directly.
BASELINE_RETENTION = {
    "HD1328-4_YES": (659362, 143652, 21.7865, 380888542.0, 99.9404),
    "HD1328-7_YES": (287220, 123271, 42.9187, 282495932.0, 99.9809),
    "HD1328-8_YES": (287186, 129152, 44.9716, 416707741.0, 99.9843),
}

# Positive-value counts per sample/timepoint from the archived binned era; the
# script must see exactly these many strictly-positive values per group.
BASELINE_POSITIVE_COUNTS = {
    ("HD1328-4_YES", "YES0"): 293306,
    ("HD1328-4_YES", "YES1"): 212152,
    ("HD1328-4_YES", "YES2"): 215304,
    ("HD1328-4_YES", "YES3"): 211483,
    ("HD1328-4_YES", "YES4"): 220481,
}

TIMEPOINTS = ["YES0", "YES1", "YES2", "YES3", "YES4"]


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def real_input_paths() -> list[Path]:
    """Paths to the real merged read-count tables."""
    paths = sorted(MERGED_DIR.glob("*.merged.tsv"))
    if not paths:
        pytest.skip(f"No merged tables found in {MERGED_DIR}")
    return paths


@pytest.fixture
def output_stem(tmp_path: Path) -> Path:
    """Temporary output stem for test artifacts."""
    return tmp_path / "test_read_count_distribution"


# =============================================================================
# TESTS
# =============================================================================
def test_baseline_retention_statistics(real_input_paths: list[Path]) -> None:
    """Assert cutoff retention recomputed from raw tables matches the archive bit-for-bit."""
    _, retention_lines = assemble_distribution(real_input_paths, "YES0", 8.0)

    assert len(retention_lines) == len(BASELINE_RETENTION), (
        f"Expected {len(BASELINE_RETENTION)} samples, got {len(retention_lines)}"
    )

    by_sample = {line.split(":")[0]: line for line in retention_lines}
    for sample, (orig_rows, rows_kept, pct_rows, _, pct_counts) in BASELINE_RETENTION.items():
        line = by_sample[sample]
        # The caption embeds both row and count retention percentages.
        assert f"{rows_kept:,}/{orig_rows:,} rows kept ({pct_rows:.1f}%)" in line, f"{sample}: row retention"
        assert f"({pct_counts:.1f}%)" in line, f"{sample}: count retention"


def test_baseline_positive_value_counts(real_input_paths: list[Path]) -> None:
    """Assert each group holds exactly as many positive values as the archived bins summed to."""
    df, _ = assemble_distribution(real_input_paths, "YES0", 8.0)
    sums = df[df["value"] > 0].groupby(["sample", "timepoint"]).size()

    for key, expected in BASELINE_POSITIVE_COUNTS.items():
        assert sums[key] == expected, f"{key}: expected {expected} positive values, got {sums[key]}"


def test_long_format_group_structure(real_input_paths: list[Path]) -> None:
    """Assert the assembled frame is long format with one timepoint column per group."""
    df, _ = assemble_distribution(real_input_paths, "YES0", 8.0)

    assert set(df.columns) == {"sample", "timepoint", "value"}
    assert sorted(df["timepoint"].unique()) == TIMEPOINTS
    assert len(df["sample"].unique()) == 3


def test_dual_artifacts_created(real_input_paths: list[Path], output_stem: Path) -> None:
    """Assert end-to-end rendering creates non-empty PDF and PNG artifacts."""
    df, retention_lines = assemble_distribution(real_input_paths, "YES0", 8.0)

    render_grouped_histogram_figure(
        df,
        output_stem,
        value_column="value",
        row_key="sample",
        col_key="timepoint",
        bins=50,
        log_scale=True,
        xlabel=X_LABEL,
        ylabel=Y_LABEL,
        marker_value=8.0,
        marker_label="Cutoff = 8",
        marker_on_col_value="YES0",
    )

    pdf_path = output_stem.parent / f"{output_stem.name}.pdf"
    png_path = output_stem.parent / f"{output_stem.name}.review.png"

    assert pdf_path.exists(), f"PDF artifact not created: {pdf_path}"
    assert png_path.exists(), f"PNG artifact not created: {png_path}"
    assert pdf_path.stat().st_size > 0, "PDF artifact is empty"
    assert png_path.stat().st_size > 0, "PNG artifact is empty"


def test_log_scale_logspace_binning(tmp_path: Path) -> None:
    """Assert log-scale mode sets a true log axis and partitions decades into distinct bins.

    The claim of the beta redesign: np.logspace edges on a log axis are
    geometrically identical to equal-width bins over log10(values). Three
    values one decade apart with 3 bins must each occupy its own bin.
    """
    df = pd.DataFrame({
        "sample": ["S1"] * 3,
        "timepoint": ["T0"] * 3,
        "value": [1.0, 100.0, 10000.0],
    })

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    render_grouped_histogram_figure(
        df,
        tmp_path / "logscale",
        value_column="value",
        row_key="sample",
        col_key="timepoint",
        bins=3,
        log_scale=True,
    )
    fig = plt.gcf()
    ax = fig.axes[0]
    assert ax.get_xscale() == "log", "Log-scale mode must set a true log x-axis"
    heights = [patch.get_height() for patch in ax.patches]
    assert sum(1 for h in heights if h > 0) == 3, "Each of 3 distinct decades should occupy its own bin"
    plt.close(fig)


def test_no_valid_data_panel(tmp_path: Path) -> None:
    """Assert a group whose values are all non-positive renders without crashing."""
    marker_df = pd.DataFrame({
        "sample": ["S1", "S1", "S1"],
        "timepoint": ["YES0", "YES1", "YES1"],
        "value": [10.0, -1.0, -5.0],
    })

    render_grouped_histogram_figure(
        marker_df,
        tmp_path / "marker_test",
        value_column="value",
        row_key="sample",
        col_key="timepoint",
        log_scale=True,
    )

    assert (tmp_path / "marker_test.pdf").exists(), "PDF not created for no-valid-data input"


def test_empty_data_handling(tmp_path: Path) -> None:
    """Assert an empty frame is handled gracefully without crashing."""
    empty_df = pd.DataFrame(columns=["sample", "timepoint", "value"])

    render_grouped_histogram_figure(
        empty_df,
        tmp_path / "empty_test",
        value_column="value",
        row_key="sample",
        col_key="timepoint",
        log_scale=True,
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
