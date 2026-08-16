#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Insertion-Level Depletion Analysis (Replicated Samples)
===================================================================

Unit tests for the dispersion-data and MA-values TSV writers that feed the
cnsplots rendering layer. Fitting a real ``DeseqDataSet`` on the full HD_DIT_HAP
counts matrix takes minutes, so these tests use small in-memory stand-ins that
expose only the attributes the writers actually read (``var.index``, ``varm``,
``results_df``); a full real-data run is verified manually instead (see the
migration task notes, not part of this automated suite).

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-16
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str(SCRIPT_DIR.resolve()))
from insertion_level_depletion_analysis_has_replicates import (  # noqa: E402
    write_dispersion_data_tsv,
    write_ma_values_tsv,
)


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def fake_dds() -> SimpleNamespace:
    """Minimal stand-in for a fitted DeseqDataSet exposing var.index and varm."""
    index = pd.Index(["chr1=100=+=geneA", "chr1=200=-=geneB", "chr2=50=+=geneC"])
    varm = {
        "_normed_means": np.array([10.0, 20.0, 30.0]),
        "genewise_dispersions": np.array([0.1, 0.2, 0.3]),
        "dispersions": np.array([0.15, 0.25, 0.35]),
        "fitted_dispersions": np.array([0.12, 0.22, 0.32]),
    }
    return SimpleNamespace(var=SimpleNamespace(index=index), varm=varm)


@pytest.fixture
def fake_stat_res() -> dict[str, SimpleNamespace]:
    """Minimal stand-in for {timepoint: DeseqStats} exposing results_df."""
    index = pd.Index(["chr1=100=+=geneA", "chr1=200=-=geneB", "chr2=50=+=geneC"])
    results = {}
    for tp, offset in [("YES1", 0.0), ("YES2", 1.0)]:
        results[tp] = SimpleNamespace(
            results_df=pd.DataFrame(
                {
                    "baseMean": [10.0, 20.0, 30.0],
                    "log2FoldChange": [0.5 + offset, -0.3 + offset, 0.0 + offset],
                    "padj": [0.01, 0.5, 0.99],
                },
                index=index,
            )
        )
    return results


# =============================================================================
# TESTS
# =============================================================================
def test_write_dispersion_data_tsv_builds_multiindex_and_columns(fake_dds: SimpleNamespace, tmp_path: Path) -> None:
    """Assert the dispersion TSV has the expected 4-level index and dispersion columns."""
    output_path = tmp_path / "dispersion_data.tsv"
    write_dispersion_data_tsv(fake_dds, output_path)

    df = pd.read_csv(output_path, sep="\t", index_col=[0, 1, 2, 3])

    assert list(df.index.names) == ["Chr", "Coordinate", "Strand", "Target"]
    assert list(df.columns) == ["normed_mean", "genewise_dispersion", "MAP_dispersion", "fitted_dispersion"]
    assert len(df) == 3
    assert df.loc[("chr1", 100, "+", "geneA"), "normed_mean"] == 10.0
    assert df.loc[("chr2", 50, "+", "geneC"), "MAP_dispersion"] == 0.35


def test_write_ma_values_tsv_concatenates_timepoints_without_negation(
    fake_stat_res: dict[str, SimpleNamespace], tmp_path: Path
) -> None:
    """Assert the MA-values TSV stacks all timepoints and preserves pre-negation sign."""
    output_path = tmp_path / "ma_values_replicates.tsv"
    write_ma_values_tsv(fake_stat_res, output_path)

    df = pd.read_csv(output_path, sep="\t")

    assert list(df.columns) == ["timepoint", "baseMean", "log2FoldChange", "padj"]
    assert len(df) == 6  # 3 insertions x 2 timepoints
    assert set(df["timepoint"].unique()) == {"YES1", "YES2"}

    # Pre-negation fidelity: YES1's first row keeps the raw +0.5 sign, not negated.
    yes1_first = df[df["timepoint"] == "YES1"].iloc[0]
    assert yes1_first["log2FoldChange"] == pytest.approx(0.5)


def test_write_ma_values_tsv_excludes_initial_timepoint(fake_stat_res: dict[str, SimpleNamespace], tmp_path: Path) -> None:
    """Assert only timepoints present in stat_res appear (no initial-timepoint row)."""
    output_path = tmp_path / "ma_values_replicates.tsv"
    write_ma_values_tsv(fake_stat_res, output_path)

    df = pd.read_csv(output_path, sep="\t")
    assert "YES0" not in set(df["timepoint"].unique())
