"""
Tests for Insertion-Level Depletion Analysis (Replicated Samples)
===================================================================

Unit tests for the dispersion-data TSV writer that feeds the cnsplots
rendering layer. Fitting a real ``DeseqDataSet`` on the full HD_DIT_HAP
counts matrix takes minutes, so these tests use a small in-memory stand-in
that exposes only the attributes the writer actually reads (``var.index``,
``varm``); a full real-data run is verified manually instead (see the
migration task notes, not part of this automated suite).

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-16
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

# The module under test imports pydeseq2 at module scope, which only the
# computation env provides. Skip cleanly elsewhere so collection errors do not
# abort the whole suite in the rendering env.
pytest.importorskip("pydeseq2", reason="requires the pydeseq2 computation env")

from depletion.insertion_level_replicates import write_dispersion_data_tsv


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
    assert df.loc[("chr1", 100, "+", "geneA"), "normed_mean"] == pytest.approx(10.0)
    assert df.loc[("chr2", 50, "+", "geneC"), "MAP_dispersion"] == pytest.approx(0.35)
