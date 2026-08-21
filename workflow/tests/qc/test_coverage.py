"""Tests for gene coverage tallying by viability category."""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
import pytest

from qc.coverage import ViabilityCol, compute_coverage_stats


# =============================================================================
# TESTS: compute_coverage_stats
# =============================================================================
def test_coverage_stats_two_of_three_covered():
    """2 covered of 3 viable genes gives covered=2, not_covered=1, total=3, coverage_pct=200/3.

    gene1 and gene2 are in the covered_genes set, gene3 is not. All three
    rows share the "viable" category, so the single CoverageStat reports
    total=3, covered=2. ``not_covered`` is a derived property (total-covered
    = 1) and ``coverage_pct`` is covered/total*100 = 2/3*100 = 66.666...%.
    """
    viability = pd.DataFrame({
        ViabilityCol.GENE_ID: ["gene1", "gene2", "gene3"],
        ViabilityCol.VIABILITY: ["viable", "viable", "viable"],
    })
    covered_genes = {"gene1", "gene2"}

    stats = compute_coverage_stats(viability, covered_genes)

    assert len(stats) == 1
    stat = stats[0]
    assert stat.category == "viable"
    assert stat.total == 3
    assert stat.covered == 2
    assert stat.not_covered == 1
    assert stat.coverage_pct == pytest.approx(200 / 3)


def test_coverage_stats_multiple_categories_preserve_display_order():
    """Categories are emitted in the fixed VIABILITY_ORDER, not input order.

    Input rows are given in the order inviable, viable, unknown, but
    VIABILITY_ORDER is ["viable", "inviable", "condition-dependent",
    "unknown"], so the output must reorder them to viable, inviable, unknown.
    """
    viability = pd.DataFrame({
        ViabilityCol.GENE_ID: ["geneA", "geneB", "geneC"],
        ViabilityCol.VIABILITY: ["inviable", "viable", "unknown"],
    })
    covered_genes = {"geneA", "geneB", "geneC"}

    stats = compute_coverage_stats(viability, covered_genes)

    assert [stat.category for stat in stats] == ["viable", "inviable", "unknown"]
    for stat in stats:
        assert stat.total == 1
        assert stat.covered == 1
        assert stat.coverage_pct == pytest.approx(100.0)


def test_coverage_stats_unlisted_category_is_appended():
    """A viability label absent from VIABILITY_ORDER is appended after the known ones."""
    viability = pd.DataFrame({
        ViabilityCol.GENE_ID: ["geneA", "geneB"],
        ViabilityCol.VIABILITY: ["exotic-category", "viable"],
    })
    covered_genes: set[str] = set()

    stats = compute_coverage_stats(viability, covered_genes)

    assert [stat.category for stat in stats] == ["viable", "exotic-category"]


def test_coverage_stats_zero_covered_gives_zero_pct():
    """No gene covered in a category yields covered=0 and coverage_pct=0.0, not a division error."""
    viability = pd.DataFrame({
        ViabilityCol.GENE_ID: ["gene1", "gene2"],
        ViabilityCol.VIABILITY: ["viable", "viable"],
    })
    covered_genes: set[str] = set()

    stats = compute_coverage_stats(viability, covered_genes)

    stat = stats[0]
    assert stat.covered == 0
    assert stat.not_covered == 2
    assert stat.coverage_pct == pytest.approx(0.0)


def test_coverage_stat_coverage_pct_zero_total_is_zero_not_nan():
    """CoverageStat.coverage_pct guards total==0 explicitly, returning 0.0 rather than raising."""
    from qc.coverage import CoverageStat

    stat = CoverageStat(category="empty", total=0, covered=0)
    assert stat.coverage_pct == 0.0
    assert stat.not_covered == 0
