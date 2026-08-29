"""Gene coverage analysis by viability: stats computation and table serialization.

Extracted from ``workflow/scripts/quality_control/gene_coverage_analysis.py``.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import pandas as pd
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
INSERTION_KEY = ["Chr", "Coordinate", "Strand", "Target"]

# Fixed display order for viability categories; any unlisted label is appended.
VIABILITY_ORDER = ["viable", "inviable", "condition-dependent", "unknown"]


class ViabilityCol(StrEnum):
    """Column names assigned to the headerless gene-viability TSV."""
    GENE_ID = "systematic_id"
    VIABILITY = "viability"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class CoverageStat:
    """Coverage tally for one viability category."""
    category: str
    total: int
    covered: int

    @property
    def not_covered(self) -> int:
        """Number of genes in this category with no insertion."""
        return self.total - self.covered

    @property
    def coverage_pct(self) -> float:
        """Percentage of genes in this category that are covered."""
        return self.covered / self.total * 100 if self.total > 0 else 0.0


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_covered_genes(lfc_file: Path, annotation_file: Path) -> set[str]:
    """Return the set of gene IDs with at least one surviving insertion."""
    lfc = pd.read_csv(lfc_file, sep="\t", usecols=INSERTION_KEY)
    logger.info(f"Loaded {len(lfc):,} insertions from LFC table")

    annotation = pd.read_csv(annotation_file, sep="\t", usecols=[*INSERTION_KEY, "Systematic ID"])
    logger.info(f"Loaded {len(annotation):,} annotation rows")

    merged = lfc.merge(annotation, on=INSERTION_KEY, how="inner")
    covered_genes = set(merged["Systematic ID"].dropna().astype(str))
    logger.success(f"Found {len(covered_genes):,} covered genes")

    return covered_genes


@logger.catch
def load_gene_viability(viability_file: Path) -> pd.DataFrame:
    """Load the headerless gene-viability TSV into a labelled DataFrame."""
    viability = pd.read_csv(
        viability_file,
        sep="\t",
        header=None,
        names=[ViabilityCol.GENE_ID, ViabilityCol.VIABILITY],
    )
    viability[ViabilityCol.GENE_ID] = viability[ViabilityCol.GENE_ID].astype(str)
    logger.info(f"Loaded viability for {len(viability):,} genes")

    return viability


@logger.catch
def compute_coverage_stats(viability: pd.DataFrame, covered_genes: set[str]) -> list[CoverageStat]:
    """Tally covered vs total genes per viability category in display order."""
    viability = viability.copy()
    viability["is_covered"] = viability[ViabilityCol.GENE_ID].isin(covered_genes)

    present = viability[ViabilityCol.VIABILITY].unique().tolist()
    ordered = [c for c in VIABILITY_ORDER if c in present]
    ordered += [c for c in present if c not in VIABILITY_ORDER]

    stats: list[CoverageStat] = []
    for category in ordered:
        group = viability[viability[ViabilityCol.VIABILITY] == category]
        stat = CoverageStat(category=category, total=len(group), covered=int(group["is_covered"].sum()))
        logger.info(f"{category}: {stat.covered:,}/{stat.total:,} covered ({stat.coverage_pct:.1f}%)")
        stats.append(stat)

    return stats
