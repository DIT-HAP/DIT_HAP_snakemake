#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "loguru",
# ]
# ///

"""
Gene Coverage Analysis by Viability
====================================

Assess how thoroughly the transposon insertion library covers protein-coding
genes, broken down by gene viability (essentiality). A gene is "covered" when
at least one insertion maps to it in the insertion-level results.

The biological expectation is that inviable (essential) genes are covered at a
lower rate than viable genes: insertions that disrupt essential genes are lost
during outgrowth, so surviving insertions are depleted from those loci. A clear
gap between the viability groups is therefore a sanity check on both library
quality and the depletion signal.

Covered genes are derived by joining the insertion-level LFC table (which
insertions survived) to the per-insertion annotation table (which gene each
insertion sits in) on the (Chr, Coordinate, Strand, Target) key. Each gene is
labelled with its PomBase viability category, and coverage is tallied per
category.

This is the computation layer: it writes only the coverage statistics table.
Rendering (bar chart + donut charts) lives in
``workflow/scripts/figures/plot_gene_coverage.py``.

Input
-----
- ``-i`` insertion-level LFC TSV, row MultiIndex (Chr, Coordinate, Strand, Target).
- ``-a`` per-insertion annotation TSV with a ``Systematic ID`` column and the
  same 4-column insertion key.
- ``-v`` PomBase gene viability TSV, two columns, no header:
  ``systematic_id`` and ``viability`` (viable / inviable / condition-dependent / unknown).

Output
------
- ``-o`` a TSV with one row per viability category, in PomBase display order
  (viable, inviable, condition-dependent, unknown, then any unlisted labels).
  Columns: category, covered, not_covered, total, coverage_pct.

Usage
-----
    python gene_coverage_analysis.py -i LFC.tsv -a annotations.tsv -v gene_viability.tsv -o gene_coverage_stats.tsv

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402

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
class InputOutputConfig:
    """Input/output paths for the gene coverage analysis."""
    lfc_file: Path
    annotation_file: Path
    viability_file: Path
    output_file: Path

    def __post_init__(self) -> None:
        """Validate inputs exist and ensure the output directory is present."""
        for path in (self.lfc_file, self.annotation_file, self.viability_file):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


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


@logger.catch
def write_coverage_table(stats: list[CoverageStat], output_file: Path) -> None:
    """Write the per-category coverage tally to a TSV."""
    table = pd.DataFrame(
        [
            {
                "category": stat.category,
                "covered": stat.covered,
                "not_covered": stat.not_covered,
                "total": stat.total,
                "coverage_pct": round(stat.coverage_pct, 3),
            }
            for stat in stats
        ]
    )
    table.to_csv(output_file, sep="\t", index=False)
    logger.success(f"Wrote coverage statistics for {len(stats)} categories to {output_file}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Analyze gene coverage by viability from insertion-level results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gene_coverage_analysis.py -i LFC.tsv -a annotations.tsv -v gene_viability.tsv -o gene_coverage_stats.tsv
        """,
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Insertion-level LFC TSV")
    parser.add_argument("-a", "--annotation", type=Path, required=True,
                        help="Per-insertion annotation TSV (with Systematic ID column)")
    parser.add_argument("-v", "--gene-viability", type=Path, required=True,
                        help="PomBase gene viability TSV (headerless: id, viability)")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output coverage statistics TSV path")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: load coverage, tally by viability, write the statistics TSV."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = InputOutputConfig(
            lfc_file=args.input,
            annotation_file=args.annotation,
            viability_file=args.gene_viability,
            output_file=args.output,
        )

        covered_genes = load_covered_genes(config.lfc_file, config.annotation_file)
        viability = load_gene_viability(config.viability_file)
        stats = compute_coverage_stats(viability, covered_genes)
        write_coverage_table(stats, config.output_file)

    except ValueError as exc:
        logger.error(f"Analysis failed: {exc}")
        return 1

    logger.success("Script completed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
