#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "matplotlib",
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

Input
-----
- ``-i`` insertion-level LFC TSV, row MultiIndex (Chr, Coordinate, Strand, Target).
- ``-a`` per-insertion annotation TSV with a ``Systematic ID`` column and the
  same 4-column insertion key.
- ``-v`` PomBase gene viability TSV, two columns, no header:
  ``systematic_id`` and ``viability`` (viable / inviable / condition-dependent / unknown).

Output
------
- ``-o`` a multi-page PDF: page 1 is a grouped bar chart of coverage percentage
  per viability category; each following page is a donut chart of covered vs
  not-covered genes for one category.

Usage
-----
    python gene_coverage_analysis.py -i LFC.tsv -a annotations.tsv -v gene_viability.tsv -o gene_coverage_analysis.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
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
import matplotlib.pyplot as plt
from loguru import logger
from matplotlib.backends.backend_pdf import PdfPages

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
STYLE_PATH = SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle"

INSERTION_KEY = ["Chr", "Coordinate", "Strand", "Target"]

# Fixed display order for viability categories; any unlisted label is appended.
VIABILITY_ORDER = ["viable", "inviable", "condition-dependent", "unknown"]

COVERED_COLOR = "#962955"      # deep pink-purple
NOT_COVERED_COLOR = "#7fb775"  # medium green


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
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

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
def plot_coverage_bar(stats: list[CoverageStat], ax: plt.Axes) -> None:
    """Draw a grouped bar chart of coverage percentage per viability category."""
    categories = [s.category for s in stats]
    percentages = [s.coverage_pct for s in stats]

    ax.bar(categories, percentages, color=COVERED_COLOR)
    ax.set_ylabel("Coverage (%)")
    ax.set_xlabel("Gene viability")
    ax.set_title("Gene coverage by viability")
    ax.set_ylim(0, 100)
    for idx, stat in enumerate(stats):
        ax.text(idx, stat.coverage_pct, f"{stat.coverage_pct:.1f}%", ha="center", va="bottom")


@logger.catch
def plot_coverage_donut(stat: CoverageStat, ax: plt.Axes) -> None:
    """Draw a donut chart of covered vs not-covered genes for one category."""
    ax.pie(
        [stat.covered, stat.not_covered],
        labels=["Covered", "Not covered"],
        colors=[COVERED_COLOR, NOT_COVERED_COLOR],
        autopct=lambda pct: f"{pct:.1f}%",
        wedgeprops={"width": 0.4},
        startangle=90,
    )
    ax.set_title(f"{stat.category}\n({stat.covered:,}/{stat.total:,} genes)")


@logger.catch
def write_report(stats: list[CoverageStat], output_file: Path) -> None:
    """Write the multi-page coverage PDF: bar-chart overview then per-category donuts."""
    plt.style.use(STYLE_PATH)
    ax_width, ax_height = plt.rcParams["figure.figsize"]

    with PdfPages(output_file) as pdf:
        fig, ax = plt.subplots(figsize=(ax_width, ax_height))
        plot_coverage_bar(stats, ax)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for stat in stats:
            fig, ax = plt.subplots(figsize=(ax_width, ax_height))
            plot_coverage_donut(stat, ax)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    logger.success(f"Wrote coverage report with {len(stats) + 1} pages to {output_file}")

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
  python gene_coverage_analysis.py -i LFC.tsv -a annotations.tsv -v gene_viability.tsv -o gene_coverage_analysis.pdf
        """,
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Insertion-level LFC TSV")
    parser.add_argument("-a", "--annotation", type=Path, required=True,
                        help="Per-insertion annotation TSV (with Systematic ID column)")
    parser.add_argument("-v", "--gene-viability", type=Path, required=True,
                        help="PomBase gene viability TSV (headerless: id, viability)")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output PDF path")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: load coverage, tally by viability, write the PDF report."""
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
        write_report(stats, config.output_file)

    except ValueError as exc:
        logger.error(f"Analysis failed: {exc}")
        return 1

    logger.success("Script completed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
