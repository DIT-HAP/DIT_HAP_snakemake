#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gene Coverage Figure (computation + rendering, single stage)
============================================================

One-stage QC figure: tally gene coverage by PomBase viability category directly
from the insertion-level LFC table, the per-insertion annotations, and the
gene-viability list — no intermediate coverage statistics TSV.

Computation (in-script, reusing ``qc.coverage``): genes with at least one
surviving insertion are counted as covered; the covered/total tally is
assembled per viability category in display order and rendered by the shared
composition renderer as a grouped bar chart followed by one donut per category.

Input
-----
- LFC TSV (``-l/--lfc``): insertion-level log fold change table with 4-level
  index (Chr, Coordinate, Strand, Target).
- Annotation TSV (``-a/--annotation``): per-insertion genomic annotations
  with Systematic ID column.
- Viability TSV (``-v/--viability``): PomBase gene viability categories,
  two-column headerless format (gene_id, viability).

Output
------
- ``<output>.pdf`` — journal-quality vector figure.
- ``<output>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_gene_coverage.py -l LFC.tsv -a annotations.tsv -v viability.tsv -o figures/gene_coverage
    python plot_gene_coverage.py -l LFC.tsv -a annotations.tsv -v viability.tsv -o figures/gene_coverage --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-28
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from qc.coverage import compute_coverage_stats, load_covered_genes, load_gene_viability  # noqa: E402
from figure_render.composition import render_composition_figure  # noqa: E402
from figure_render._schema import require_columns  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input paths and output stem."""
    lfc_file: Path
    annotation_file: Path
    viability_file: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate that every input exists and output directory can be created."""
        for path in (self.lfc_file, self.annotation_file, self.viability_file):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["category", "covered", "not_covered", "total", "coverage_pct"]

COVERED_LABEL = "Covered"
NOT_COVERED_LABEL = "Not covered"

# Horizontal bar: categories ride the y axis, percentages the x axis.
X_LABEL = "Coverage (%)"
Y_LABEL = "Gene viability"
TITLE = "Gene coverage by viability"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def load_coverage_data(input_path: Path) -> pd.DataFrame:
    """Load the archived coverage statistics TSV and validate its schema.

    Kept for tests and pixel-baseline tooling that operate on the archived
    table; the pipeline main() path computes the same frame in-memory via
    ``compute_coverage_frame`` instead.
    """
    df = pd.read_csv(input_path, sep="\t")
    require_columns(df, REQUIRED_COLUMNS, context=f"coverage TSV {input_path.name}")
    logger.info(f"Loaded {len(df)} viability categories")
    return df


@logger.catch(reraise=True)
def compute_coverage_frame(lfc_file: Path, annotation_file: Path, viability_file: Path) -> pd.DataFrame:
    """Tally covered vs total genes per viability category into the render-ready frame."""
    covered_genes = load_covered_genes(lfc_file, annotation_file)
    viability = load_gene_viability(viability_file)
    stats = compute_coverage_stats(viability, covered_genes)

    df = pd.DataFrame(
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
    logger.info(f"Computed coverage for {len(df)} viability categories")
    return df


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Render gene coverage figure from LFC, annotations, and viability tables")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Insertion-level LFC TSV")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Per-insertion annotation TSV (with Systematic ID column)")
    parser.add_argument("-v", "--gene-viability", type=Path, required=True, help="PomBase gene viability TSV (headerless: id, viability)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output file stem (extension will be added)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Compute coverage in memory, render the figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    output_stem = args.output.with_suffix("") if args.output.suffix == ".pdf" else args.output

    try:
        config = PlotConfig(
            lfc_file=args.input,
            annotation_file=args.annotation,
            viability_file=args.gene_viability,
            output_stem=output_stem,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Gene Coverage Figure Rendering ===")

    try:
        df = compute_coverage_frame(config.lfc_file, config.annotation_file, config.viability_file)

        render_composition_figure(
            df, config.output_stem,
            category_column="category", percentage_column="coverage_pct",
            part_column="covered", whole_column="not_covered",
            part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
            xlabel=X_LABEL, ylabel=Y_LABEL, title=TITLE,
        )
    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
