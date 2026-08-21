#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy",
#     "pandas",
#     "loguru",
# ]
# ///

"""
Gene-Level Depletion Analysis for Transposon Insertion Sequencing
=================================================================

Aggregates insertion-level log2 fold changes (LFC) to the gene level using a
pre-computed, gene-timepoint-normalised weight for each insertion. Weighting is
owned entirely upstream by ``compute_insertion_weights.py``: this script only
joins insertion LFC values onto that weight table and sums, so it has no
knowledge of whether the weights came from DESeq2 p-values or curve-fitting R2.

Because each insertion's weight already sums to 1 within its gene-timepoint
group, aggregation is a plain weighted sum rather than a weighted average —
mathematically identical to ``np.average`` with normalised weights, but
vectorised across every gene and timepoint at once.

Input
-----
- ``--lfc_path`` (TSV): insertion-level LFC, 4-level index (Chr, Coordinate,
  Strand, Target), one column per timepoint.
- ``--weights_path`` (TSV): long-format insertion weights from
  ``compute_insertion_weights.py``, with columns Chr, Coordinate, Strand,
  Target, Timepoint, Systematic ID, Weight (already gene-timepoint normalised).
- ``--annotations_path`` (TSV): genomic annotations, 4-level index, used only
  to look up each gene's ``Name``, ``FYPOviability`` and
  ``DeletionLibrary_essentiality``.

Output
------
- ``--output_path`` (TSV): gene-level statistics table (wide, per-timepoint LFC),
  indexed by ``Systematic ID``.
- ``LFC.tsv`` written alongside the output.

Usage
-----
    python gene_level_depletion_analysis.py -l lfc.tsv -a annotations.tsv -w insertion_weights.tsv -o gene_level_statistics.tsv
    python gene_level_depletion_analysis.py -l lfc.tsv -a annotations.tsv -w insertion_weights.tsv -o gene_level_statistics.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-18
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from depletion.gene_level import (  # noqa: E402
    GENE_ID,
    load_lfc_long,
    load_weights_long,
    load_gene_metadata,
    aggregate_to_gene_level,
    generate_summary,
    display_summary,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisConfig:
    """Immutable configuration holding validated input file paths."""
    lfc_path: Path
    weights_path: Path
    annotations_path: Path

    def __post_init__(self) -> None:
        for path in (self.lfc_path, self.weights_path, self.annotations_path):
            if not path.exists():
                raise ValueError(f"Input file {path} does not exist")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Gene-level depletion analysis for transposon insertion sequencing"
    )

    parser.add_argument('-l', '--lfc_path', type=Path, required=True,
                       help='Path to TSV file with insertion-level LFC')
    parser.add_argument('-a', '--annotations_path', type=Path, required=True,
                       help='Path to TSV file with annotations')
    parser.add_argument('-w', '--weights_path', type=Path, required=True,
                       help='Path to long-format insertion weights TSV')
    parser.add_argument('-o', '--output_path', type=Path, required=True,
                       help='Path for output TSV file')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable DEBUG level logging')

    return parser.parse_args()

def main() -> int:
    """Execute the gene-level depletion analysis."""
    start_time = time.time()

    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    logger.info("Starting gene-level depletion analysis")

    try:
        config = AnalysisConfig(
            lfc_path=args.lfc_path,
            weights_path=args.weights_path,
            annotations_path=args.annotations_path
        )

        args.output_path.parent.mkdir(parents=True, exist_ok=True)

        lfc_long = load_lfc_long(config.lfc_path)
        weights = load_weights_long(config.weights_path)
        gene_metadata = load_gene_metadata(config.annotations_path)

        gene_results = aggregate_to_gene_level(lfc_long, weights, gene_metadata)

        summary = generate_summary(gene_results)
        display_summary(summary)

        gene_results = gene_results.set_index(GENE_ID)
        gene_results.to_csv(args.output_path.parent / "LFC.tsv", sep="\t")
        gene_results.to_csv(args.output_path, sep="\t")

        elapsed = time.time() - start_time
        logger.info(f"Analysis completed in {elapsed:.1f}s")
        logger.info(f"Results saved to: {args.output_path}")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
