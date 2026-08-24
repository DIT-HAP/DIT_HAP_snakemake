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
Insertion-Level Aggregation Weights
===================================

Compute the per-insertion, per-timepoint weights that gene-level depletion
analysis uses to collapse insertion log2 fold changes into a gene LFC. This is
the single owner of the weighting algorithm: the aggregation script consumes the
``Weight`` column produced here and never transforms it further, so a new
weighting formula only has to be added in this file.

Which statistic the weights come from is decided by one flag,
``--has-replicates``. With replicates, DESeq2 supplies adjusted p-values and each
insertion is weighted by ``-log10(padj)``. Without replicates there are no
p-values at all, so the weights fall back on curve-fitting goodness of fit,
``-log10(1 - R2)``. Both are clipped into the open interval (1e-6, 1 - 1e-6)
before the log, which keeps a maximally-insignificant insertion at a tiny floor
weight rather than zero — a gene whose insertions all sit at that floor collapses
to an arithmetic mean instead of dividing by zero.

Weights are emitted in long format, keyed by insertion, timepoint and gene,
because a small number of insertions are annotated to two genes and the
normalising denominator is gene-specific: such an insertion carries a different
weight for each of its genes. Only in-gene insertions with an observed LFC take
part, and each gene-timepoint group's weights sum to 1.

Input
-----
- ``--stats`` (TSV): the source statistic, with a 4-level index
  (Chr, Coordinate, Strand, Target). With ``--has-replicates``, ``padj.tsv`` with
  one column per timepoint. With ``--no-has-replicates``, the curve-fitting
  statistics table carrying an ``R2`` column and one ``*_fitted`` column per
  timepoint.
- ``--lfc`` (TSV): insertion-level LFC, same 4-level index, one column per
  timepoint. Cells with no LFC are excluded from the normalising denominator.
- ``--annotations`` (TSV): genomic annotations, same 4-level index, with
  ``Type``, ``Distance_to_stop_codon`` and ``Systematic ID`` columns.

Output
------
- ``--output`` (TSV): long-format weights with columns Chr, Coordinate, Strand,
  Target, Timepoint, ``Systematic ID`` and ``Weight``.

Usage
-----
    python compute_insertion_weights.py --has-replicates \
        -s padj.tsv -l LFC.tsv -a annotations.tsv -o insertion_weights.tsv
    python compute_insertion_weights.py --no-has-replicates \
        -s insertion_level_fitting_statistics.tsv -l LFC.tsv \
        -a annotations.tsv -o insertion_weights.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-18
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from depletion.weights import (  # noqa: E402
    load_inputs,
    filter_in_gene,
    raw_weights,
    normalise_per_gene_timepoint,
    display_summary,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class WeightsConfig:
    """Validated inputs and branch selection for weight computation."""
    has_replicates: bool
    stats_file: Path
    lfc_file: Path
    annotations_file: Path
    output_file: Path

    def __post_init__(self) -> None:
        for path in (self.stats_file, self.lfc_file, self.annotations_file):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute insertion-level weights for gene-level depletion aggregation"
    )

    # Required rather than defaulted: the two branches expect different --stats
    # tables, so a silent default would mis-read whichever file it was handed.
    parser.add_argument("--has-replicates", action=argparse.BooleanOptionalAction, required=True,
                        help="Whether the upstream branch had biological replicates: "
                             "weight by -log10(padj) if so, by -log10(1 - R2) if not")
    parser.add_argument("-s", "--stats", type=Path, required=True,
                        help="Source statistic TSV: padj with replicates, curve-fitting "
                             "statistics without")
    parser.add_argument("-l", "--lfc", type=Path, required=True,
                        help="Insertion-level LFC TSV")
    parser.add_argument("-a", "--annotations", type=Path, required=True,
                        help="Genomic annotations TSV")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output long-format weights TSV")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable DEBUG level logging")

    return parser.parse_args()


def main() -> int:
    """Compute and write insertion-level aggregation weights."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = WeightsConfig(
            has_replicates=args.has_replicates,
            stats_file=args.stats,
            lfc_file=args.lfc,
            annotations_file=args.annotations,
            output_file=args.output,
        )
    except ValueError as exc:
        logger.error(f"Configuration error: {exc}")
        return 1

    try:
        branch = "with replicates" if config.has_replicates else "without replicates"
        logger.info(f"Computing insertion-level weights ({branch})")

        inputs = load_inputs(
            config.stats_file, config.lfc_file, config.annotations_file, config.has_replicates
        )
        inputs = filter_in_gene(inputs)
        weights = normalise_per_gene_timepoint(
            raw_weights(inputs, config.has_replicates), inputs
        )

        display_summary(weights)

        weights.to_csv(config.output_file, sep="\t")
        logger.success(f"Weights saved to {config.output_file}")

    except Exception as exc:
        logger.exception(f"An unexpected error occurred: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
