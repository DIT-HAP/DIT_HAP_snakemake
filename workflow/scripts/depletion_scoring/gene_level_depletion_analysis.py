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
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

TIMEPOINT_AXIS = "Timepoint"
LFC_COLUMN = "LFC"
WEIGHT_COLUMN = "Weight"
CONTRIBUTION_COLUMN = "Contribution"

GENE_ID = "Systematic ID"
GENE_METADATA_COLUMNS = ["Name", "FYPOviability", "DeletionLibrary_essentiality"]
GENE_GROUP_COLUMNS = [GENE_ID, *GENE_METADATA_COLUMNS]

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
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
# --- Data Loading ---
@logger.catch
def load_lfc_long(lfc_path: Path) -> pd.DataFrame:
    """Load insertion-level LFC and reshape to long format (site, Timepoint, LFC)."""
    logger.info(f"Loading LFC data from {lfc_path}")
    lfc_df = pd.read_csv(lfc_path, index_col=INSERTION_INDEX_COLUMNS, sep="\t")
    logger.info(f"Loaded {lfc_df.shape[0]:,} insertions")
    return lfc_df.rename_axis(TIMEPOINT_AXIS, axis=1).stack().to_frame(LFC_COLUMN)


@logger.catch
def load_weights_long(weights_path: Path) -> pd.DataFrame:
    """Load the pre-normalised, long-format insertion weight table."""
    logger.info(f"Loading insertion weights from {weights_path}")
    weights = pd.read_csv(weights_path, sep="\t")
    logger.info(f"Loaded {len(weights):,} insertion-timepoint-gene weight rows")
    return weights


@logger.catch
def load_gene_metadata(annotations_path: Path) -> pd.DataFrame:
    """Load one metadata row per gene (Name, FYPOviability, DeletionLibrary_essentiality)."""
    logger.info(f"Loading annotations from {annotations_path}")
    annotations_df = pd.read_csv(annotations_path, index_col=INSERTION_INDEX_COLUMNS, sep="\t")

    # PomBase releases may leave DeletionLibrary_essentiality blank for genes that were
    # previously annotated as "Not_determined". A bare NaN causes pandas groupby to
    # silently drop those genes (default dropna=True). Treat missing values the same
    # way as the explicit "Not_determined" label so no gene is lost.
    na_ess = annotations_df["DeletionLibrary_essentiality"].isna().sum()
    if na_ess > 0:
        logger.warning(
            f"Found {na_ess} insertions with missing DeletionLibrary_essentiality; "
            "filling with 'Not_determined' to prevent silent gene loss in groupby"
        )
        annotations_df["DeletionLibrary_essentiality"] = (
            annotations_df["DeletionLibrary_essentiality"].fillna("Not_determined")
        )

    return annotations_df[GENE_GROUP_COLUMNS].drop_duplicates(GENE_ID)

# --- Gene-Level Aggregation ---
@logger.catch
def aggregate_to_gene_level(lfc_long: pd.DataFrame, weights: pd.DataFrame,
                           gene_metadata: pd.DataFrame) -> pd.DataFrame:
    """Weighted-sum insertion LFC to gene level; each gene-timepoint's weights already sum to 1."""
    merged = weights.merge(
        lfc_long, on=[*lfc_long.index.names[:-1], TIMEPOINT_AXIS], how="inner"
    )
    logger.info(f"Matched {len(merged):,} of {len(weights):,} weight rows to an observed LFC")

    merged = merged.merge(gene_metadata, on=GENE_ID, how="left")
    merged[CONTRIBUTION_COLUMN] = merged[LFC_COLUMN] * merged[WEIGHT_COLUMN]

    gene_lfc = merged.groupby([*GENE_GROUP_COLUMNS, TIMEPOINT_AXIS])[CONTRIBUTION_COLUMN].sum(min_count=1)
    gene_wide = gene_lfc.unstack(TIMEPOINT_AXIS).round(3).dropna(how="all")

    logger.info(f"Completed analysis for {len(gene_wide):,} genes")
    return gene_wide.reset_index()

# --- Summary Statistics ---
def generate_summary(gene_df: pd.DataFrame) -> dict[str, int]:
    """Generate summary statistics for the analysis."""
    return {
        'Total genes analyzed': len(gene_df),
        'FYPOviability: Essential genes': len(gene_df[gene_df['FYPOviability'] == 'inviable']),
        'FYPOviability: Non-essential genes': len(gene_df[gene_df['FYPOviability'] == 'viable']),
        'DeletionLibrary_essentiality: Essential genes': len(gene_df[gene_df['DeletionLibrary_essentiality'] == 'E']),
        'DeletionLibrary_essentiality: Non-essential genes': len(gene_df[gene_df['DeletionLibrary_essentiality'] == 'V'])
    }

def display_summary(stats: dict[str, int]) -> None:
    """Display summary statistics."""
    logger.info("\n" + "="*60)
    logger.info("GENE-LEVEL DEPLETION ANALYSIS SUMMARY")
    logger.info("="*60)

    for key, value in stats.items():
        logger.info(f"{key:<40}: {value}")

    logger.info("="*60)

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
