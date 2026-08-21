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
#     "pybedtools",
# ]
# ///

"""
Genomic Feature Annotation of Transposon Insertion Sites
========================================================

Annotates transposon insertion sites with genomic features including genes,
intergenic regions, and coding sequences. Insertion coordinates are intersected
against a genome-region BED annotation using pybedtools, then distances to the
start/stop codon, the affected amino-acid residue, the reading frame, and the
insertion direction relative to each gene are computed.

Boundary duplicates (insertions falling exactly on a region edge, or matching
both a coding and an intergenic region) are resolved so that each insertion
keeps a single, most-specific annotation.

Input
-----
- Insertion file (TSV or CSV) with columns: Chr, Coordinate, Strand, Target.
- Genome-region BED file (TSV) with region intervals and their metadata.

Output
------
- Tab-separated annotation table with per-insertion genomic features, codon
  distances, affected residues, and insertion direction.

Usage
-----
    python annotate_genomic_features.py -i insertions.tsv -g genome_region.bed -o annotated.tsv
    python annotate_genomic_features.py --input insertions.tsv --genome-region genome_region.bed --output annotated.tsv --verbose

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
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from read_processing.annotation import (  # noqa: E402
    AnalysisResult,
    load_insertion_data,
    load_genome_regions,
    annotate_insertions,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths for the annotation workflow."""
    input_file: Path
    genome_region_file: Path
    output_file: Path

    def __post_init__(self) -> None:
        # Validation and output-directory creation (no attribute assignment,
        # so this is safe on a frozen dataclass).
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        if not self.genome_region_file.exists():
            raise ValueError(f"Input file does not exist: {self.genome_region_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Annotate insertion sites with genomic features")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input insertion file")
    parser.add_argument("-g", "--genome-region", type=Path, required=True, help="Path to the genome region BED file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output annotation file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point: validate configuration and run the annotation workflow."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    logger.info(f"Pandas version: {pd.__version__}")
    logger.info(f"NumPy version: {np.__version__}")

    try:
        # Create and validate configuration
        config = InputOutputConfig(
            input_file=args.input,
            genome_region_file=args.genome_region,
            output_file=args.output,
        )

        # Run the core analysis/logic
        logger.info(f"Starting annotation workflow for {config.input_file}")
        logger.info(f"Genome regions: {config.genome_region_file}")
        logger.info(f"Output file: {config.output_file}")

        # Load data
        insertions_df = load_insertion_data(config.input_file)
        regions_df = load_genome_regions(config.genome_region_file)

        # Annotate insertions
        annotated_df, results = annotate_insertions(insertions_df, regions_df)

        # Save results
        logger.info(f"Saving annotations to {config.output_file}")

        # Save with proper formatting
        annotated_df.to_csv(
            config.output_file,
            index=False,
            header=True,
            float_format="%.3f",
            sep="\t",
        )

        # Display statistics
        logger.info("=" * 60)
        logger.info("ANNOTATION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Total insertions: {results.total_insertions:,}")
        logger.success(f"Annotated insertions: {results.annotated_insertions:,}")

        logger.info("\nFeature distribution:")
        logger.info(f"  Coding regions: {results.coding_insertions:,} ({results.coding_percentage:.1f}%)")
        logger.info(f"  Intergenic regions: {results.intergenic_insertions:,}")

        if results.unique_genes > 0:
            logger.info("\nGene impact:")
            logger.info(f"  Unique genes affected: {results.unique_genes:,}")
            logger.info(f"  Forward insertions: {results.forward_insertions:,}")
            logger.info(f"  Reverse insertions: {results.reverse_insertions:,}")

        logger.success(f"Annotations saved to {config.output_file}")

        logger.success("Annotation workflow completed successfully!")

        logger.success(f"Annotation completed successfully! Processed {results.total_insertions} insertions.")
        return 0

    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
