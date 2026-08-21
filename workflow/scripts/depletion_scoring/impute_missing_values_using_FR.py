#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
#     "pandas",
# ]
# ///

"""
Impute Missing Values Using Forward-Reverse (FR) Complementation
================================================================

Fills missing insertion count values by borrowing data from the complementary
strand. Transposon insertions occur on both the forward (+) and reverse (-)
strands at the same genomic locus; when one strand's replicate counts are all
missing while the opposite-strand insertion has complete data, the missing
values are imputed from that opposite-strand insertion.

Only in-gene insertions are imputed. Complete intergenic insertions are carried
through unchanged and concatenated with the imputed in-gene set to form the
final complete-data table used by downstream depletion analysis.

Input
-----
- Counts TSV with a 4-level row MultiIndex (Chr, Coordinate, Strand, Target) and
  a 2-level column MultiIndex (sample, timepoint); read with header=[0, 1].
- Insertion annotation TSV with the same 4-level row MultiIndex, providing the
  genomic feature Type and distance columns used to separate in-gene from
  intergenic insertions.

Output
------
- Imputed counts TSV (same index/column structure as the input) written to the
  --output path.
- imputation_statistics.tsv written alongside the output, recording the number
  of imputed replicate values per retained insertion.

Usage
-----
    python impute_missing_values_using_FR.py -i counts.tsv -a annotations.tsv -o imputed_counts.tsv
    python impute_missing_values_using_FR.py -i counts.tsv -a annotations.tsv -o imputed_counts.tsv --verbose

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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402
from depletion.imputation import (  # noqa: E402
    ImputationResult,
    filter_insertions,
    transfer_FR_index,
    impute_missing_values,
    calculate_imputation_statistics,
    print_imputation_statistics,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ImputationConfig:
    """Validated input/output paths for the imputation run."""
    input_file: Path
    annotation_file: Path
    output_file: Path

    def __post_init__(self) -> None:
        """Validate that inputs exist as TSVs and ensure the output directory exists."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        if self.input_file.suffix != ".tsv":
            raise ValueError(f"Input file must be a TSV file: {self.input_file}")
        if not self.annotation_file.exists():
            raise ValueError(f"Annotation file does not exist: {self.annotation_file}")
        if self.annotation_file.suffix != ".tsv":
            raise ValueError(f"Annotation file must be a TSV file: {self.annotation_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments for the imputation script."""
    parser = argparse.ArgumentParser(description="Impute missing values using Forward-Reverse complementation.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input counts data file")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Path to the insertion annotation file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output imputed counts file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point of the script for imputing missing values using FR complementation."""
    args = parse_args()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(log_level)

    # Validate input and output paths using the dataclass config
    try:
        config = ImputationConfig(
            input_file=args.input,
            annotation_file=args.annotation,
            output_file=args.output
        )

        logger.info(f"Starting imputation of {config.input_file}")

        # Load input data
        counts_df = pd.read_csv(config.input_file, index_col=[0, 1, 2, 3], sep="\t", header=[0, 1])
        insertion_annotations = pd.read_csv(config.annotation_file, index_col=[0, 1, 2, 3], sep="\t")

        timepoints = counts_df.columns.get_level_values(1).unique()

        # Filter insertions by genomic location
        _, in_gene_insertions = filter_insertions(insertion_annotations)

        # Separate in-gene and intergenic insertions
        in_gene_counts_df = counts_df[counts_df.index.isin(in_gene_insertions)].copy()
        intergenic_counts_df = counts_df[~counts_df.index.isin(in_gene_insertions)].copy().dropna(axis=0, how="any")

        logger.info(f"Insertions with all replicates available in intergenic regions: {intergenic_counts_df.shape[0]}")
        logger.info(f"Insertions with at least one replicate available in coding genes: {in_gene_counts_df.shape[0]}")
        logger.info(f"Insertions with all replicates available in coding genes: {in_gene_counts_df.dropna(axis=0, how='any').shape[0]}")

        # Perform imputation
        imputed_in_gene_counts_df, in_gene_has_complementary_idxs = impute_missing_values(in_gene_counts_df)

        # Remove remaining missing values after imputation
        imputed_in_gene_counts_df_noNA = imputed_in_gene_counts_df.dropna(axis=0, how="any")
        logger.info(f"Insertions with all replicates available in coding genes after imputation: {imputed_in_gene_counts_df_noNA.shape[0]}")

        # Concatenate complete datasets
        imputed_counts = pd.concat([intergenic_counts_df, imputed_in_gene_counts_df_noNA], axis=0)

        imputation_statistics = counts_df.loc[imputed_counts.index].xs(timepoints[0], level=1, axis=1).isna().sum(axis=1).astype(int).rename("num_of_imputed_insertions")
        imputation_statistics.to_csv(config.output_file.parent / "imputation_statistics.tsv", index=True, sep="\t")
        logger.success(f"Number of imputed insertions saved to {config.output_file.parent / 'imputation_statistics.tsv'}")

        logger.info(f"Total insertions with all replicates available after imputation: {imputed_counts.value_counts()}")

        # only imputate less than once
        # imputated_LTonce_idx = imputation_statistics[imputation_statistics <= 1].index

        # Save imputed datas
        # imputed_counts.loc[imputated_LTonce_idx].to_csv(config.output_file, index=True, sep="\t")
        imputed_counts.to_csv(config.output_file, index=True, sep="\t")
        logger.success(f"Imputation complete. Results saved to {config.output_file}")

        # Calculate and print statistics
        stats = calculate_imputation_statistics(
            counts_df, in_gene_insertions, intergenic_counts_df, in_gene_counts_df,
            imputed_in_gene_counts_df_noNA, in_gene_has_complementary_idxs
        )
        print_imputation_statistics(stats)

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
