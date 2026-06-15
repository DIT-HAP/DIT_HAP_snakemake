"""
The Python script for imputing missing values using Forward-Reverse (FR) complementation.

This script performs imputation of missing values by leveraging complementary strand
data from genomic insertions, improving data completeness for downstream analysis.

Typical Usage:
    python impute_missing_values_using_FR.py --input counts_data.tsv --annotation insertion_annotations.tsv --output imputed_counts.tsv

Input: Counts data file with insertion information
Output: Imputed counts data with missing values filled using FR complementation
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Tuple, Dict
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd


# =============================== Constants ===============================
# Configuration constants for filtering
INTERGENIC_DISTANCE_THRESHOLD = 500
DISTANCE_TO_STOP_CODON_THRESHOLD = 3


# =============================== Configuration & Models ===============================
class ImputationConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_file: Path = Field(..., description="Path to the input counts data file")
    annotation_file: Path = Field(..., description="Path to the insertion annotation file")
    output_file: Path = Field(..., description="Path to the output imputed counts file")

    @field_validator('input_file')
    def validate_input_file(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        if not v.suffix == '.tsv':
            raise ValueError(f"Input file must be a TSV file: {v}")
        return v
    
    @field_validator('annotation_file')
    def validate_annotation_file(cls, v):
        if not v.exists():
            raise ValueError(f"Annotation file does not exist: {v}")
        if not v.suffix == '.tsv':
            raise ValueError(f"Annotation file must be a TSV file: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True


class ImputationResult(BaseModel):
    """Pydantic model to hold and validate the results of the imputation analysis."""
    total_insertions: int = Field(..., ge=0, description="Total number of insertions processed")
    in_gene_insertions: int = Field(..., ge=0, description="Number of insertions in coding genes")
    intergenic_insertions: int = Field(..., ge=0, description="Number of insertions in intergenic regions")
    complete_insertions_before: int = Field(..., ge=0, description="Number of insertions with all replicates available before imputation")
    complete_in_gene_before: int = Field(..., ge=0, description="Number of complete insertions in coding genes before imputation")
    complete_intergenic: int = Field(..., ge=0, description="Number of complete insertions in intergenic regions")
    complete_in_gene_after: int = Field(..., ge=0, description="Number of complete insertions in coding genes after imputation")
    imputed_insertions: int = Field(..., ge=0, description="Number of insertions imputed using FR complementation")
    complementarity_used: int = Field(..., ge=0, description="Number of complementary indices used for imputation")


# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# =============================== Core Functions ===============================
@logger.catch
def filter_insertions(insertion_annotations: pd.DataFrame) -> Tuple[pd.Index, pd.Index]:
    """
    Filter insertions into intergenic and in-gene categories based on genomic annotations.
    
    Args:
        insertion_annotations: DataFrame containing insertion annotations with genomic features
        
    Returns:
        Tuple of (intergenic_insertions_filtered, in_gene_insertions) indices
    """
    intergenic_insertions_filtered = insertion_annotations[
        (insertion_annotations["Type"] == "Intergenic region") & 
        (insertion_annotations["Distance_to_region_start"] > INTERGENIC_DISTANCE_THRESHOLD) & 
        (insertion_annotations["Distance_to_region_end"] > INTERGENIC_DISTANCE_THRESHOLD)
    ].index

    in_gene_insertions = insertion_annotations.query(
        f"Type != 'Intergenic region' & Distance_to_stop_codon > {DISTANCE_TO_STOP_CODON_THRESHOLD}"
    ).index

    return intergenic_insertions_filtered, in_gene_insertions


@logger.catch
def transfer_FR_index(idxs: tuple) -> tuple:
    """
    Transfer Forward-Reverse index by flipping the strand orientation.
    
    Args:
        idxs: Tuple representing the insertion index (chromosome, position, strand, ...)
        
    Returns:
        Tuple with the strand orientation flipped
    """
    idxs = list(idxs)
    idxs[2] = "+" if idxs[2] == "-" else "-"
    return tuple(idxs)


@logger.catch
def impute_missing_values(in_gene_counts_df: pd.DataFrame) -> Tuple[pd.DataFrame, List[tuple]]:
    """
    Impute missing values using Forward-Reverse complementation.
    
    Args:
        in_gene_counts_df: DataFrame containing insertion counts for coding genes
        
    Returns:
        Tuple of (imputed_dataframe, complementary_indices_used)
    """
    stacked_df = in_gene_counts_df.stack(level=0, dropna=False)
    stacked_dropna_df = in_gene_counts_df.stack(level=0, dropna=True)
    
    # Find indices with missing values
    in_gene_isna_idx = stacked_df[stacked_df.isna().all(axis=1)].index
    in_gene_complementary_idx = [transfer_FR_index(idx) for idx in in_gene_isna_idx]
    
    # Find indices that have complementary data available
    in_gene_index_for_imputation = list(set(in_gene_complementary_idx) & set(stacked_dropna_df.index))
    in_gene_has_complementary_idxs = [transfer_FR_index(idx) for idx in in_gene_index_for_imputation]

    # Perform imputation by transferring complementary data
    stacked_df.loc[in_gene_has_complementary_idxs, :] = stacked_df.loc[in_gene_index_for_imputation, :].values

    # Return unstacked dataframe with proper level ordering
    imputed_df = stacked_df.unstack().reorder_levels([1, 0], axis=1)
    
    return imputed_df, in_gene_has_complementary_idxs


@logger.catch
def calculate_imputation_statistics(counts_df: pd.DataFrame, in_gene_insertions: pd.Index, 
                                   intergenic_counts_df: pd.DataFrame, in_gene_counts_df: pd.DataFrame,
                                   imputed_in_gene_counts_df_noNA: pd.DataFrame,
                                   in_gene_has_complementary_idxs: List[tuple]) -> ImputationResult:
    """
    Calculate comprehensive statistics for the imputation process.
    
    Args:
        counts_df: Original counts dataframe
        in_gene_insertions: Insertions in coding genes
        intergenic_counts_df: Complete intergenic counts dataframe
        in_gene_counts_df: In-gene counts dataframe
        imputed_in_gene_counts_df_noNA: Imputed in-gene counts without missing values
        in_gene_has_complementary_idxs: List of complementary indices used for imputation
        
    Returns:
        ImputationResult object containing all statistics
    """
    insertion_num = counts_df.shape[0]
    ingene_num = in_gene_counts_df.shape[0]
    intergenic_num = counts_df[~counts_df.index.isin(in_gene_insertions)].shape[0]

    noNA_insertion_num = counts_df.dropna(axis=0, how="any").shape[0]
    noNA_ingene_num = in_gene_counts_df.dropna(axis=0, how="any").shape[0]
    noNA_intergenic_num = intergenic_counts_df.shape[0]

    noNA_imputed_ingene_num = imputed_in_gene_counts_df_noNA.shape[0]
    increased_ingene_num = noNA_imputed_ingene_num - noNA_ingene_num

    return ImputationResult(
        total_insertions=insertion_num,
        in_gene_insertions=ingene_num,
        intergenic_insertions=intergenic_num,
        complete_insertions_before=noNA_insertion_num,
        complete_in_gene_before=noNA_ingene_num,
        complete_intergenic=noNA_intergenic_num,
        complete_in_gene_after=noNA_imputed_ingene_num,
        imputed_insertions=increased_ingene_num,
        complementarity_used=len(in_gene_has_complementary_idxs)
    )


@logger.catch
def print_imputation_statistics(result: ImputationResult) -> None:
    """
    Print detailed statistics of the imputation process.
    
    Args:
        result: ImputationResult containing all statistics
    """
    logger.info("### Impute missing values using FR completed ###")
    
    logger.info(f"*** Total insertions: {result.total_insertions}")
    logger.info(f"*** Insertions in coding genes: {result.in_gene_insertions} ({result.in_gene_insertions/result.total_insertions*100:.2f}%)")
    logger.info(f"*** Insertions in intergenic regions: {result.intergenic_insertions} ({result.intergenic_insertions/result.total_insertions*100:.2f}%)")

    logger.info(f"*** Insertions with all replicates available: {result.complete_insertions_before} ({result.complete_insertions_before/result.total_insertions*100:.2f}%)")
    logger.info(f"*** Insertions with all replicates available in coding genes: {result.complete_in_gene_before} ({result.complete_in_gene_before/result.complete_insertions_before*100:.2f}%) - Compared with insertions with all replicates available")
    logger.info(f"*** Insertions with all replicates available in intergenic regions: {result.complete_intergenic} ({result.complete_intergenic/result.complete_insertions_before*100:.2f}%) - Compared with insertions with all replicates available")

    logger.info(f"*** Insertions with all replicates available in coding genes after imputation: {result.complete_in_gene_after} ({result.complete_in_gene_after/result.complete_insertions_before*100:.2f}%)")
    logger.info(f"*** Increase in insertions with all replicates available in coding genes: {result.imputed_insertions} ({result.imputed_insertions/result.complete_insertions_before*100:.2f}%)")
    logger.info(f"*** Insertions with all replicates available in coding genes after imputation (in all in-gene insertions): {result.complete_in_gene_after} ({result.complete_in_gene_after/result.in_gene_insertions*100:.2f}%)")
    logger.info(f"*** Complementary indices used for imputation: {result.complementarity_used}")


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments for the imputation script."""
    parser = argparse.ArgumentParser(description="Impute missing values using Forward-Reverse complementation.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input counts data file")
    parser.add_argument("-a", "--annotation", type=Path, required=True, help="Path to the insertion annotation file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output imputed counts file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


@logger.catch
def main():
    """Main entry point of the script for imputing missing values using FR complementation."""

    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
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

        samples = counts_df.columns.get_level_values(0).unique()
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
        imputation_statistics.to_csv(config.output_file.parent/"imputation_statistics.tsv", index=True, sep="\t")
        logger.success(f"Number of imputed insertions saved to {config.output_file.parent/'imputation_statistics.tsv'}")

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
        sys.exit(1)

if __name__ == "__main__":
    main()