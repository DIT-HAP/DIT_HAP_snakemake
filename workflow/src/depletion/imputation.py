"""Impute missing insertion counts using Forward-Reverse (FR) complementation.

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
- ``ImputationResult`` dataclass and the imputed counts DataFrame, as prepared
  by the invoking script.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass

import pandas as pd
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# Configuration constants for filtering
INTERGENIC_DISTANCE_THRESHOLD = 500
DISTANCE_TO_STOP_CODON_THRESHOLD = 3

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ImputationResult:
    """Summary statistics describing the FR imputation outcome."""
    total_insertions: int
    in_gene_insertions: int
    intergenic_insertions: int
    complete_insertions_before: int
    complete_in_gene_before: int
    complete_intergenic: int
    complete_in_gene_after: int
    imputed_insertions: int
    complementarity_used: int


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def filter_insertions(insertion_annotations: pd.DataFrame) -> tuple[pd.Index, pd.Index]:
    """Split insertions into intergenic and in-gene index sets based on annotations."""
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
    """Return the insertion index with its strand orientation flipped (+ <-> -)."""
    idxs = list(idxs)
    idxs[2] = "+" if idxs[2] == "-" else "-"
    return tuple(idxs)


@logger.catch
def impute_missing_values(in_gene_counts_df: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple]]:
    """Impute all-missing insertion rows from their complementary-strand counterpart."""
    stacked_df = in_gene_counts_df.stack(level=0)
    stacked_dropna_df = stacked_df.dropna(how="all")

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
                                    in_gene_has_complementary_idxs: list[tuple]) -> ImputationResult:
    """Compute summary counts describing insertions before and after imputation."""
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
    """Log detailed statistics of the imputation process."""
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
