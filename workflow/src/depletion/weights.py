"""Insertion-Level Aggregation Weights

Compute the per-insertion, per-timepoint weights that gene-level depletion
analysis uses to collapse insertion log2 fold changes into a gene LFC. This is
the single owner of the weighting algorithm: the aggregation script consumes the
``Weight`` column produced here and never transforms it further, so a new
weighting formula only has to be added in this file.

Which statistic the weights come from is decided by one thing — whether the
upstream branch had biological replicates, passed in as ``has_replicates``. With
replicates, DESeq2 supplies adjusted p-values and each insertion is weighted by
``-log10(padj)``. Without replicates there are no p-values at all, so the weights
fall back on curve-fitting goodness of fit, ``-log10(1 - R2)``. Both are clipped
into the open interval (1e-6, 1 - 1e-6) before the log, which keeps a
maximally-insignificant insertion at a tiny floor weight rather than zero — a
gene whose insertions all sit at that floor collapses to an arithmetic mean
instead of dividing by zero.

Weights are emitted in long format, keyed by insertion, timepoint and gene,
because a small number of insertions are annotated to two genes and the
normalising denominator is gene-specific: such an insertion carries a different
weight for each of its genes. Only in-gene insertions with an observed LFC take
part, and each gene-timepoint group's weights sum to 1.

Input
-----
- ``stats`` (DataFrame): the source statistic, with a 4-level index
  (Chr, Coordinate, Strand, Target). With replicates, ``padj.tsv`` with one
  column per timepoint. Without replicates, the curve-fitting statistics table
  carrying an ``R2`` column and one ``*_fitted`` column per timepoint.
- ``lfc`` (DataFrame): insertion-level LFC, same 4-level index, one column per
  timepoint. Cells with no LFC are excluded from the normalising denominator.
- ``annotations`` (DataFrame): genomic annotations, same 4-level index, with
  ``Type``, ``Distance_to_stop_codon`` and ``Systematic ID`` columns.

Output
------
- Long-format weights with columns Chr, Coordinate, Strand, Target, Timepoint,
  ``Systematic ID`` and ``Weight``.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-18
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from io_tables import read_insertion_table

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
INDEX_COLUMNS = [0, 1, 2, 3]

# Clip bounds for every quantity that goes through -log10. The upper bound is
# what keeps a maximally-insignificant insertion (padj == 1, or R2 <= 0) at a
# tiny positive floor weight instead of exactly zero, so a gene whose insertions
# all sit there becomes an arithmetic mean rather than a 0/0 division.
PROB_FLOOR = 1e-6
PROB_CEILING = 1 - 1e-6

TIMEPOINT_AXIS = "Timepoint"
WEIGHT_COLUMN = "Weight"
FITTED_SUFFIX = "_fitted"


class AnnotCol(StrEnum):
    """Annotation columns used for in-gene filtering and gene grouping."""
    TYPE = "Type"
    DISTANCE_TO_STOP = "Distance_to_stop_codon"
    GENE_ID = "Systematic ID"


class StatsCol(StrEnum):
    """Columns read from the curve-fitting statistics table."""
    R2 = "R2"

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class WeightInputs:
    """Source tables sharing a 4-level insertion index."""
    stats: pd.DataFrame
    lfc: pd.DataFrame
    annotations: pd.DataFrame

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
# --- Data Loading ---
@logger.catch
def load_inputs(stats_file: Path, lfc_file: Path, annotations_file: Path, has_replicates: bool) -> WeightInputs:
    """Load the source statistic, insertion LFC and annotations."""
    statistic = "padj" if has_replicates else "curve-fitting R2"
    logger.info(f"Loading {statistic} source statistic from {stats_file}")
    logger.info(f"Loading insertion LFC from {lfc_file}")
    logger.info(f"Loading annotations from {annotations_file}")

    inputs = WeightInputs(
        stats=read_insertion_table(stats_file),
        lfc=read_insertion_table(lfc_file),
        annotations=read_insertion_table(annotations_file),
    )

    logger.info(
        f"Loaded {len(inputs.lfc):,} insertions, "
        f"{len(inputs.annotations):,} annotation rows"
    )
    return inputs


@logger.catch
def filter_in_gene(inputs: WeightInputs) -> WeightInputs:
    """Restrict the LFC and stats tables to in-gene insertions.

    The annotations table is deliberately left FULL and unfiltered. Gene-level
    depletion analysis has historically joined the (LFC, Weight) stack against
    every annotation row sharing an insertion's index, not just the row(s) that
    individually pass the in-gene query. That reproduces two long-standing
    quirks this script must match bit-for-bit: (1) a site with more than one
    annotation row for the same gene (e.g. one row per transcript feature) has
    its weight counted once per row, and (2) a site shared by two genes is
    admitted whenever ANY of its rows passes the filter, even for the gene
    whose own row does not. See docs/plans/2026-08-18-decouple-insertion-weights-design.md.
    """
    in_gene_index = inputs.annotations.query(
        f"(`{AnnotCol.TYPE}` != 'Intergenic region') "
        f"and (`{AnnotCol.DISTANCE_TO_STOP}` > 4)"
    ).index

    filtered = WeightInputs(
        stats=inputs.stats[inputs.stats.index.isin(in_gene_index)].copy(),
        lfc=inputs.lfc[inputs.lfc.index.isin(in_gene_index)].copy(),
        annotations=inputs.annotations,
    )

    logger.info(
        f"In-gene insertions retained: {len(filtered.lfc):,} of {len(inputs.lfc):,}"
    )
    return filtered

# --- Weight Formulas ---
def neg_log10(values: pd.DataFrame) -> pd.DataFrame:
    """Clip a probability-like frame into (0, 1) and take ``-log10``.

    fillna runs before clip on purpose: reversed, a NaN would skip the ceiling
    and -log10(1) would give it a zero weight, silencing the insertion instead
    of leaving it at the floor.
    """
    return -np.log10(values.fillna(1).clip(lower=PROB_FLOOR, upper=PROB_CEILING))


@logger.catch
def r2_confidence_frame(stats: pd.DataFrame) -> pd.DataFrame:
    """Broadcast each insertion's clipped R2 across the fitted timepoint columns."""
    # Timepoint names come from the "*_fitted" columns; removesuffix is exact
    # (unlike rstrip, which strips any trailing character from the set).
    timepoints = [
        column.removesuffix(FITTED_SUFFIX)
        for column in stats.filter(regex=rf".*{FITTED_SUFFIX}$").columns
    ]
    if not timepoints:
        raise ValueError(f"No '*{FITTED_SUFFIX}' columns found in the fitting statistics")
    logger.info(f"Timepoints detected from fitted columns: {timepoints}")

    # A negative R2 (fit worse than a flat mean) clips to the floor, so those
    # insertions land on the minimum weight rather than a negative one.
    confidence = stats[StatsCol.R2].clip(lower=PROB_FLOOR, upper=PROB_CEILING)
    negative_r2 = int((stats[StatsCol.R2] < 0).sum())
    if negative_r2:
        logger.info(f"{negative_r2:,} insertions have R2 < 0 and clip to the floor")

    return pd.DataFrame(
        {timepoint: 1 - confidence for timepoint in timepoints},
        index=stats.index,
    )


@logger.catch
def raw_weights(inputs: WeightInputs, has_replicates: bool) -> pd.DataFrame:
    """Compute un-normalised insertion-by-timepoint weights.

    With replicates the stats table is already one padj column per timepoint, so
    it goes straight through ``-log10``. Without replicates there is a single R2
    per insertion, which first has to be broadcast across the timepoints.
    """
    if has_replicates:
        weights = neg_log10(inputs.stats)
        label = "padj"
    else:
        weights = neg_log10(r2_confidence_frame(inputs.stats))
        label = "1 - R2"

    values = weights.to_numpy()
    logger.info(
        f"-log10({label}): raw weights mean={np.nanmean(values):.4g}, "
        f"min={np.nanmin(values):.4g}, max={np.nanmax(values):.4g}"
    )
    return weights

# --- Normalisation ---
@logger.catch
def normalise_per_gene_timepoint(weights: pd.DataFrame, inputs: WeightInputs) -> pd.DataFrame:
    """Reshape weights to long format and normalise within each gene-timepoint."""
    # Only cells with an observed LFC take part, so the denominator covers
    # exactly the rows gene-level aggregation will consume.
    weights = weights.where(inputs.lfc.notna())

    long = weights.rename_axis(TIMEPOINT_AXIS, axis=1).stack().to_frame(WEIGHT_COLUMN)
    merged = long.join(inputs.annotations[[AnnotCol.GENE_ID.value]], how="left")
    merged = merged[merged[WEIGHT_COLUMN].notna() & merged[AnnotCol.GENE_ID].notna()]

    group_keys = [AnnotCol.GENE_ID.value, TIMEPOINT_AXIS]
    weight_sums = merged.groupby(group_keys)[WEIGHT_COLUMN].transform("sum")

    # A zero sum would divide by zero. The clip ceiling in neg_log10 makes this
    # unreachable today (every weight is strictly positive), so treat it as a
    # loud failure rather than silently substituting uniform weights: a future
    # weighting formula needs to decide its own fallback deliberately.
    zero_sum_groups = int(merged.loc[weight_sums == 0].groupby(group_keys).ngroups)
    if zero_sum_groups:
        raise ValueError(
            f"{zero_sum_groups:,} gene-timepoint groups have a zero weight sum; "
            f"the weighting formula needs an explicit fallback before it can be normalised"
        )

    merged[WEIGHT_COLUMN] /= weight_sums

    logger.info(
        f"Normalised {len(merged):,} insertion-timepoint-gene rows across "
        f"{merged.groupby(group_keys).ngroups:,} gene-timepoints"
    )
    return merged

# --- Summary ---
def display_summary(weights: pd.DataFrame) -> None:
    """Log the shape and per-gene spread of the normalised weights."""
    per_group = weights.groupby([AnnotCol.GENE_ID.value, TIMEPOINT_AXIS])[WEIGHT_COLUMN]

    logger.info("=" * 60)
    logger.info("INSERTION-LEVEL WEIGHTS SUMMARY")
    logger.info("=" * 60)
    logger.info(f"{'Rows (insertion x timepoint x gene)':<40}: {len(weights):,}")
    logger.info(f"{'Genes':<40}: {weights[AnnotCol.GENE_ID].nunique():,}")
    logger.info(f"{'Timepoints':<40}: {weights.index.get_level_values(TIMEPOINT_AXIS).nunique()}")
    logger.info(f"{'Gene-timepoint groups':<40}: {per_group.ngroups:,}")
    logger.info(f"{'Mean insertions per group':<40}: {per_group.size().mean():.1f}")
    logger.info(f"{'Max single-insertion weight':<40}: {weights[WEIGHT_COLUMN].max():.4f}")
    logger.info("=" * 60)
