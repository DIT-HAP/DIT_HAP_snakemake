"""Hard filtering of insertion reads by read-count threshold.

Filter insertion reads by a minimum read-count threshold measured at a chosen
initial timepoint. Each sample is processed independently: an insertion is
retained for a sample only when its read count at the initial timepoint meets
or exceeds the cutoff, so a sparsely-covered insertion in one sample does not
suppress the same locus in a well-covered sample.

The core algorithm groups the wide, multi-indexed count matrix by ``Sample``
(the first column level), builds a boolean mask from the ``>= cutoff`` test on
the initial-timepoint column, keeps the masked rows per sample, and concatenates
the surviving per-sample blocks back into a single matrix. Retention statistics
(total, retained, removed, retention rate, samples processed) are logged for
reproducibility.

Input
-----
- A tab-separated counts matrix with a 4-level row MultiIndex (columns 0-3) and
  a 2-level column MultiIndex ``(Sample, Timepoint)`` on header rows 0 and 1.

Output
------
- A tab-separated matrix of retained insertions, written with the row index and
  the 2-level column header preserved (``sep="\\t"``, ``header=True``,
  ``index=True``). Same structure as the input, filtered by row.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Results of the filtering analysis."""
    total_insertions: int
    retained_insertions: int
    removed_insertions: int
    retention_rate: float
    samples_processed: int

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_insertion_data(input_file: Path) -> pd.DataFrame:
    """Load insertion reads from a TSV file with multi-index rows and columns."""
    logger.info(f"Loading insertion data from: {input_file}")

    df = pd.read_csv(
        input_file,
        sep="\t",
        index_col=[0, 1, 2, 3],
        header=[0, 1],
    )
    logger.success(f"Loaded {df.shape[0]:,} insertions with {df.shape[1]} columns")

    if not isinstance(df.columns, pd.MultiIndex):
        raise ValueError("Expected MultiIndex columns with (Sample, Timepoint) structure")

    if len(df.columns.levels) != 2:
        raise ValueError("Expected columns to have 2 levels: Sample and Timepoint")

    return df


@logger.catch
def validate_timepoint_exists(df: pd.DataFrame, timepoint: str) -> None:
    """Validate that the specified timepoint exists in the data."""
    available_timepoints = df.columns.get_level_values(1).unique()
    if timepoint not in available_timepoints:
        logger.error(f"Timepoint '{timepoint}' not found in data")
        logger.error(f"Available timepoints: {list(available_timepoints)}")
        raise ValueError(f"Invalid timepoint: {timepoint}")

    logger.debug(f"Validated timepoint '{timepoint}' exists in data")


@logger.catch
def apply_hard_filtering(
    df: pd.DataFrame, initial_timepoint: str, cutoff_threshold: int
) -> tuple[pd.DataFrame, AnalysisResult]:
    """Apply hard filtering across all samples based on the read-count threshold."""
    logger.info("Starting hard filtering process...")

    # Validate timepoint exists
    validate_timepoint_exists(df, initial_timepoint)

    # Display initial data info
    logger.info("=" * 60)
    logger.info("INITIAL DATA SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total insertions: {df.shape[0]:,}")
    logger.info(f"Total samples: {len(df.columns.get_level_values(0).unique())}")
    logger.info(f"Timepoints: {list(df.columns.get_level_values(1).unique())}")
    logger.info(f"Initial timepoint: {initial_timepoint}")
    logger.info(f"Cutoff threshold: {cutoff_threshold}")

    # Process each sample
    filtered_samples = {}
    sample_results = []
    total_insertions = 0
    retained_insertions = 0

    for sample_name, sample_data_t in df.T.groupby(level="Sample"):
        sample_data = sample_data_t.T
        logger.debug(f"Processing sample: {sample_name}")

        sample_total = len(sample_data)
        mask = sample_data[(sample_name, initial_timepoint)] >= cutoff_threshold
        sample_retained = mask.sum()

        retention_rate = (sample_retained / sample_total * 100
                          if sample_total > 0 else 0)

        total_insertions += sample_total
        retained_insertions += sample_retained

        logger.debug(
            f"Sample {sample_name}: {sample_retained:,}/{sample_total:,} "
            f"insertions retained ({retention_rate:.2f}%)"
        )

        if sample_retained > 0:
            filtered_samples[sample_name] = sample_data[mask]

        sample_results.append({
            "sample_name": sample_name,
            "total_insertions": sample_total,
            "retained_insertions": sample_retained,
            "retention_rate": retention_rate,
        })

    # Combine filtered samples
    if not filtered_samples:
        logger.warning("No samples retained any insertions after filtering")
        filtered_df = pd.DataFrame(columns=df.columns)
    else:
        filtered_df = pd.concat(filtered_samples.values(), axis=1)

    # Calculate overall statistics
    removed_insertions = total_insertions - retained_insertions
    retention_rate = (retained_insertions / total_insertions * 100
                      if total_insertions > 0 else 0)

    result = AnalysisResult(
        total_insertions=total_insertions,
        retained_insertions=retained_insertions,
        removed_insertions=removed_insertions,
        retention_rate=retention_rate,
        samples_processed=len(sample_results),
    )

    return filtered_df, result
