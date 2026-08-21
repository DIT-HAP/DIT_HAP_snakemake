"""Insertion-level depletion analysis for non-replicated samples.

Perform insertion-level depletion analysis on transposon sequencing data that
lacks replicates. Loads insertion counts and a set of control insertions,
normalises counts against the control-insertion medians, and computes
log-fold changes (LFC / M values) relative to an initial timepoint.

Input
-----
- Counts TSV: 4 index columns and a two-level column header (sample, timepoint).
- Control insertions TSV: 4 index columns.
- Initial timepoint label.

Output
------
- ``AnalysisResult`` dataclass and the LFC/normalized-counts/baseMean
  DataFrames, as prepared by the invoking script.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Results container for insertion-level depletion analysis."""
    status: str
    message: str


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_and_preprocess_data(
    counts_file: Path, control_insertions_file: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load counts and control insertions data from TSV files."""
    counts_df = pd.read_csv(
        counts_file, sep="\t", index_col=[0, 1, 2, 3], header=[0, 1]
    )
    control_insertions_df = pd.read_csv(
        control_insertions_file, sep="\t", index_col=[0, 1, 2, 3]
    )
    return counts_df, control_insertions_df


@logger.catch
def perform_median_normalization(
    counts_df: pd.DataFrame, control_insertions_df: pd.DataFrame
) -> pd.DataFrame:
    """Normalize counts using median values from control insertions."""
    median_values = counts_df.loc[control_insertions_df.index].median()
    min_median_values = median_values.min()
    normalized_counts = counts_df.mul(min_median_values).div(median_values)
    return normalized_counts


@logger.catch
def calculate_MA_values(
    normalized_counts: pd.DataFrame, init_timepoint: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate M (log-fold change) and A (average abundance) values for MA plots."""
    M_values = (
        -(normalized_counts + 1)
        .div((normalized_counts.xs(init_timepoint, level=1, axis=1) + 1), axis=0)
        .map(np.log2)
    )

    A_values = (normalized_counts + 1).mul(
        (normalized_counts.xs(init_timepoint, level=1, axis=1) + 1), axis=0
    ).map(np.log2) * 0.5

    return M_values, A_values
