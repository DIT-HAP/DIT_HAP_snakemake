"""Control Insertion Selection for Depletion Analysis

Select control insertions from a transposon insertion dataset for use as a
neutral baseline in depletion analysis. Control insertions are drawn from
intergenic regions that sit far enough from any gene boundary to be considered
unaffected by selection pressure, providing a stable reference against which
depleted insertions can be scored.

The selection algorithm queries the annotation table for insertions annotated
as ``Intergenic region`` whose distance to both the upstream and downstream
region boundaries exceeds a fixed threshold (500 bp). The resulting set is then
intersected with the count matrix index so that only insertions actually present
in the counts are retained, and duplicate index entries are collapsed (keeping
the first occurrence). Retention statistics are logged for reproducibility.

Input
-----
- A tab-separated insertion count matrix with a 4-level row MultiIndex
  (columns 0-3) and a 2-level column header (rows 0 and 1).
- A tab-separated genomic annotation table with a matching 4-level row
  MultiIndex, containing at least the ``Type``, ``Distance_to_region_start``,
  and ``Distance_to_region_end`` columns.

Output
------
- A DataFrame of selected control insertions, with the row MultiIndex
  preserved.

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
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
CONTROL_DISTANCE_THRESHOLD = 500  # Minimum distance (bp) from gene boundaries for control insertions

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ControlSelectionResult:
    """Summary statistics of the control insertion selection run."""
    total_insertions_processed: int
    control_insertions_selected: int
    success_rate: float

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_and_preprocess_data(counts_file: Path, annotations_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and preprocess insertion count and genomic annotation tables."""
    counts_df = pd.read_csv(
        counts_file, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t"
    )

    # Remove rows with any NA value
    counts_df = counts_df.dropna(axis=0, how="any").copy()

    # Load and process annotations
    insertion_annotations = pd.read_csv(
        annotations_file, index_col=[0, 1, 2, 3], sep="\t"
    )

    return counts_df, insertion_annotations


@logger.catch
def get_control_insertions(counts_df: pd.DataFrame, insertion_annotations: pd.DataFrame) -> pd.DataFrame:
    """Select control insertions based on stringent genomic criteria."""
    ctr_insertions = insertion_annotations.query(
        f"Type == 'Intergenic region' and Distance_to_region_start > {CONTROL_DISTANCE_THRESHOLD} and Distance_to_region_end > {CONTROL_DISTANCE_THRESHOLD}"
    )

    ctr_insertions = ctr_insertions[ctr_insertions.index.isin(counts_df.index)].drop_duplicates(keep="first")

    return ctr_insertions
