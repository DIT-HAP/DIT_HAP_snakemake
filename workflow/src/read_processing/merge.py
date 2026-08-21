"""Merge operations on insertion tables.

Holds merge logic for two related but independent scripts: merging PBL/PBR
strand-specific insertion counts into a single strand-resolved table
(``merge_strand_insertions.py``), and merging similar timepoint columns in a
wide count matrix (``merge_similar_timepoints.py``). Both operate on insertion
tables and their functions do not collide, so they share one module.

``merge_similar_timepoints.py``'s only function, ``merge_timepoints``, is pure
orchestration plus file I/O with no domain logic to extract (per Phase 4
convention 2), so it dissolves entirely into that script's ``main()`` and
contributes nothing here.

Input
-----
- PBL and PBR insertion DataFrames, as prepared by the invoking script.

Output
------
- ``MergeResult`` dataclass and the merged insertion DataFrame.

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
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class MergeResult:
    """Summary statistics produced by the merge operation."""
    total_sites_processed: int
    total_reads_merged: int
    coordinate_strand_pairs: int

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def merge_insertion_data(pbl_df: pd.DataFrame, pbr_df: pd.DataFrame) -> pd.DataFrame:
    """Merge PBL and PBR insertion data into a strand-resolved table."""
    logger.info("Merging PBL and PBR insertion data")

    merged_df = pd.merge(
        pbl_df,
        pbr_df,
        how="outer",
        on=["Chr", "Coordinate"],
        suffixes=("_PBL", "_PBR"),
    ).fillna(0)

    # Create plus strand data
    plus_df = merged_df[["Chr", "Coordinate", "-_PBL", "+_PBR"]].copy()
    plus_df["Strand"] = "+"
    plus_df.rename(columns={"-_PBL": "PBL", "+_PBR": "PBR"}, inplace=True)

    # Create minus strand data
    minus_df = merged_df[["Chr", "Coordinate", "+_PBL", "-_PBR"]].copy()
    minus_df["Strand"] = "-"
    minus_df.rename(columns={"+_PBL": "PBL", "-_PBR": "PBR"}, inplace=True)

    # Combine and finalize
    final_df = pd.concat([plus_df, minus_df], axis=0)
    final_df = final_df.set_index(["Chr", "Coordinate", "Strand"])
    final_df = final_df.astype(int).sort_index()
    final_df["Reads"] = final_df["PBL"] + final_df["PBR"]

    logger.success(f"Merged data: {len(final_df):,} coordinate-strand pairs")
    return final_df
