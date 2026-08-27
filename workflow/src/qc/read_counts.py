"""Read count cutoff retention statistics.

Once held the full pre-binned distribution pipeline (binning, loading,
sample-name parsing); those parts died when the read count figure became a
single-stage script that histograms raw values directly. Only the cutoff
retention statistic survives — the footer annotations of that figure still
consume it.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
from loguru import logger


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
STATS_COLUMNS = [
    "sample",
    "original_rows",
    "rows_kept",
    "pct_rows_kept",
    "original_counts",
    "counts_kept",
    "pct_counts_kept",
]


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def calculate_cutoff_statistics(df: pd.DataFrame, initial_time_point: str, cutoff: float) -> dict[str, float | int]:
    """Compute row and count retention after applying the cutoff to the initial time point."""
    if initial_time_point not in df.columns:
        raise ValueError(f"Initial time point column '{initial_time_point}' not found in {list(df.columns)}")

    if not pd.api.types.is_numeric_dtype(df[initial_time_point]):
        raise ValueError(f"Initial time point column '{initial_time_point}' is not numeric")

    original_rows = len(df)
    original_counts = float(df[initial_time_point].sum())

    kept = df.loc[df[initial_time_point] >= cutoff, initial_time_point]
    rows_kept = len(kept)
    counts_kept = float(kept.sum())

    return {
        "original_rows": original_rows,
        "rows_kept": rows_kept,
        "pct_rows_kept": (rows_kept / original_rows) * 100.0 if original_rows > 0 else 0.0,
        "original_counts": original_counts,
        "counts_kept": counts_kept,
        "pct_counts_kept": (counts_kept / original_counts) * 100.0 if original_counts > 0 else 0.0,
    }
