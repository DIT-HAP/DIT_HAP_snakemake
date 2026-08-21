"""Read count distribution computation: histogram binning and cutoff statistics.

Extracted from ``workflow/scripts/quality_control/read_count_distribution_analysis.py``.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
DISTRIBUTION_COLUMNS = ["sample", "timepoint", "bin_left", "bin_right", "count"]
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
def load_and_validate_data(file_path: Path) -> pd.DataFrame:
    """Load a TSV file with its 4-level row index and reject an empty table."""
    df = pd.read_csv(file_path, sep="\t", engine="python", index_col=[0, 1, 2, 3])

    if df.empty:
        raise ValueError("Empty DataFrame after reading")

    return df


@logger.catch
def parse_sample_name(file_path: Path) -> str:
    """Derive the sample label from a filename by dropping every dotted suffix."""
    return file_path.name.split(".")[0]


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


@logger.catch
def compute_binned_distribution(df: pd.DataFrame, sample: str, bins: int) -> pd.DataFrame:
    """Bin log10 of strictly positive values for every numeric column into a long-format frame."""
    numeric_cols = df.select_dtypes(include=np.number).columns
    if numeric_cols.empty:
        raise ValueError("No numeric columns found")

    frames: list[pd.DataFrame] = []

    for col_name in numeric_cols:
        col_data = df[col_name].dropna()
        positive = col_data[col_data > 0]

        # Marker row keeps the group present in the output so the renderer can
        # still allocate a "No valid data" panel for it.
        if positive.empty:
            logger.warning(f"  {sample} / {col_name}: no positive values, emitting empty marker row")
            frames.append(
                pd.DataFrame(
                    {
                        "sample": [sample],
                        "timepoint": [str(col_name)],
                        "bin_left": [np.nan],
                        "bin_right": [np.nan],
                        "count": [0],
                    }
                )
            )
            continue

        counts, edges = np.histogram(np.log10(positive.to_numpy()), bins=bins)

        frames.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "timepoint": str(col_name),
                    "bin_left": edges[:-1],
                    "bin_right": edges[1:],
                    "count": counts,
                }
            )
        )
        logger.info(f"  {sample} / {col_name}: {len(positive)} positive values across {bins} bins")

    return pd.concat(frames, ignore_index=True)[DISTRIBUTION_COLUMNS]


@logger.catch
def log_summary_table(stats_df: pd.DataFrame) -> None:
    """Log the per-sample retention statistics as a plain-text table."""
    if stats_df.empty:
        logger.info("No statistics to display.")
        return

    logger.info("--- Processing Summary ---")
    for line in stats_df.to_string(index=False).split("\n"):
        logger.info(line)
