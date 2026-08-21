"""Insertion orientation pairs extraction: strand-pair computation.

Extracted from ``workflow/scripts/quality_control/insertion_orientation_analysis.py``.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
OUTPUT_COLUMNS = ["sample", "timepoint", "plus_count", "minus_count"]


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def extract_strand_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Pair +/- strand counts per Sample x Timepoint and keep rows with both strands positive."""
    required_levels = {"Strand", "Sample", "Timepoint"}
    available_levels = set(df.index.names) | set(df.columns.names)
    missing_levels = required_levels - available_levels
    if missing_levels:
        raise ValueError(f"Missing required index/column levels: {sorted(missing_levels)}")

    # Stack both column levels then pivot Strand back out, so each row holds the
    # +/- pair. dropna(axis=0) removes keys lacking one strand, exactly as before.
    plus_minus_pair = (
        df.stack(future_stack=True).stack(future_stack=True).unstack("Strand").dropna(axis=0)
    )
    logger.info(f"Paired {len(plus_minus_pair)} strand rows before positivity filtering")

    missing_strands = [strand for strand in ("+", "-") if strand not in plus_minus_pair.columns]
    if missing_strands:
        raise ValueError(f"Missing strand columns after unstacking: {missing_strands}")

    # Preserve the original min(axis=1) > 0 semantics: both strands strictly
    # positive, which is what makes the log-log correlation well defined.
    filtered = plus_minus_pair[plus_minus_pair.min(axis=1) > 0]
    logger.info(f"Retained {len(filtered)} rows with both strands strictly positive")

    if filtered.empty:
        logger.warning("No valid strand pairs after positivity filtering!")
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    pairs = (
        filtered[["+", "-"]]
        .rename(columns={"+": "plus_count", "-": "minus_count"})
        .reset_index()
        .rename(columns={"Sample": "sample", "Timepoint": "timepoint"})
    )

    return pairs[OUTPUT_COLUMNS]
