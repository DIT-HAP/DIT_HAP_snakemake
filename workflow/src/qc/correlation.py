"""PBL-PBR pairs extraction: TSV loading and filename parsing.

Extracted from ``workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py``.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pandas as pd
from loguru import logger


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def read_tsv_file(file_path: Path) -> pd.DataFrame | None:
    """Read a TSV file and return Chr/Coordinate/Strand plus strictly-positive PBL/PBR pairs, or None if invalid."""
    logger.info(f"Reading TSV file: {file_path}")

    df = pd.read_csv(file_path, sep='\t', index_col=[0, 1, 2])

    # Check if PBL and PBR columns exist
    required_cols = ['PBL', 'PBR']
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        logger.warning(f"Warning: Missing columns {missing_cols} in {file_path}")
        return None

    # Remove rows with missing values in PBL or PBR
    df_clean = df[['PBL', 'PBR']].dropna()

    # Remove zero or negative values for log scaling
    df_clean = df_clean[(df_clean['PBL'] > 0) & (df_clean['PBR'] > 0)]

    if df_clean.empty:
        logger.warning(f"Warning: No valid data points in {file_path}")
        return None

    # Restore Chr/Coordinate/Strand from the row index as plain columns
    return df_clean.reset_index()


@logger.catch
def parse_filename(file_path: Path) -> tuple[str, str, str] | None:
    """Parse filename stem into sample, timepoint, condition. Return None if format is invalid."""
    stem = file_path.stem
    parts = stem.split('_')

    if len(parts) != 3:
        logger.warning(f"Filename {file_path.name} does not match expected pattern {{sample}}_{{timepoint}}_{{condition}}.tsv (got {len(parts)} parts)")
        return None

    sample, timepoint, condition = parts
    return sample, timepoint, condition
