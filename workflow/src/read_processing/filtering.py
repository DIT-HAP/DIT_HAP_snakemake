"""Threshold-based filtering of aligned read pairs.

Filters aligned read pairs from BAM-derived tables using thresholds loaded
from a YAML configuration file. R1 and R2 reads are filtered independently
with separate MAPQ, NCIGAR, and NM thresholds, plus optional rejection of
reads that carry supplementary (SA) or secondary (XA) alignments and an
optional proper-pair requirement.

``filter_aligned_reads.py``'s ``filter_read_pairs`` is the chunk loop, reader,
and writer and stays in that script per Phase 4 convention 3.

Input
-----
- A chunk (DataFrame) of read-pair data with columns ``R1_MAPQ``, ``R2_MAPQ``,
  ``R1_NCIGAR``, ``R2_NCIGAR``, ``R1_NM``, ``R2_NM``, ``R1_SA``, ``R2_SA``,
  ``R1_XA``, ``R2_XA``, and ``Is_Proper_Pair``.
- YAML config with an ``aligned_read_filtering`` section containing
  ``read_1_filtering``, ``read_2_filtering``, and ``require_proper_pair``.

Output
------
- ``FilterThresholds`` and ``AnalysisResult`` dataclasses, a per-chunk boolean
  filter mask, and the filtered chunk.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import yaml
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
NA_VALUES = ["N/A", "NA", ""]

# SAM optional-field type codes: these tags are numeric ('i'/'f') per the SAM
# spec (MC/MD/SA/XA are strings and need no dtype coercion below).
NUMERIC_TAGS = frozenset({"AS", "MQ", "NM", "XS"})

# Base (non-tag) fields that are always populated by parse_bam_to_tsv.py and
# never carry "N/A" (mapq/n_cigar/flag default to 0, never "N/A" -- see
# ReadInfo).
ALWAYS_INT_FIELDS = frozenset({"MAPQ", "NCIGAR", "Flag", "FLAG"})

# Base fields plus numeric tags that CAN be "N/A" -- coerced to nullable float
# to match the dtype pd.read_csv would infer for a column containing NaN.
NULLABLE_NUMERIC_FIELDS = frozenset({"LEN", "Pos", "Ref_Start", "Ref_End"}) | NUMERIC_TAGS

# Columns checked via .isna() in build_filter_mask -- need the "N/A"/"NA"/""
# string sentinels turned into real NaN. Every other string column is left
# untouched: pd.to_numeric(errors="coerce") already treats those sentinels as
# unparseable on its own, and no other string column is null-checked.
ISNA_CHECKED_COLUMNS = ("R1_SA", "R1_XA", "R2_SA", "R2_XA")

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FilterThresholds:
    """Per-read filtering thresholds for MAPQ, NCIGAR, NM, and SA/XA alignments."""
    mapq_threshold: float | None = None
    ncigar_value: int | None = None
    nm_threshold: int | None = None
    no_sa: bool = False
    no_xa: bool = False


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Statistics from a filtering run."""
    total_rows: int
    filtered_rows: int
    removed_rows: int
    retention_rate: float
    chunks_processed: int

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
def load_config_from_yaml(config_file) -> dict[str, Any]:
    """Load filtering thresholds from the ``aligned_read_filtering`` section of a YAML file."""
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Extract aligned_read_filtering configuration
    if "aligned_read_filtering" not in config:
        raise ValueError("'aligned_read_filtering' section not found in config file")

    filtering_config = config["aligned_read_filtering"]

    # Convert YAML config to internal format
    internal_config: dict[str, Any] = {}
    for read in ["read_1_filtering", "read_2_filtering"]:
        internal_config[read] = FilterThresholds(
            mapq_threshold=filtering_config.get(read, {}).get("mapq_threshold"),
            ncigar_value=filtering_config.get(read, {}).get("ncigar_value"),
            nm_threshold=filtering_config.get(read, {}).get("nm_threshold"),
            no_sa=filtering_config.get(read, {}).get("no_sa"),
            no_xa=filtering_config.get(read, {}).get("no_xa"),
        )

    internal_config["require_proper_pair"] = filtering_config.get("require_proper_pair")

    logger.info(f"Loaded configuration from: {config_file}")
    logger.debug(f"R1 filters: MAPQ={internal_config['read_1_filtering'].mapq_threshold}, NCIGAR={internal_config['read_1_filtering'].ncigar_value}, NM={internal_config['read_1_filtering'].nm_threshold}")
    logger.debug(f"R2 filters: MAPQ={internal_config['read_2_filtering'].mapq_threshold}, NCIGAR={internal_config['read_2_filtering'].ncigar_value}, NM={internal_config['read_2_filtering'].nm_threshold}")
    logger.debug(f"Pair filters: no_sa={internal_config['read_1_filtering'].no_sa}, no_xa={internal_config['read_1_filtering'].no_xa}, proper_pair={internal_config['require_proper_pair']}")

    return internal_config


def strip_read_prefix(column: str) -> str:
    """Strip a leading R1_/R2_ prefix from a column name, if present."""
    if column.startswith("R1_") or column.startswith("R2_"):
        return column[3:]
    return column


def coerce_column_dtypes(chunk: pd.DataFrame) -> pd.DataFrame:
    """Coerce an all-string Parquet chunk to the dtypes pd.read_csv would infer from the equivalent TSV."""
    for column in ISNA_CHECKED_COLUMNS:
        if column in chunk.columns:
            chunk[column] = chunk[column].replace(NA_VALUES, np.nan)
    for column in chunk.columns:
        field = strip_read_prefix(column)
        if field in ALWAYS_INT_FIELDS:
            chunk[column] = chunk[column].astype("int64")
        elif field in NULLABLE_NUMERIC_FIELDS:
            chunk[column] = pd.to_numeric(chunk[column], errors="coerce").astype("float64")
    return chunk


@logger.catch
def build_filter_mask(
    chunk: pd.DataFrame,
    r1_filters: FilterThresholds,
    r2_filters: FilterThresholds,
    require_proper_pair: bool,
) -> pd.Series:
    """Build a boolean mask selecting read pairs that pass every configured criterion."""
    # Initialize mask with all True
    filter_mask = pd.Series([True] * len(chunk), index=chunk.index)

    # Apply R1 filters
    if r1_filters.mapq_threshold is not None:
        filter_mask &= (chunk["R1_MAPQ"] >= r1_filters.mapq_threshold)

    if r1_filters.ncigar_value is not None:
        filter_mask &= (chunk["R1_NCIGAR"] <= r1_filters.ncigar_value)

    if r1_filters.nm_threshold is not None:
        filter_mask &= (chunk["R1_NM"] <= r1_filters.nm_threshold)

    if r1_filters.no_sa:
        filter_mask &= (chunk["R1_SA"].isna() | (chunk["R1_SA"] == "N/A"))

    if r1_filters.no_xa:
        filter_mask &= (chunk["R1_XA"].isna() | (chunk["R1_XA"] == "N/A"))

    # Apply R2 filters
    if r2_filters.mapq_threshold is not None:
        filter_mask &= (chunk["R2_MAPQ"] >= r2_filters.mapq_threshold)

    if r2_filters.ncigar_value is not None:
        filter_mask &= (chunk["R2_NCIGAR"] <= r2_filters.ncigar_value)

    if r2_filters.nm_threshold is not None:
        filter_mask &= (chunk["R2_NM"] <= r2_filters.nm_threshold)

    if r2_filters.no_sa:
        filter_mask &= (chunk["R2_SA"].isna() | (chunk["R2_SA"] == "N/A"))

    if r2_filters.no_xa:
        filter_mask &= (chunk["R2_XA"].isna() | (chunk["R2_XA"] == "N/A"))

    # Apply proper pair filter
    if require_proper_pair:
        filter_mask &= (chunk["Is_Proper_Pair"].str.capitalize() == "Yes")

    return filter_mask


@logger.catch
def process_chunk(
    chunk: pd.DataFrame,
    chunk_num: int,
    r1_filters: FilterThresholds,
    r2_filters: FilterThresholds,
    require_proper_pair: bool,
    first_chunk: bool,
) -> tuple[pd.DataFrame, bool]:
    """Filter a single chunk and log first-chunk diagnostics and periodic progress."""
    chunk = coerce_column_dtypes(chunk)
    chunk_rows_before = len(chunk)

    # Display info for first chunk
    if first_chunk:
        logger.info("=" * 60)
        logger.info("ORIGINAL DATA INFORMATION")
        logger.info("=" * 60)
        logger.info(f"Columns: {len(chunk.columns)}")
        logger.info(f"First chunk size: {chunk_rows_before:,} rows")

        logger.debug("Column Data Types:")
        for col, dtype in chunk.dtypes.items():
            logger.debug(f"  {col}: {dtype}")

    # Build and apply filter mask
    filter_mask = build_filter_mask(
        chunk,
        r1_filters,
        r2_filters,
        require_proper_pair,
    )

    filtered_chunk = chunk[filter_mask]
    chunk_filtered_rows = len(filtered_chunk)

    # Log progress
    if chunk_num == 1 or chunk_num % 10 == 0:
        retention_rate = (chunk_filtered_rows / chunk_rows_before * 100
                          if chunk_rows_before > 0 else 0)
        logger.info(
            f"Chunk {chunk_num}: {chunk_filtered_rows:,}/{chunk_rows_before:,} "
            f"rows retained ({retention_rate:.1f}%)"
        )

    return filtered_chunk, False
