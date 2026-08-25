"""Pure computation for mapping/filtering statistics extraction.

Moved out of
``workflow/scripts/quality_control/extract_mapping_filtering_statistics.py``
so the log-parsing and aggregation logic (PBL/PBR filtering summaries) can be
imported and unit tested independently of the CLI entrypoint.
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import re
from dataclasses import asdict, dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
SUMMARY_PATTERN = re.compile(
    r".*\| ============================================================\s*\n"
    r".*\| FILTERING SUMMARY\s*\n"
    r".*\| ============================================================\s*\n"
    r".*\| Total chunks processed: (\d+)\s*\n"
    r".*\| Original read pairs: ([\d,]+)\s*\n"
    r".*\| Filtered read pairs: ([\d,]+)\s*\n"
    r".*\| Removed read pairs: ([\d,]+)\s*\n"
    r".*\| Overall retention rate: ([\d.]+)%\s*\n"
    r".*\| Output written to: (.+?\.(?:PBL|PBR)\.filtered\.parquet)"
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FilteringStatistics:
    """Per-sample filtering counts and retention rates for PBL and PBR read pairs."""
    chunks_processed_pbl: int | None = None
    original_read_pairs_pbl: int | None = None
    filtered_read_pairs_pbl: int | None = None
    removed_read_pairs_pbl: int | None = None
    retention_rate_pbl: float | None = None
    output_file_pbl: str | None = None
    chunks_processed_pbr: int | None = None
    original_read_pairs_pbr: int | None = None
    filtered_read_pairs_pbr: int | None = None
    removed_read_pairs_pbr: int | None = None
    retention_rate_pbr: float | None = None
    output_file_pbr: str | None = None
    total_original_pairs: int | None = None
    total_filtered_pairs: int | None = None
    overall_retention_rate: float | None = None


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Summary of how many samples and log files were processed."""
    total_samples_processed: int
    total_log_files: int
    output_path: Path


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def parse_log_file(log_file: Path) -> dict[str, FilteringStatistics]:
    """Parse a single log file and extract filtering statistics."""
    sample_name = log_file.stem
    logger.info(f"Processing: {sample_name}")

    try:
        with open(log_file, "r") as f:
            content = f.read()
    except Exception as e:
        logger.error(f"Error reading {log_file}: {str(e)}")
        return {}

    matches = SUMMARY_PATTERN.findall(content)

    if not matches:
        logger.warning(f"No filtering summary sections found in: {sample_name}")
        return {}

    logger.debug(f"Found {len(matches)} filtering summary sections in {sample_name}")

    stats_dict = {}

    for match in matches:
        chunks_processed = int(match[0])
        original_pairs = int(match[1].replace(",", ""))
        filtered_pairs = int(match[2].replace(",", ""))
        removed_pairs = int(match[3].replace(",", ""))
        retention_rate = float(match[4]) / 100
        output_path = match[5]

        # Determine if this is PBL or PBR based on output path
        if ".PBL.filtered.parquet" in output_path:
            suffix = "pbl"
        elif ".PBR.filtered.parquet" in output_path:
            suffix = "pbr"
        else:
            logger.warning(f"Could not determine PBL/PBR from output path: {output_path}")
            continue

        # Update statistics
        stats_dict.update({
            f"chunks_processed_{suffix}": chunks_processed,
            f"original_read_pairs_{suffix}": original_pairs,
            f"filtered_read_pairs_{suffix}": filtered_pairs,
            f"removed_read_pairs_{suffix}": removed_pairs,
            f"retention_rate_{suffix}": retention_rate,
            f"output_file_{suffix}": output_path,
        })

        logger.debug(f"  {suffix.upper()}: {original_pairs:,} -> {filtered_pairs:,} ({retention_rate*100:.2f}% retained)")

    return {sample_name: FilteringStatistics(**stats_dict)}


@logger.catch
def extract_summary_data(log_files: list[Path]) -> dict[str, FilteringStatistics]:
    """Extract filtering statistics from multiple log files."""
    logger.info(f"Found {len(log_files)} log files with filtering statistics")

    all_statistics = {}

    for log_file in log_files:
        file_stats = parse_log_file(log_file)
        all_statistics.update(file_stats)

    return all_statistics


@logger.catch
def create_dataframe(statistics: dict[str, FilteringStatistics]) -> pd.DataFrame:
    """Create a pandas DataFrame from filtering statistics dictionary."""
    if not statistics:
        logger.error("No statistics extracted from any log files")
        return pd.DataFrame()

    # Convert to DataFrame (exclude unset/None fields, matching model_dump(exclude_none=True))
    df = pd.DataFrame.from_dict(
        {
            sample: {k: v for k, v in asdict(stats).items() if v is not None}
            for sample, stats in statistics.items()
        },
        orient="index",
    )

    # Sort columns for better readability
    pbl_cols = [col for col in df.columns if col.endswith("_pbl")]
    pbr_cols = [col for col in df.columns if col.endswith("_pbr")]
    all_cols = sorted(pbl_cols) + sorted(pbr_cols)

    # Ensure all expected columns exist
    for col in all_cols:
        if col not in df.columns:
            df[col] = None

    df = df.reindex(columns=all_cols)

    # Calculate totals
    if "original_read_pairs_pbl" in df.columns and "original_read_pairs_pbr" in df.columns:
        df["total_original_pairs"] = df[["original_read_pairs_pbl", "original_read_pairs_pbr"]].sum(axis=1)

    if "filtered_read_pairs_pbl" in df.columns and "filtered_read_pairs_pbr" in df.columns:
        df["total_filtered_pairs"] = df[["filtered_read_pairs_pbl", "filtered_read_pairs_pbr"]].sum(axis=1)
        df["pbl_pbr_ratio"] = df["filtered_read_pairs_pbl"] / df["filtered_read_pairs_pbr"].replace(0, pd.NA)

    # Calculate overall retention rate
    if "total_original_pairs" in df.columns and "total_filtered_pairs" in df.columns:
        df["overall_retention_rate"] = (
            df["total_filtered_pairs"] / df["total_original_pairs"]
        ).round(4)

    return df
