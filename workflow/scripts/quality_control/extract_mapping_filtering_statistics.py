#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
#     "loguru",
# ]
# ///

"""
Extract Read-Pair Filtering Statistics from Log Files
=====================================================

Parses log files produced during read-pair filtering to extract per-sample
filtering statistics for PBL (left) and PBR (right) read pairs. For each log
it reads the "FILTERING SUMMARY" blocks, pulls out chunk counts, original /
filtered / removed read-pair counts and retention rates, then aggregates the
PBL and PBR figures into combined totals and an overall retention rate.

The parser matches summary blocks of the following shape (one per read type)::

    ============================================================
    FILTERING SUMMARY
    ============================================================
    Total chunks processed: 12345
    Original read pairs: 1,234,567
    Filtered read pairs: 1,123,456
    Removed read pairs: 111,111
    Overall retention rate: 91.05%
    Output written to: sample_name.PBL.filtered.parquet

Input
-----
- One or more filtering log files, each containing "FILTERING SUMMARY"
  sections whose "Output written to" line ends in ``.PBL.filtered.parquet`` or
  ``.PBR.filtered.parquet``. The sample name is taken from each log file's stem.

Output
------
- A tab-separated statistics table indexed by sample (index label ``Sample``)
  with per-read-type counts and retention rates plus the aggregated columns
  ``total_original_pairs``, ``total_filtered_pairs`` and
  ``overall_retention_rate``. Floats are written with ``%.2f`` precision.

Usage
-----
    python extract_mapping_filtering_statistics.py -i sample1.log sample2.log -o filtering_statistics.tsv
    python extract_mapping_filtering_statistics.py -i sample1.log -o stats.tsv -v

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import re
import sys
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
class InputOutputConfig:
    """Validated input log paths and the output statistics path."""
    input_files: list[Path]
    output_file: Path

    def __post_init__(self) -> None:
        """Validate that every input file exists and create the output parent directory."""
        for file_path in self.input_files:
            if not file_path.exists():
                raise ValueError(f"Input file does not exist: {file_path}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


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
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


setup_logger()

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

    # Calculate overall retention rate
    if "total_original_pairs" in df.columns and "total_filtered_pairs" in df.columns:
        df["overall_retention_rate"] = (
            df["total_filtered_pairs"] / df["total_original_pairs"]
        ).round(4)

    return df


@logger.catch
def save_results(df: pd.DataFrame, output_file: Path) -> None:
    """Save the filtering statistics DataFrame to a TSV file."""
    if df.empty:
        logger.error("Cannot save empty DataFrame")
        return

    # Sort by sample name
    df = df.rename_axis("Sample", axis=0).sort_index()

    # Save to file
    df.to_csv(output_file, sep="\t", index=True, float_format="%.2f")
    logger.success(f"Statistics saved to: {output_file}")

    # Display summary
    summary_cols = ["total_original_pairs", "total_filtered_pairs", "overall_retention_rate"]
    available_cols = [col for col in summary_cols if col in df.columns]

    if available_cols:
        logger.info("Summary statistics:")
        logger.info(f"\n{df[available_cols].describe()}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Extract filtering statistics from read pair filtering log files"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        nargs="+",
        required=True,
        help="Path to the log files",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Path to save the output statistics TSV",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate inputs, extract statistics, and write the TSV."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        # Validate inputs
        config = InputOutputConfig(
            input_files=args.input,
            output_file=args.output,
        )

        logger.info("Starting filtering statistics extraction")

        # Extract statistics
        statistics = extract_summary_data(config.input_files)

        if not statistics:
            logger.error("No statistics extracted from any log files")
            return 1

        # Create DataFrame
        df = create_dataframe(statistics)

        # Save results
        save_results(df, config.output_file)

        # Create analysis result
        result = AnalysisResult(
            total_samples_processed=len(statistics),
            total_log_files=len(config.input_files),
            output_path=config.output_file,
        )

        logger.success("Extraction completed successfully!")
        logger.info(f"Processed {result.total_samples_processed} samples from {result.total_log_files} log files")
        logger.info(f"Output saved to: {result.output_path}")

    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
