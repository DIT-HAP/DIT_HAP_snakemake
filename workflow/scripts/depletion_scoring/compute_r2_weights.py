#!/usr/bin/env python3

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
Compute R²-Based Confidence Weights from Curve Fitting
======================================================

Reads insertion-level curve fitting statistics and converts each insertion's R²
goodness-of-fit into a confidence score. R² values are clipped to the open
interval (1e-6, 1 - 1e-6) to avoid degenerate 0/1 weights, and the resulting
confidence is turned into a per-timepoint weight of ``1 - confidence``.

The output weight matrix is consumed by gene-level depletion analysis when
biological replicates are unavailable, letting poorly-fit insertions contribute
less to downstream aggregation.

Input
-----
- TSV with a (Chr, Coordinate, Strand, Target) MultiIndex, an ``R2`` column, and
  one ``*_fitted`` column per timepoint.

Output
------
- TSV with the same MultiIndex and one column per timepoint, each holding
  ``1 - confidence`` (tab-separated).

Usage
-----
    python compute_r2_weights.py \
        -i insertion_level_fitting_statistics.tsv \
        -o insertions_LFC_fitted_with_r_square_as_weights.tsv
    python compute_r2_weights.py -i stats.tsv -o weights.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
R2_CLIP_MIN = 1e-6
R2_CLIP_MAX = 1 - 1e-6

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class Config:
    """Configuration for computing R² weights."""
    input_file: Path
    output_file: Path

    def __post_init__(self) -> None:
        """Validate the input exists and create the output directory."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

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

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def compute_weights(config: Config) -> None:
    """Compute per-insertion R²-based confidence weights from curve fitting."""
    fitting_res = pd.read_csv(config.input_file, sep="\t", index_col=[0, 1, 2, 3])
    logger.info(f"Loaded fitting results: {fitting_res.shape[0]:,} insertions")

    fitting_res["confidence"] = fitting_res["R2"].clip(lower=R2_CLIP_MIN, upper=R2_CLIP_MAX)

    # Extract timepoint names by stripping the "_fitted" suffix from column names.
    # Uses rstrip (character set, not suffix) — correct for timepoints like YES0/YES24
    # that contain no characters from the set {_, f, i, t, e, d}.
    tp_cols = [col.rstrip("_fitted")
               for col in fitting_res.filter(regex=r".*_fitted$").columns.tolist()]
    logger.info(f"Timepoints detected: {tp_cols}")

    weights = pd.DataFrame(index=fitting_res.index, columns=tp_cols)
    for col in tp_cols:
        weights[col] = 1 - fitting_res["confidence"]

    weights.to_csv(config.output_file, sep="\t")
    logger.success(f"Weights for {len(tp_cols):,} timepoints → {config.output_file}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Compute R²-based confidence weights from insertion-level curve fitting",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compute_r2_weights.py \\
      -i insertion_level_fitting_statistics.tsv \\
      -o insertions_LFC_fitted_with_r_square_as_weights.tsv
        """,
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Insertion-level curve fitting statistics TSV")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output weights TSV")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()

def main() -> int:
    """Main orchestrator: parse args, build config, compute weights."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = Config(input_file=args.input, output_file=args.output)
        logger.info(f"Computing R² weights from {config.input_file}")
        compute_weights(config)
        logger.success("Script completed successfully!")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
