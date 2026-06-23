"""
Compute per-insertion R²-based confidence weights from curve fitting results.

Reads insertion-level curve fitting statistics, clips R² values to (1e-6, 1-1e-6)
to define a confidence score, then outputs a weight matrix where each weight is
1 - confidence. These weights are used by gene-level depletion analysis when
biological replicates are unavailable.

Typical Usage:
    python compute_r2_weights.py \
        -i insertion_level_fitting_statistics.tsv \
        -o insertions_LFC_fitted_with_r_square_as_weights.tsv

Input: TSV with (Chr, Coordinate, Strand, Target) MultiIndex and R2 column plus
       *_fitted columns for each timepoint
Output: TSV with same index and one column per timepoint containing (1 - confidence)
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Optional
from dataclasses import dataclass

from loguru import logger
import pandas as pd


# =============================== Constants ===============================
R2_CLIP_MIN = 1e-6
R2_CLIP_MAX = 1 - 1e-6


# =============================== Configuration & Models ===============================

@dataclass
class Config:
    """Configuration for computing R² weights."""
    input_file: Path
    output_file: Path

    def __post_init__(self):
        """Validate input exists and create output dir."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)


# =============================== Setup Logging ===============================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# =============================== Core Functions ===============================

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


# =============================== Main Function ===============================

def parse_arguments() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute R²-based confidence weights from insertion-level curve fitting",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compute_r2_weights.py \\
      -i insertion_level_fitting_statistics.tsv \\
      -o insertions_LFC_fitted_with_r_square_as_weights.tsv
        """
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Insertion-level curve fitting statistics TSV")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output weights TSV")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


@logger.catch
def main() -> None:
    """Main entry point of the script."""
    args = parse_arguments()
    setup_logging("DEBUG" if args.verbose else "INFO")

    config = Config(input_file=args.input, output_file=args.output)

    logger.info(f"Computing R² weights from {config.input_file}")
    compute_weights(config)
    logger.success("Script completed successfully!")


if __name__ == "__main__":
    main()
