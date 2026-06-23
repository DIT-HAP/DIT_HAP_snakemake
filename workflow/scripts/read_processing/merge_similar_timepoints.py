"""
Merge similar timepoints by summing read counts across specified timepoint columns.

Reads a wide-format insertion count TSV with timepoint columns, sums the specified
similar timepoints into a single merged column, drops the original columns, and
re-sorts columns alphabetically.

Typical Usage:
    python merge_similar_timepoints.py \
        -i counts.tsv -o merged.tsv \
        -s 0h YES0 -m YES0 -d 0h

Input: TSV with (Chr, Coordinate, Strand, Target) row MultiIndex, timepoint columns
Output: Same format with merged timepoint column replacing similar_timepoints
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Optional
from dataclasses import dataclass

from loguru import logger
import pandas as pd


# =============================== Configuration & Models ===============================

@dataclass
class Config:
    """Configuration for merging similar timepoints."""
    input_file: Path
    output_file: Path
    similar_timepoints: list
    merged_timepoint: str
    drop_columns: list

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
def merge_timepoints(config: Config) -> None:
    """Merge similar timepoints by summing read counts."""
    df = pd.read_csv(config.input_file, sep="\t", index_col=[0, 1, 2, 3])
    logger.info(f"Loaded {df.shape[0]:,} rows, columns: {list(df.columns)}")

    logger.info(f"Summing {config.similar_timepoints} → '{config.merged_timepoint}'")
    df[config.merged_timepoint] = df[config.similar_timepoints].sum(axis=1)

    logger.info(f"Dropping columns: {config.drop_columns}")
    df.drop(columns=config.drop_columns, inplace=True)
    df.sort_index(axis=1, inplace=True)

    df.to_csv(config.output_file, sep="\t", index=True)
    logger.success(f"Merged → {df.shape[0]:,} rows, columns: {list(df.columns)}")
    logger.success(f"Saved to {config.output_file}")


# =============================== Main Function ===============================

def parse_arguments() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge similar timepoints by summing read counts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python merge_similar_timepoints.py \\
      -i counts.tsv -o merged.tsv \\
      -s 0h YES0 -m YES0 -d 0h
        """
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Input counts TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output merged counts TSV file")
    parser.add_argument("-s", "--similar-timepoints", nargs="+", required=True,
                        help="Timepoint column names to sum together")
    parser.add_argument("-m", "--merged-timepoint", type=str, required=True,
                        help="Name for the merged timepoint column")
    parser.add_argument("-d", "--drop-columns", nargs="+", required=True,
                        help="Column names to drop after merging")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


@logger.catch
def main() -> None:
    """Main entry point of the script."""
    args = parse_arguments()
    setup_logging("DEBUG" if args.verbose else "INFO")

    config = Config(
        input_file=args.input,
        output_file=args.output,
        similar_timepoints=args.similar_timepoints,
        merged_timepoint=args.merged_timepoint,
        drop_columns=args.drop_columns,
    )

    logger.info(f"Merging timepoints in {config.input_file}")
    merge_timepoints(config)
    logger.success("Script completed successfully!")


if __name__ == "__main__":
    main()
