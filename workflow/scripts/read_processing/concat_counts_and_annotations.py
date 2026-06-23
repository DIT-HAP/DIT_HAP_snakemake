"""
Concatenate per-sample insertion count and annotation TSV files into single matrices.

Reads multiple per-sample count TSVs (one per sample/condition combination) and
multiple annotation TSVs, concatenates them into a wide count matrix with a
(Sample, Timepoint) MultiIndex on columns and a single deduplicated annotation table.

Typical Usage:
    python concat_counts_and_annotations.py \
        -i s1.tsv s2.tsv -a a1.tsv a2.tsv \
        -oc raw_reads.tsv -oa annotations.tsv

Input: Per-sample TSVs with (Chr, Coordinate, Strand, Target) row MultiIndex
Output: Combined counts TSV (wide) and deduplicated annotations TSV
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Optional
from dataclasses import dataclass, field

from loguru import logger
import pandas as pd


# =============================== Configuration & Models ===============================

@dataclass
class Config:
    """Configuration for concatenating per-sample counts and annotations."""
    counts_files: list
    annotation_files: list
    output_counts: Path
    output_annotations: Path

    def __post_init__(self):
        """Validate all inputs exist and create output dirs."""
        for f in [Path(x) for x in self.counts_files + self.annotation_files]:
            if not f.exists():
                raise ValueError(f"Input file does not exist: {f}")
        self.output_counts.parent.mkdir(parents=True, exist_ok=True)
        self.output_annotations.parent.mkdir(parents=True, exist_ok=True)


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
def concatenate(config: Config) -> None:
    """Concatenate per-sample count and annotation TSV files."""
    counts_df = {}
    for f in config.counts_files:
        key = Path(f).name.split(".")[0]
        counts_df[key] = pd.read_csv(f, sep="\t", index_col=[0, 1, 2, 3])
        logger.info(f"Loaded counts {key}: {counts_df[key].shape[0]:,} rows")
    counts = pd.concat(counts_df, axis=1).rename_axis(["Sample", "Timepoint"], axis=1)
    logger.info(f"Combined counts shape: {counts.shape[0]:,} rows × {counts.shape[1]} columns")

    annotations_df = {}
    for f in config.annotation_files:
        key = Path(f).name.split(".")[0]
        annotations_df[key] = pd.read_csv(f, index_col=[0, 1, 2, 3], sep="\t")
        logger.info(f"Loaded annotations {key}: {annotations_df[key].shape[0]:,} rows")
    annotations = (
        pd.concat(list(annotations_df.values()), axis=0)
        .reset_index()
        .drop_duplicates()
        .set_index(["Chr", "Coordinate", "Strand", "Target"])
    )

    counts.to_csv(config.output_counts, sep="\t", index=True)
    annotations.to_csv(config.output_annotations, sep="\t", index=True)
    logger.success(f"Counts: {counts.shape[0]:,} rows → {config.output_counts}")
    logger.success(f"Annotations: {annotations.shape[0]:,} rows → {config.output_annotations}")


# =============================== Main Function ===============================

def parse_arguments() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Concatenate per-sample insertion counts and annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python concat_counts_and_annotations.py \\
      -i s1_PBL.tsv s1_EMM.tsv -a s1_PBL.annotated.tsv s1_EMM.annotated.tsv \\
      -oc raw_reads.tsv -oa annotations.tsv
        """
    )
    parser.add_argument("-i", "--counts", nargs="+", type=Path, required=True,
                        help="Per-sample count TSV files")
    parser.add_argument("-a", "--annotations", nargs="+", type=Path, required=True,
                        help="Per-sample annotation TSV files")
    parser.add_argument("-oc", "--output-counts", type=Path, required=True,
                        help="Output combined counts TSV")
    parser.add_argument("-oa", "--output-annotations", type=Path, required=True,
                        help="Output combined annotations TSV")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


@logger.catch
def main() -> None:
    """Main entry point of the script."""
    args = parse_arguments()
    setup_logging("DEBUG" if args.verbose else "INFO")

    config = Config(
        counts_files=[str(f) for f in args.counts],
        annotation_files=[str(f) for f in args.annotations],
        output_counts=args.output_counts,
        output_annotations=args.output_annotations,
    )

    logger.info(f"Concatenating {len(config.counts_files):,} count files and "
                f"{len(config.annotation_files):,} annotation files")
    concatenate(config)
    logger.success("Script completed successfully!")


if __name__ == "__main__":
    main()
