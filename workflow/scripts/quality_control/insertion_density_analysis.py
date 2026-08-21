#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
#     "numpy",
#     "pandas",
# ]
# ///
"""
Insertion Density Analysis for Transposon Insertion Sequencing
==============================================================

Analyzes insertion density patterns in transposon sequencing data, comparing
the initial and final timepoints to reveal the classic piggyBac phenotyping
signal: insertion sites well-covered at the initial timepoint that drop out
by the final timepoint mark depleted (typically essential) genes.

Loads raw insertion read counts (summed across biological replicates) and
genomic annotations, filters for in-gene insertions using established
criteria (non-intergenic and distance to stop codon > 4), and computes
per-gene density metrics at both timepoints plus their change.

For each gene it derives insertion-site density per kilobase at the initial
and final timepoint (and their difference / log2 fold-change), gap-distribution
metrics (including a Gini coefficient of insertion location), read-depth
inequality (Gini coefficient of depth) at both timepoints, and strand-preference
measures.

This is the computation layer: it writes only the density statistics table.
Rendering (histograms, initial-vs-final scatter panels) lives in
``workflow/scripts/figures/plot_insertion_density.py``.

Input
-----
- Insertion data (``-i``): tab-separated file with a 4-level row index and a
  2-level column header (Sample, Timepoint). Read counts are summed across
  samples (replicates) at the initial (``-t``) and final (``-f``) timepoints.
- Annotation data (``-a``): tab-separated file with a 4-level index and gene
  annotation columns (Type, Distance_to_stop_codon, Systematic ID, Name,
  FYPOviability, Chr_Interval, Strand_Interval, ParentalRegion_start,
  ParentalRegion_end, ParentalRegion_length, Insertion_direction).

Output
------
- Density statistics table (``-o``): tab-separated, one row per gene.
  Columns: Systematic ID (index), Name, Chr, Start, End, Length, Strand,
  FYPOviability, total_insertions, unique_sites_initial, unique_sites_final,
  gene_length, insertion_density_per_kb_initial, insertion_density_per_kb_final,
  insertion_density_change, insertion_density_log2fc, num_gaps, largest_gap,
  largest_gap_fraction, smallest_gap, smallest_gap_fraction, mean_gap_length,
  mean_gap_length_fraction, median_gap_length, median_gap_length_fraction,
  gap_length_sd, gap_length_sd_fraction, all_gap_lengths, all_gap_lengths_fraction,
  gini_coefficient_of_location, total_reads_initial, total_reads_final,
  mean_reads_per_insertion_initial, mean_reads_per_insertion_final,
  gini_coefficient_of_depth_initial, gini_coefficient_of_depth_final,
  total_reads_log2fc, forward_insertions, reverse_insertions, forward_preference,
  reverse_preference, strand_bias, paired_sites, paired_sites_fraction.

Usage
-----
    python insertion_density_analysis.py -i raw_reads.filtered.tsv -a annotations.tsv -o density_stats.tsv -t YES0 -f YES4

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  3.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402
from qc.density import (  # noqa: E402
    AnalysisResult,
    load_insertion_data,
    load_annotation_data,
    filter_in_gene_insertions,
    calculate_insertion_statistics,
    calculate_gap_statistics,
    calculate_gini_coefficient,
    calculate_read_statistics,
    calculate_strand_statistics,
    analyze_gene_insertions,
    generate_summary_statistics,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths and the initial/final timepoint column names."""

    insertion_data_path: Path
    annotations_path: Path
    output_path: Path
    initial_timepoint: str
    final_timepoint: str

    def __post_init__(self) -> None:
        for path in (self.insertion_data_path, self.annotations_path):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
            if not path.is_file():
                raise ValueError(f"Input path is not a file: {path}")
            if path.suffix.lower() not in ['.csv', '.tsv', '.txt']:
                raise ValueError(f"Input file must be CSV or TSV format. Got: {path.suffix}")
            if path.stat().st_size == 0:
                raise ValueError(f"Input file is empty: {path}")

        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        if self.output_path.suffix.lower() not in ['.csv', '.tsv', '.txt']:
            raise ValueError(f"Output file must be CSV or TSV format. Got: {self.output_path.suffix}")

        if not self.initial_timepoint:
            raise ValueError("initial_timepoint must be a non-empty string")
        if not self.final_timepoint:
            raise ValueError("final_timepoint must be a non-empty string")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Insertion density analysis script")
    parser.add_argument("-i", "--insertion_data_path", type=Path, required=True, help="Input TSV file with raw insertion read counts (4-level index, 2-level (Sample, Timepoint) header)")
    parser.add_argument("-a", "--annotations_path", type=Path, required=True, help="Input TSV file with annotations")
    parser.add_argument("-o", "--output_path", type=Path, required=True, help="Output TSV file with density statistics")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint column name")
    parser.add_argument("-f", "--final_timepoint", type=str, required=True, help="Final timepoint column name")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point of the script."""
    args = parse_args()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(log_level)

    # Validate input and output paths using the config dataclass
    start_time = time.time()

    try:
        config = InputOutputConfig(
            insertion_data_path=args.insertion_data_path,
            annotations_path=args.annotations_path,
            output_path=args.output_path,
            initial_timepoint=args.initial_timepoint,
            final_timepoint=args.final_timepoint,
        )

        # Get file sizes for result tracking
        insertion_file_size = config.insertion_data_path.stat().st_size
        annotation_file_size = config.annotations_path.stat().st_size

        logger.info("Starting insertion density analysis")
        logger.info(f"Insertion data file: {config.insertion_data_path} ({insertion_file_size:,} bytes)")
        logger.info(f"Annotations file: {config.annotations_path} ({annotation_file_size:,} bytes)")
        logger.info(f"Initial timepoint: {config.initial_timepoint}")
        logger.info(f"Final timepoint: {config.final_timepoint}")

        # Load data
        insertion_data = load_insertion_data(
            config.insertion_data_path, config.initial_timepoint, config.final_timepoint
        )
        annotations = load_annotation_data(config.annotations_path)

        # Filter for in-gene insertions
        in_gene_insertions = filter_in_gene_insertions(insertion_data, annotations)
        logger.info(f"Analyzing {in_gene_insertions['Systematic ID'].nunique()} genes")

        # Analyze each gene via a single grouped pass (avoids re-filtering the full
        # frame once per gene, which was quadratic in the number of genes)
        results_df = in_gene_insertions.groupby('Systematic ID', sort=False).apply(
            lambda group: pd.Series(analyze_gene_insertions(group.name, group)),
            include_groups=False,
        )
        results_df = results_df.drop(columns='Systematic ID')
        results_df.index.name = 'Systematic ID'

        # Generate summary statistics
        stats = generate_summary_statistics(results_df)

        # Update result object with final statistics
        end_time = time.time()
        analysis_duration = end_time - start_time

        res = AnalysisResult(
            total_genes_analyzed=stats['total_genes_analyzed'],
            total_insertions_analyzed=stats['total_insertions_analyzed'],
            mean_insertion_density_per_kb_initial=stats['mean_insertion_density_per_kb_initial'],
            mean_insertion_density_per_kb_final=stats['mean_insertion_density_per_kb_final'],
            mean_insertion_density_log2fc=stats['mean_insertion_density_log2fc'],
            mean_gini_coefficient_of_depth_initial=stats['mean_gini_coefficient_of_depth_initial'],
            mean_strand_bias=stats['mean_strand_bias'],
        )

        # Save results
        results_df.to_csv(config.output_path, index=True, sep="\t")

        # Final summary
        logger.success("Analysis completed successfully")
        logger.success(f"Results saved to {config.output_path}")
        logger.success(f"Analyzed {len(results_df)} genes with insertion data")

        # Log summary statistics
        summary = asdict(res)
        logger.info(f"Analysis summary: {summary}")
        logger.info(f"Performance: {analysis_duration:.2f} seconds for {res.total_genes_analyzed} genes")

    except ValueError as e:
        logger.error(f"Validation error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

