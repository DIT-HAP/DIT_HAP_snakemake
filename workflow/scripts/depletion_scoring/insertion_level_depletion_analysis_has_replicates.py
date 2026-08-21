#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
#     "numpy",
#     "pandas",
#     "pydeseq2",
# ]
# ///

"""
Insertion-Level Depletion Analysis (Replicated Samples)
=======================================================

Differential-abundance analysis of transposon insertion counts using PyDESeq2
(DESeq2). The script identifies insertions that deplete across time points
relative to an initial time point, using replicated samples.

Counts are read with a four-level row MultiIndex (chromosome, coordinate,
strand, target) and a two-level column MultiIndex (group, condition). Size
factors are estimated from a supplied set of control insertions, dispersions
and log2 fold changes are fitted, Cook's-distance outliers are refit, and a
Wald test is run for each non-initial time point versus the initial time point.
Log2 fold changes and Wald statistics are negated so that depletion is reported
as a negative fold change. Per-statistic tables (baseMean, log2FoldChange,
lfcSE, stat, pvalue, padj) and normalized/raw count matrices are written out.

Input
-----
- Counts TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex, ``header=[0, 1]`` column
  MultiIndex, tab-separated integer insertion counts.
- Control-insertions TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex, tab-separated;
  its index selects the control insertions used for size-factor estimation.

Output
------
- ``log2FoldChange.tsv`` (the ``-o`` path): primary log2 fold-change table.
- Alongside it in the same directory: ``insertion_level_statistics.tsv``,
  ``baseMean.tsv``, ``lfcSE.tsv``, ``stat.tsv``, ``pvalue.tsv``, ``padj.tsv``,
  ``normed_counts.tsv``, ``count_X.tsv``, ``cooks.tsv``.
- Dispersion figure-data TSV for the rendering layer, at the explicit
  ``--dispersion_data`` path: per-insertion normed mean with genewise/MAP/fitted
  dispersions.

Usage
-----
    python insertion_level_depletion_analysis_has_replicates.py \
        -i counts.tsv -c control_insertions.tsv -t 0h -o output_dir/log2FoldChange.tsv \
        --dispersion_data figure_dir/dispersion_data.tsv
    python insertion_level_depletion_analysis_has_replicates.py \
        -i counts.tsv -c control_insertions.tsv -t 0h -o output_dir/log2FoldChange.tsv \
        --dispersion_data figure_dir/dispersion_data.tsv --verbose

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
import time
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402
from depletion.insertion_level_replicates import (  # noqa: E402
    AnalysisResult,
    load_and_preprocess_data,
    create_deseq_dataset,
    perform_differential_analysis,
    write_dispersion_data_tsv,
    concatenate_results,
    transform_index_to_multiindex,
)
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths and run flags for the analysis."""
    counts_file: Path
    control_insertions_file: Path
    initial_timepoint: str
    output_file: Path
    dispersion_data_file: Path
    verbose: bool = False

    def __post_init__(self) -> None:
        for path in (self.counts_file, self.control_insertions_file):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
        for path in (self.output_file, self.dispersion_data_file):
            path.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perform differential expression analysis on insertion counts.")
    parser.add_argument("-i", "--counts_file", type=Path, required=True, help="Path to the counts file.")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint to analyze.")
    parser.add_argument("-c", "--control_insertions_file", type=Path, required=True, help="Path to the control insertions file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file.")
    parser.add_argument("--dispersion_data", type=Path, required=True, help="Path to the dispersion figure-data TSV.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point for insertion-level depletion analysis."""
    start_time = time.time()

    args = parse_args()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(log_level)

    try:
        # Validate input and output paths using the config dataclass
        config = InputOutputConfig(
            counts_file=args.counts_file,
            control_insertions_file=args.control_insertions_file,
            initial_timepoint=args.initial_timepoint,
            output_file=args.output,
            dispersion_data_file=args.dispersion_data,
            verbose=args.verbose
        )

        logger.info("Starting insertion-level depletion analysis")
        logger.info(f"Counts file: {config.counts_file}")
        logger.info(f"Control insertions file: {config.control_insertions_file}")
        logger.info(f"Initial timepoint: {config.initial_timepoint}")
        logger.info(f"Output file: {config.output_file}")

        # Load and preprocess data
        counts_df, metadata, counts_df_index_names, counts_df_columns_names, control_insertions = load_and_preprocess_data(
            config.counts_file, config.control_insertions_file
        )

        logger.info("Metadata for analysis:")
        logger.info(f"\n{metadata}")
        logger.info(f"Number of control insertions: {len(control_insertions)}")

        timepoints = metadata["condition"].unique().tolist()
        timepoints.remove(config.initial_timepoint)
        logger.info(f"Control timepoint: {config.initial_timepoint}")
        logger.info(f"Timepoints for analysis: {timepoints}")

        # Create DESeq2 dataset
        dds = create_deseq_dataset(counts_df, metadata, control_insertions, config.initial_timepoint)
        write_dispersion_data_tsv(dds, config.dispersion_data_file)
        logger.info(f"Dispersion data saved to {config.dispersion_data_file}")

        # Transform normalized counts
        logger.info("Transforming normalized counts to multi-index...")
        normalized_counts = transform_index_to_multiindex(dds, "normed_counts")
        normalized_counts = normalized_counts.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
        normalized_counts.to_csv(config.output_file.parent / "normed_counts.tsv", index=True, float_format="%.3f", sep="\t")

        # Transform count matrix
        logger.info("Transforming count matrix to multi-index...")
        count_X = pd.DataFrame(dds.X, index=dds.obs.index.tolist(), columns=dds.var.index.tolist()).T
        count_X.index = pd.MultiIndex.from_tuples(count_X.index.str.split("=").tolist())
        count_X.columns = pd.MultiIndex.from_tuples(count_X.columns.str.split("#").tolist())
        count_X = count_X.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
        count_X.to_csv(config.output_file.parent / "count_X.tsv", index=True, float_format="%.3f", sep="\t")

        # Transform Cook's distances
        cooks_df = transform_index_to_multiindex(dds, "cooks")
        cooks_df = cooks_df.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
        cooks_df.to_csv(config.output_file.parent / "cooks.tsv", index=True, float_format="%.3f", sep="\t")

        # Perform differential analysis
        logger.info("Performing differential analysis...")
        stat_res = perform_differential_analysis(dds, timepoints, config.initial_timepoint)

        # Concatenate results
        logger.info("Concatenating results...")
        concated_results = concatenate_results(stat_res, timepoints)
        concated_results = concated_results.rename_axis(counts_df_index_names, axis=0).rename_axis(["Timepoint", "Statistic"], axis=1)

        # Add the metrics for the initial timepoint
        logger.info("Adding metrics for the initial timepoint...")
        logger.info("All timepoints share the same baseMean")
        baseMean_initial = concated_results.xs("baseMean", axis=1, level="Statistic").iloc[:,0]
        logger.info("Set the initial log2FoldChange to 0")
        log2FoldChange_initial = 0
        logger.info("Set the initial lfcSE to NaN")
        lfcSE_initial = np.nan
        logger.info("Set the initial stat to NaN")
        stat_initial = np.nan
        logger.info("Set the initial pvalue to 1")
        pvalue_initial = 1
        logger.info("Set the initial padj to 1")
        padj_initial = 1
        logger.info("Insert the initial timepoint metrics to the first column...")
        concated_results.insert(0, (config.initial_timepoint,"padj"), padj_initial)
        concated_results.insert(0, (config.initial_timepoint,"pvalue"), pvalue_initial)
        concated_results.insert(0, (config.initial_timepoint,"stat"), stat_initial)
        concated_results.insert(0, (config.initial_timepoint,"lfcSE"), lfcSE_initial)
        concated_results.insert(0, (config.initial_timepoint,"log2FoldChange"), log2FoldChange_initial)
        concated_results.insert(0, (config.initial_timepoint,"baseMean"), baseMean_initial)

        numeric_columns = {"baseMean": 3, "log2FoldChange": 3, "lfcSE": 3, "stat": 3, "pvalue": 6, "padj": 6}
        logger.info("Rounding numeric columns...")
        for stat_name, decimal_places in numeric_columns.items():
            stat_columns = concated_results.xs(stat_name, axis=1, level="Statistic", drop_level=False)
            concated_results[stat_columns.columns] = stat_columns.round(decimal_places)
        logger.info("Saving results...")
        concated_results.to_csv(config.output_file.parent / "insertion_level_statistics.tsv", sep="\t")

        # Save individual statistic files
        baseMean_df = concated_results.xs("baseMean", axis=1, level="Statistic")
        baseMean_df.to_csv(config.output_file.parent/"baseMean.tsv", index=True, sep="\t")

        LFC_df = concated_results.xs("log2FoldChange", axis=1, level="Statistic")
        LFC_df.to_csv(config.output_file, index=True, sep="\t")

        lfcSE_df = concated_results.xs("lfcSE", axis=1, level="Statistic")
        lfcSE_df.to_csv(config.output_file.parent/"lfcSE.tsv", index=True, sep="\t")

        stat_df = concated_results.xs("stat", axis=1, level="Statistic")
        stat_df.to_csv(config.output_file.parent/"stat.tsv", index=True, sep="\t")

        pvalue_df = concated_results.xs("pvalue", axis=1, level="Statistic")
        pvalue_df.to_csv(config.output_file.parent/"pvalue.tsv", index=True, sep="\t")

        padj_df = concated_results.xs("padj", axis=1, level="Statistic")
        padj_df.to_csv(config.output_file.parent/"padj.tsv", index=True, sep="\t")

        # Create analysis result
        end_time = time.time()
        execution_time = end_time - start_time
        result = AnalysisResult(
            total_insertions_analyzed=len(counts_df.columns),
            timepoints_processed=len(timepoints),
            control_insertions_count=len(control_insertions),
            execution_time=execution_time
        )

        logger.success(f"Analysis completed in {execution_time:.2f} seconds")
        logger.info(f"Analyzed {result.total_insertions_analyzed} insertions across {result.timepoints_processed} timepoints")
        logger.info(f"Used {result.control_insertions_count} control insertions for normalization")
        logger.success(f"Results saved to {config.output_file.parent}")

    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
