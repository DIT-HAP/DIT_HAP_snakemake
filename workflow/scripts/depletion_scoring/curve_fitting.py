#!/usr/bin/env python3

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "joblib",
#     "loguru",
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "scipy",
# ]
# ///

"""
Sigmoid Curve Fitting for Depletion Analysis
============================================

Fit Gompertz-type sigmoid growth curves to depletion time-series data from
transposon insertion sequencing experiments. Each dataset (gene or insertion)
is fitted independently by minimising a Huber loss with an L1 penalty on the
lag parameter, subject to smoothness and range constraints via SciPy's
``minimize``. Fitted parameters and derived metrics (R2, RMSE, AIC, BIC,
inflection times, AUC) are written to a tab-separated table alongside the
per-timepoint fitted values and residuals.

Input
-----
- TSV file with one or more gene/insertion identifier columns followed by one
  column per time point (log fold-change values). Passed via ``-i/--input``.
- Optional TSV weight file with matching index columns (``-w/--weight``).

Output
------
- Main TSV of fitted parameters and metrics (``-o/--output``).
- ``fitting_LFCs.tsv`` and ``fitting_results.tsv`` written next to the output.
- A ``*_fitted_curves.pdf`` path is derived for optional plotting.

Usage
-----
    python curve_fitting.py -i data.tsv -t 0 2 4 6 8 10 12 14 -o results.tsv
    python curve_fitting.py -i data.tsv -t 0 2 4 6 8 -o results.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

# Force BLAS/OpenMP threads to 1 BEFORE numpy is imported. Curve fitting is
# parallelised across datasets with joblib (process pool), so per-worker BLAS
# threads add nothing but oversubscription. More importantly, multi-threaded
# BLAS makes scipy.minimize non-deterministic on ill-conditioned fits: the same
# input yields different local optima depending on thread count, so results
# would vary with how the script is launched. Snakemake exports
# OMP_NUM_THREADS/OPENBLAS_NUM_THREADS = the rule's `threads` (see snakemake
# shell.py), so setdefault() would NOT override it and the pin would silently
# fail under snakemake. Assign unconditionally to guarantee reproducible,
# thread-count-independent output. loky re-imports this module in each worker.
for _thr_var in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "GOTO_NUM_THREADS",
):
    os.environ[_thr_var] = "1"

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from joblib import Parallel, delayed
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from depletion.curve_fitting import (  # noqa: E402
    fit_and_augment,
    generate_summary_statistics,
    display_summary_table,
    process_depletion_data,
)
from depletion.curve_model import sigmoid_function  # noqa: E402

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class CurveFittingConfig:
    """Validated curve fitting configuration."""
    input_file: Path
    output_file: Path
    time_points: list[float]
    weight_file: Path | None = None
    jobs: int = 1
    verbose: bool = False

    def __post_init__(self) -> None:
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        if len(self.time_points) < 3:
            raise ValueError("At least 3 time points are required")
        if self.jobs == 0 or self.jobs < -1:
            raise ValueError(f"jobs must be -1 or a positive integer, got {self.jobs}")

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Fit sigmoid curves to depletion time-series data",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--input", type=Path, required=True,
                       help="Path to input TSV file with depletion data")
    parser.add_argument("-w", "--weight", type=Path, required=False, default=None,
                       help="Path to weight TSV file")
    parser.add_argument("-t", "--time_points", required=True, nargs='+',
                       type=float, help="Time points for the experiment")
    parser.add_argument("-o", "--output", type=Path, required=True,
                       help="Path to output TSV file for fitted parameters")
    parser.add_argument("-j", "--jobs", type=int, default=1,
                       help="Parallel worker processes for curve fitting "
                            "(1 = serial, -1 = all cores; default: 1)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose logging")

    return parser.parse_args()


def main() -> int:
    """Main entry point of the script."""
    start_time = time.time()

    # Parse arguments and setup
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Validate configuration
    try:
        config = CurveFittingConfig(
            input_file=args.input,
            output_file=args.output,
            time_points=args.time_points,
            weight_file=args.weight,
            jobs=args.jobs,
            verbose=args.verbose
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    try:
        logger.info("Starting sigmoid curve fitting analysis")
        logger.info(f"Input file: {config.input_file}")
        logger.info(f"Time points: {config.time_points}")

        # Process data
        x_values, y_values, weight_values, IDs, index_columns, timepoint_columns = process_depletion_data(
            config.input_file, config.time_points, config.weight_file
        )
        t_last = x_values[-1]

        # Fit curves in parallel across datasets. Each dataset is independent;
        # joblib preserves input order so `all_results` matches the previous
        # serial loop element-for-element.
        logger.info(f"Fitting sigmoid curves with {config.jobs} worker(s)...")
        all_results = Parallel(n_jobs=config.jobs)(
            delayed(fit_and_augment)(
                x_values, y_data, weight_values[i], ID, t_last, timepoint_columns
            )
            for i, (y_data, ID) in enumerate(zip(y_values, IDs))
        )

        # Create results DataFrame
        results_df = pd.DataFrame(all_results)
        results_df.insert(1, 'time_points', [",".join(map(str, list(x_values)))] * len(results_df))

        # Round numeric columns
        numeric_columns = {
            'A':3, 'DR':3, 'DL':3, 't10':3, 't50':3, 't90':3, 't_window':3, 't_inflection':3, 'y_inflection':3, 'auc':3, 'R2':6, 'RMSE':3, 'normalized_RMSE':6, 'AIC':3, 'BIC':3,
        }
        results_df[list(numeric_columns.keys())] = results_df[list(numeric_columns.keys())].round(numeric_columns)

        # Set multi-level index
        results_df.set_index("ID", inplace=True)
        multiple_index = pd.MultiIndex.from_tuples([idx.split("=") for idx in results_df.index.tolist()])
        results_df.index = multiple_index
        results_df.rename_axis(index_columns, inplace=True)

        # Save results
        results_df.to_csv(config.output_file, index=True, sep="\t")

        # fitted_LFCs
        fitting_LFCs = results_df.filter(like="fitted")
        fitting_LFCs.columns = fitting_LFCs.columns.str.replace("_fitted", "")
        fitting_LFCs.to_csv(config.output_file.parent/"fitting_LFCs.tsv", index=True, sep="\t")
        # fitted_results
        results_df[list(numeric_columns.keys())].to_csv(config.output_file.parent/"fitting_results.tsv", index=True, sep="\t")

        # Plotting moved to separate rendering script (plot_curve_fitting.py)

        # Calculate and display statistics
        stats = generate_summary_statistics(results_df)
        display_summary_table(stats)

        # Final summary
        elapsed_time = time.time() - start_time
        logger.success(f"Analysis completed in {elapsed_time:.1f} seconds")
        logger.success(f"Results saved to: {config.output_file}")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
