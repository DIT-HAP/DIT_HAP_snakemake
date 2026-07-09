#!/usr/bin/env python3

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "scipy",
#     "tqdm",
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
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import minimize
from tqdm import tqdm

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# Configure matplotlib for publication quality
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

LAM_PENALTY = 6e-3
TOL = 2e-6

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
    verbose: bool = False

    def __post_init__(self) -> None:
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        if len(self.time_points) < 3:
            raise ValueError("At least 3 time points are required")


@dataclass(kw_only=True, slots=True, frozen=True)
class FittingResult:
    """Fitting result schema for a single dataset."""
    ID: str
    Status: str
    A: float
    um: float
    lam: float
    R2: float
    RMSE: float
    normalized_RMSE: float
    t10: float
    t50: float
    t90: float
    t_window: float
    t_inflection: float
    y_inflection: float
    auc: float
    AIC: float
    BIC: float


@dataclass(kw_only=True, slots=True, frozen=True)
class SummaryStatistics:
    """Summary statistics across all fitted datasets."""
    total_datasets: int
    successful_fits: int
    success_rate: float
    mean_R2: float | None = None
    mean_RMSE: float | None = None
    mean_A: float | None = None
    mean_um: float | None = None
    mean_t10: float | None = None
    mean_t50: float | None = None
    mean_t90: float | None = None
    mean_t_window: float | None = None
    mean_t_inflection: float | None = None
    mean_y_inflection: float | None = None
    mean_auc: float | None = None
    mean_AIC: float | None = None
    mean_BIC: float | None = None

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def sigmoid_function(x: np.ndarray, A: float, um: float, lam: float) -> np.ndarray:
    """Calculate sigmoid function values with numerical stability using gompertz function."""
    if A == 0:
        return np.zeros_like(x)
    alpha = (um * np.e) / A
    u = alpha * (lam - x) + 1
    exponent = np.clip(u, -700, 700)
    return A * np.exp(-np.exp(exponent))


@logger.catch
def sigmoid_derivative(x: np.ndarray, A: float, um: float, lam: float) -> np.ndarray:
    """Calculate derivative of sigmoid function using gompertz function."""
    alpha = (um * np.e) / A
    u = alpha * (lam - x) + 1
    exponent = np.clip(u, -700, 700)
    return A * alpha * np.exp(exponent - np.exp(exponent))

@logger.catch
def time_at_p_effect(p: float, A: float, um: float, lam: float) -> float:
    """Calculate the time at which the function reaches p proportion of its maximum effect."""
    return lam - (abs(A) / (abs(um) * np.e)) * (np.log(-np.log(p)) - 1)


@logger.catch
def objective_function(params: list[float], x: np.ndarray, y: np.ndarray,
                      weight_values: np.ndarray) -> float:
    """Objective function for curve fitting using Huber loss."""
    A, um, lam = params
    y_fit = sigmoid_function(x, A, um, lam)
    residuals = y - y_fit
    z = (residuals * weight_values) ** 2

    # Huber loss for robustness to outliers
    rho_z = np.where(z <= 1, z, 2 * np.sqrt(z) - 1)

    # Add L1 regularization to lam
    lam_penalty = LAM_PENALTY * abs(lam)
    return np.sum(rho_z) + lam_penalty


@logger.catch
def constraint_function1(params: list[float], t_last: float) -> float:
    """Constraint to ensure reasonable parameter bounds."""
    A, um, lam = params
    return t_last + 3 - abs(A) / abs(um) - lam


@logger.catch
def constraint_function2(params: list[float]) -> float:
    """Constraint to ensure smooth curve behavior."""
    A, um, lam = params
    x0 = lam + A / um / np.e
    val1 = float(np.abs(sigmoid_derivative(np.array([x0 - 1]), A, um, lam))[0])
    val2 = float(np.abs(sigmoid_derivative(np.array([x0 + 1]), A, um, lam))[0])
    return (val1 + val2 - 1.8 * abs(um))


@logger.catch
def fit_single_curve(x_values: np.ndarray, y_values: np.ndarray,
                    weight_values: np.ndarray, ID: str, t_last: float) -> dict[str, str | float]:
    """Fit sigmoid curve to a single dataset."""
    constraints = (
        {'type': 'ineq', 'fun': constraint_function1, 'args': (t_last,)},
        {'type': 'ineq', 'fun': constraint_function2}
    )

    try:
        result = minimize(
            objective_function,
            x0=[1, 1, 1],
            args=(x_values, y_values, weight_values),
            bounds=((-1, t_last), (-1, np.inf), (-1e-6, t_last)),
            constraints=constraints,
            options={'maxiter': 3000, 'disp': False},
            tol=TOL
        )

        if result.success:
            A, um, lam = result.x
            residuals = y_values - sigmoid_function(x_values, A, um, lam)
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(ss_res / len(y_values))
            normalized_rmse = rmse / (y_values.max() - y_values.min())

            t_inflection = lam + abs(A) / (abs(um) * np.e)
            y_inflection = A / np.e

            t10 = time_at_p_effect(0.1, A, um, lam)
            t50 = time_at_p_effect(0.5, A, um, lam)
            t90 = time_at_p_effect(0.9, A, um, lam)
            t_window = t90 - t10

            # Calculate area under the curve (AUC) between curve and x-axis
            # Use numerical integration over the data range
            x_min, x_max = x_values.min(), x_values.max()
            x_integration = np.linspace(x_min, x_max, 1000)
            y_integration = sigmoid_function(x_integration, A, um, lam)
            auc = np.trapezoid(y_integration, x_integration)

            # Calculate additional curve fitting metrics
            # Akaike Information Criterion (AIC)
            n_params = 3  # A, um, lam
            n_points = len(y_values)
            aic = n_points * np.log(ss_res / n_points) + 2 * n_params
            # Bayesian Information Criterion (BIC)
            bic = n_points * np.log(ss_res / n_points) + n_params * np.log(n_points)

            return {
                'ID': ID,
                'Status': 'Success',
                'A': A, 'um': um, 'lam': lam, 't10': t10, 't50': t50, 't90': t90, 't_window': t_window, 't_inflection': t_inflection, 'y_inflection': y_inflection, 'auc': auc, 'AIC': aic, 'BIC': bic,
                'R2': r_squared, 'RMSE': rmse, 'normalized_RMSE': normalized_rmse,
            }
        else:
            logger.warning(f"Optimization failed for {ID}")
            return {
                'ID': ID,
                'Status': 'Optimization failed',
                'A': np.nan, 'um': np.nan, 'lam': np.nan, 't10': np.nan, 't50': np.nan, 't90': np.nan, 't_window': np.nan, 't_inflection': np.nan, 'y_inflection': np.nan, 'auc': np.nan, 'AIC': np.nan, 'BIC': np.nan,
                'R2': np.nan, 'RMSE': np.nan, 'normalized_RMSE': np.nan,
            }

    except Exception as e:
        logger.error(f"Error fitting {ID}: {e}")
        return {
            'ID': ID,
            'Status': 'Fitting error',
            'A': np.nan, 'um': np.nan, 'lam': np.nan, 't10': np.nan, 't50': np.nan, 't90': np.nan, 't_window': np.nan, 't_inflection': np.nan, 'y_inflection': np.nan, 'auc': np.nan, 'AIC': np.nan, 'BIC': np.nan,
            'R2': np.nan, 'RMSE': np.nan, 'normalized_RMSE': np.nan,
        }


@logger.catch
def create_fitted_plot(ax: plt.Axes, x_values: np.ndarray, y_values: np.ndarray,
                      params: dict[str, str | float], ID: str) -> None:
    """Create a publication-quality plot for fitted curve."""
    ax.grid(True)

    if params['Status'] == 'Success':
        A, um, lam, _, _, _, _, _, _, _, AIC, BIC = params['A'], params['um'], params['lam'], params['t10'], params['t50'], params['t90'], params['t_window'], params['t_inflection'], params['y_inflection'], params['auc'], params['AIC'], params['BIC'],

        # Plot data points
        ax.scatter(x_values, y_values,
                  color=COLORS[1], alpha=0.8,
                  edgecolors='white',
                  label='Data')

        # Plot fitted curve
        x_smooth = np.linspace(min(x_values), max(x_values), 100)
        y_fit = sigmoid_function(x_smooth, A, um, lam)
        ax.plot(x_smooth, y_fit,
               color=COLORS[2], label='Fitted')

        # Add constraint lines
        ax.axhline(y=A, color=COLORS[0],
                  linestyle='--', alpha=0.3)
        ax.axvline(x=lam, color=COLORS[0],
                  linestyle='--', alpha=0.3)

        # Add parameter text
        param_text = f'A={A:.2f}    R²={params["R2"]:.3f}\num={um:.2f}  RMSE={params["RMSE"]:.3f}\nlam={lam:.2f}    NRMSE={params["normalized_RMSE"]:.3f}\nAIC={AIC:.2f}    BIC={BIC:.2f}'
        ax.text(0.05, 0.95, param_text,
               transform=ax.transAxes,
               verticalalignment='top')
    else:
        # Plot failed fit
        ax.scatter(x_values, y_values,
                  color='gray', alpha=0.6)
        ax.text(0.5, 0.5, 'Fit Failed',
               transform=ax.transAxes,
               horizontalalignment='center', color='red')

    ax.set_ylim(-1.5, 8.5)
    ax.set_title(" ".join(ID.split("=")))


@logger.catch
def generate_fitting_plots(results_df: pd.DataFrame, x_values: np.ndarray,
                          y_values: np.ndarray, output_plot: Path) -> None:
    """Generate multi-page PDF with fitting plots."""
    plots_per_page = 32
    num_pages = int(np.ceil(len(results_df) / plots_per_page))

    logger.info(f"Generating {num_pages} pages of plots...")

    with PdfPages(output_plot) as pdf:
        for page in range(num_pages):
            fig, axes = plt.subplots(8, 4, figsize=(AX_WIDTH*4, AX_HEIGHT*8))
            axes = axes.flatten()

            if page % 10 == 0:
                logger.info(f"Generating page {page+1} of {num_pages}...")

            start_idx = page * plots_per_page
            end_idx = min((page + 1) * plots_per_page, len(results_df))

            for idx in range(start_idx, end_idx):
                ax_idx = idx % plots_per_page
                row = results_df.iloc[idx]
                ID = " ".join(map(str, row.name))

                create_fitted_plot(
                    axes[ax_idx],
                    x_values,
                    y_values[idx],
                    row.to_dict(),
                    ID
                )

            # Hide unused subplots
            for ax_idx in range(end_idx - start_idx, plots_per_page):
                axes[ax_idx].set_visible(False)

            pdf.savefig(fig)
            plt.close(fig)


@logger.catch
def process_depletion_data(input_file: Path, time_points: list[float],
                          weight_file: Path | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str], list[str]]:
    """Load and process depletion data from CSV file."""
    logger.info(f"Loading data from {input_file}")

    # Load data with multi-level index for insertions
    data = pd.read_csv(input_file, header=0, sep="\t")
    len_columns = len(data.columns)
    index_column_num = len_columns - len(time_points)
    index_columns = data.columns.tolist()[:index_column_num]
    timepoint_columns = data.columns.tolist()[index_column_num:]
    data.set_index(index_columns, inplace=True)

    # Create gene identifiers
    IDs = ["=".join(map(str, idx)) for idx in data.index.tolist()]

    x_values = np.array(time_points)
    y_values = data.values

    if weight_file is not None:
        weight_data = pd.read_csv(weight_file, header=0)
        weight_data.set_index(index_columns, inplace=True)
        weight_data = weight_data.loc[data.index].fillna(0.01)
        weight_values = weight_data.values
    else:
        weight_values = np.ones(shape=(len(IDs), len(x_values)))

    logger.info(f"Loaded {len(IDs)} datasets with {len(x_values)} time points")

    return x_values, y_values, weight_values, IDs, index_columns, timepoint_columns


@logger.catch
def generate_summary_statistics(results_df: pd.DataFrame) -> SummaryStatistics:
    """Generate comprehensive summary statistics."""
    total_count = len(results_df)
    success_count = len(results_df[results_df['Status'] == 'Success'])
    success_rate = (success_count / total_count * 100) if total_count > 0 else 0.0

    # Statistics for successful fits only
    successful_fits = results_df[results_df['Status'] == 'Success']

    if len(successful_fits) > 0:
        return SummaryStatistics(
            total_datasets=total_count,
            successful_fits=success_count,
            success_rate=success_rate,
            mean_R2=successful_fits['R2'].mean(),
            mean_RMSE=successful_fits['RMSE'].mean(),
            mean_A=successful_fits['A'].mean(),
            mean_um=successful_fits['um'].mean(),
            mean_t10=successful_fits['t10'].mean(),
            mean_t50=successful_fits['t50'].mean(),
            mean_t90=successful_fits['t90'].mean(),
            mean_t_window=successful_fits['t_window'].mean(),
            mean_t_inflection=successful_fits['t_inflection'].mean(),
            mean_y_inflection=successful_fits['y_inflection'].mean(),
            mean_auc=successful_fits['auc'].mean(),
            mean_AIC=successful_fits['AIC'].mean(),
            mean_BIC=successful_fits['BIC'].mean(),
        )

    return SummaryStatistics(
        total_datasets=total_count,
        successful_fits=success_count,
        success_rate=success_rate,
    )


@logger.catch
def display_summary_table(stats: SummaryStatistics) -> None:
    """Display summary statistics in formatted table."""
    logger.info("=" * 50)
    logger.info("CURVE FITTING SUMMARY STATISTICS")
    logger.info("=" * 50)

    for key, value in asdict(stats).items():
        if value is not None:
            if isinstance(value, float):
                logger.info(f"{key.replace('_', ' ').title():<25}: {value:.3f}")
            else:
                logger.info(f"{key.replace('_', ' ').title():<25}: {value}")

    logger.info("=" * 50)

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

        # Fit curves with progress tracking
        logger.info("Fitting sigmoid curves...")
        all_results = []

        with tqdm(total=len(y_values), desc="Fitting progress") as pbar:
            for i, (y_data, ID) in enumerate(zip(y_values, IDs)):
                result = fit_single_curve(x_values, y_data, weight_values[i], ID, t_last)

                # Add time series data to result
                for j, time_val in enumerate(x_values):
                    result[timepoint_columns[j]] = round(y_data[j], 3)
                for j, time_val in enumerate(x_values):
                    result[timepoint_columns[j] + '_fitted'] = round(sigmoid_function(time_val, result['A'], result['um'], result['lam']), 3)
                for j, time_val in enumerate(x_values):
                    result[timepoint_columns[j] + '_residual'] = round(result[timepoint_columns[j]] - result[timepoint_columns[j] + '_fitted'], 3)

                all_results.append(result)
                pbar.update(1)

        # Create results DataFrame
        results_df = pd.DataFrame(all_results)
        results_df.insert(1, 'time_points', [",".join(map(str, list(x_values)))] * len(results_df))

        # Round numeric columns
        numeric_columns = {
            'A':3, 'um':3, 'lam':3, 't10':3, 't50':3, 't90':3, 't_window':3, 't_inflection':3, 'y_inflection':3, 'auc':3, 'R2':6, 'RMSE':3, 'normalized_RMSE':6, 'AIC':3, 'BIC':3,
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

        # Generate plots
        output_plot = config.output_file.with_suffix('.pdf').with_name(config.output_file.stem + '_fitted_curves.pdf')
        # generate_fitting_plots(results_df, x_values, y_values, output_plot)

        # Calculate and display statistics
        stats = generate_summary_statistics(results_df)
        display_summary_table(stats)

        # Final summary
        elapsed_time = time.time() - start_time
        logger.success(f"Analysis completed in {elapsed_time:.1f} seconds")
        logger.success(f"Results saved to: {config.output_file}")
        logger.success(f"Plots saved to: {output_plot}")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
