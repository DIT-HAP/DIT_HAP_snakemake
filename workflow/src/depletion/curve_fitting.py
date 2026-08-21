"""Sigmoid curve fitting domain logic for depletion analysis.

Fit Gompertz-type sigmoid growth curves to depletion time-series data from
transposon insertion sequencing experiments. Each dataset (gene or insertion)
is fitted independently by minimising a Huber loss with an L1 penalty on the
lag parameter, subject to smoothness and range constraints via SciPy's
``minimize``.

Input
-----
- Time/y-value arrays and weights, as prepared by the invoking script.

Output
------
- ``FittingResult`` / ``SummaryStatistics`` dataclasses and plain dicts of
  fitted parameters and derived metrics.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger
from scipy.optimize import minimize

from depletion.curve_model import sigmoid_function

# =============================================================================
# CONSTANTS
# =============================================================================
DL_PENALTY = 6e-3
TOL = 2e-6

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class FittingResult:
    """Fitting result schema for a single dataset."""
    ID: str
    Status: str
    A: float
    DR: float
    DL: float
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
    mean_DR: float | None = None
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
# CORE LOGIC
# =============================================================================
@logger.catch
def sigmoid_derivative(x: np.ndarray, A: float, DR: float, DL: float) -> np.ndarray:
    """Calculate derivative of sigmoid function using gompertz function."""
    alpha = (DR * np.e) / A
    u = alpha * (DL - x) + 1
    exponent = np.clip(u, -700, 700)
    return A * alpha * np.exp(exponent - np.exp(exponent))

@logger.catch
def time_at_p_effect(p: float, A: float, DR: float, DL: float) -> float:
    """Calculate the time at which the function reaches p proportion of its maximum effect."""
    return DL - (abs(A) / (abs(DR) * np.e)) * (np.log(-np.log(p)) - 1)


@logger.catch
def objective_function(params: list[float], x: np.ndarray, y: np.ndarray,
                      weight_values: np.ndarray) -> float:
    """Objective function for curve fitting using Huber loss."""
    A, DR, DL = params
    y_fit = sigmoid_function(x, A, DR, DL)
    residuals = y - y_fit
    z = (residuals * weight_values) ** 2

    # Huber loss for robustness to outliers
    rho_z = np.where(z <= 1, z, 2 * np.sqrt(z) - 1)

    # Add L1 regularization to DL
    dl_penalty = DL_PENALTY * abs(DL)
    return np.sum(rho_z) + dl_penalty


@logger.catch
def constraint_function1(params: list[float], t_last: float) -> float:
    """Constraint to ensure reasonable parameter bounds."""
    A, DR, DL = params
    return t_last + 3 - abs(A) / abs(DR) - DL


@logger.catch
def constraint_function2(params: list[float]) -> float:
    """Constraint to ensure smooth curve behavior."""
    A, DR, DL = params
    x0 = DL + A / DR / np.e
    val1 = float(np.abs(sigmoid_derivative(np.array([x0 - 1]), A, DR, DL))[0])
    val2 = float(np.abs(sigmoid_derivative(np.array([x0 + 1]), A, DR, DL))[0])
    return (val1 + val2 - 1.8 * abs(DR))


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
            x0=[-1, -1, 1],
            args=(x_values, y_values, weight_values),
            # bounds=((-1, t_last), (-1, np.inf), (-1e-6, t_last)),
            bounds=((-t_last, 1), (-np.inf, 1), (-1e-6, t_last)),
            constraints=constraints,
            options={'maxiter': 3000, 'disp': False},
            tol=TOL
        )

        if result.success:
            A, DR, DL = result.x
            residuals = y_values - sigmoid_function(x_values, A, DR, DL)
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            rmse = np.sqrt(ss_res / len(y_values))
            normalized_rmse = rmse / (y_values.max() - y_values.min())

            t_inflection = DL + abs(A) / (abs(DR) * np.e)
            y_inflection = A / np.e

            t10 = time_at_p_effect(0.1, A, DR, DL)
            t50 = time_at_p_effect(0.5, A, DR, DL)
            t90 = time_at_p_effect(0.9, A, DR, DL)
            t_window = t90 - t10

            # Calculate area under the curve (AUC) between curve and x-axis
            # Use numerical integration over the data range
            x_min, x_max = x_values.min(), x_values.max()
            x_integration = np.linspace(x_min, x_max, 1000)
            y_integration = sigmoid_function(x_integration, A, DR, DL)
            auc = np.trapezoid(y_integration, x_integration)

            # Calculate additional curve fitting metrics
            # Akaike Information Criterion (AIC)
            n_params = 3  # A, DR, DL
            n_points = len(y_values)
            aic = n_points * np.log(ss_res / n_points) + 2 * n_params
            # Bayesian Information Criterion (BIC)
            bic = n_points * np.log(ss_res / n_points) + n_params * np.log(n_points)

            return {
                'ID': ID,
                'Status': 'Success',
                'A': A, 'DR': DR, 'DL': DL, 't10': t10, 't50': t50, 't90': t90, 't_window': t_window, 't_inflection': t_inflection, 'y_inflection': y_inflection, 'auc': auc, 'AIC': aic, 'BIC': bic,
                'R2': r_squared, 'RMSE': rmse, 'normalized_RMSE': normalized_rmse,
            }
        else:
            logger.warning(f"Optimization failed for {ID}")
            return {
                'ID': ID,
                'Status': 'Optimization failed',
                'A': np.nan, 'DR': np.nan, 'DL': np.nan, 't10': np.nan, 't50': np.nan, 't90': np.nan, 't_window': np.nan, 't_inflection': np.nan, 'y_inflection': np.nan, 'auc': np.nan, 'AIC': np.nan, 'BIC': np.nan,
                'R2': np.nan, 'RMSE': np.nan, 'normalized_RMSE': np.nan,
            }

    except Exception as e:
        logger.error(f"Error fitting {ID}: {e}")
        return {
            'ID': ID,
            'Status': 'Fitting error',
            'A': np.nan, 'DR': np.nan, 'DL': np.nan, 't10': np.nan, 't50': np.nan, 't90': np.nan, 't_window': np.nan, 't_inflection': np.nan, 'y_inflection': np.nan, 'auc': np.nan, 'AIC': np.nan, 'BIC': np.nan,
            'R2': np.nan, 'RMSE': np.nan, 'normalized_RMSE': np.nan,
        }


def fit_and_augment(x_values: np.ndarray, y_data: np.ndarray,
                    weight_row: np.ndarray, ID: str, t_last: float,
                    timepoint_columns: list[str]) -> dict[str, str | float]:
    """Fit one dataset and attach its per-timepoint observed/fitted/residual columns.

    Module-level (picklable) worker so joblib's process pool can distribute it.
    Each dataset is independent, so results are identical to serial execution;
    joblib preserves input order, keeping the assembled DataFrame byte-for-byte
    equivalent to the previous single-core loop.
    """
    result = fit_single_curve(x_values, y_data, weight_row, ID, t_last)

    for j, _ in enumerate(x_values):
        result[timepoint_columns[j]] = round(y_data[j], 3)
    for j, time_val in enumerate(x_values):
        result[timepoint_columns[j] + '_fitted'] = round(
            sigmoid_function(time_val, result['A'], result['DR'], result['DL']), 3
        )
    for j, _ in enumerate(x_values):
        result[timepoint_columns[j] + '_residual'] = round(
            result[timepoint_columns[j]] - result[timepoint_columns[j] + '_fitted'], 3
        )

    return result


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
            mean_DR=successful_fits['DR'].mean(),
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
