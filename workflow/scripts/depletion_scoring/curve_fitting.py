"""
Sigmoid Curve Fitting for Depletion Analysis

This script fits sigmoid growth curves to depletion time-series data from
transposon insertion sequencing experiments. It processes multiple datasets
simultaneously and generates publication-quality plots with fitted parameters.

Typical Usage:
    python curve_fitting.py -i data.csv -t 0 2 4 6 8 10 12 14 -o results.csv

Input: CSV file with gene/insertion data as rows and time points as columns
Output: CSV file with fitted parameters and PDF with visualization plots
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple, Union
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
from sklearn.metrics import max_error
from tqdm import tqdm


# =============================== Constants ===============================
# Configure matplotlib for publication quality
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

LAM_PENALTY = 6e-3
TOL = 2e-6


# =============================== Configuration & Models ===============================
class CurveFittingConfig(BaseModel):
    """Pydantic model for validating curve fitting configuration."""
    input_file: Path = Field(..., description="Path to input CSV file with depletion data")
    output_file: Path = Field(..., description="Path to output CSV file for fitted parameters")
    time_points: List[float] = Field(..., description="Time points for the experiment")
    weight_file: Optional[Path] = Field(None, description="Path to weight CSV file")
    verbose: bool = Field(False, description="Enable verbose logging")

    @field_validator('input_file')
    def validate_input_file(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v

    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True)
        return v

    @field_validator('time_points')
    def validate_time_points(cls, v):
        if len(v) < 3:
            raise ValueError("At least 3 time points are required")
        return v

    class Config:
        frozen = True


class FittingResult(BaseModel):
    """Pydantic model for validating fitting results."""
    ID: str = Field(..., description="Gene/insertion identifier")
    Status: str = Field(..., description="Fitting status")
    A: float = Field(..., description="Maximum depletion level (asymptote)")
    um: float = Field(..., description="Maximum depletion rate")
    lam: float = Field(..., description="Lag time parameter")
    R2: float = Field(..., description="R-squared value")
    RMSE: float = Field(..., description="Root mean square error")
    normalized_RMSE: float = Field(..., description="Normalized RMSE")
    t10: float = Field(..., description="Time to reach 10% of the maximum depletion level")
    t50: float = Field(..., description="Time to reach 50% of the maximum depletion level")
    t90: float = Field(..., description="Time to reach 90% of the maximum depletion level")
    t_window: float = Field(..., description="Time window between t10 and t90")
    t_inflection: float = Field(..., description="Time of the inflection point")
    y_inflection: float = Field(..., description="Depletion level at the inflection point")
    auc: float = Field(..., description="Area under the curve")
    AIC: float = Field(..., description="Akaike Information Criterion")
    BIC: float = Field(..., description="Bayesian Information Criterion")


class SummaryStatistics(BaseModel):
    """Pydantic model for validating summary statistics."""
    total_datasets: int = Field(..., ge=0, description="Total number of datasets")
    successful_fits: int = Field(..., ge=0, description="Number of successful fits")
    success_rate: float = Field(..., ge=0.0, le=100.0, description="Success rate percentage")
    mean_R2: Optional[float] = Field(None, ge=0.0, le=1.0, description="Mean R-squared value")
    mean_RMSE: Optional[float] = Field(None, ge=0.0, description="Mean RMSE value")
    mean_A: Optional[float] = Field(None, description="Mean A parameter")
    mean_um: Optional[float] = Field(None, description="Mean um parameter")
    mean_t10: Optional[float] = Field(None, description="Mean t10 parameter")
    mean_t50: Optional[float] = Field(None, description="Mean t50 parameter")
    mean_t90: Optional[float] = Field(None, description="Mean t90 parameter")
    mean_t_window: Optional[float] = Field(None, description="Mean t_window parameter")
    mean_t_inflection: Optional[float] = Field(None, description="Mean t_inflection parameter")
    mean_y_inflection: Optional[float] = Field(None, description="Mean y_inflection parameter")
    mean_auc: Optional[float] = Field(None, description="Mean auc parameter")
    mean_AIC: Optional[float] = Field(None, description="Mean AIC value")
    mean_BIC: Optional[float] = Field(None, description="Mean BIC value")

# =============================== Setup Logging ===============================
def setup_logging(verbose: bool = False) -> None:
    """Configure loguru for the application."""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=level,
        colorize=False
    )


# =============================== Core Functions ===============================
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
def objective_function(params: List[float], x: np.ndarray, y: np.ndarray, 
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
def constraint_function1(params: List[float], t_last: float) -> float:
    """Constraint to ensure reasonable parameter bounds."""
    A, um, lam = params
    return t_last + 3 - abs(A) / abs(um) - lam


@logger.catch
def constraint_function2(params: List[float]) -> float:
    """Constraint to ensure smooth curve behavior."""
    A, um, lam = params
    x0 = lam + A / um / np.e
    val1 = float(np.abs(sigmoid_derivative(np.array([x0 - 1]), A, um, lam))[0])
    val2 = float(np.abs(sigmoid_derivative(np.array([x0 + 1]), A, um, lam))[0])
    return (val1 + val2 - 1.8 * abs(um))


@logger.catch
def fit_single_curve(x_values: np.ndarray, y_values: np.ndarray, 
                    weight_values: np.ndarray, ID: str, t_last: float) -> Dict[str, Union[str, float]]:
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
                      params: Dict[str, Union[str, float]], ID: str) -> None:
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
def process_depletion_data(input_file: Path, time_points: List[float], 
                          weight_file: Optional[Path] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], List[str]]:
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
    success_rate = (success_count / total_count * 100) if total_count > 0 else 0
    
    # Statistics for successful fits only
    successful_fits = results_df[results_df['Status'] == 'Success']
    
    stats = SummaryStatistics(
        total_datasets=total_count,
        successful_fits=success_count,
        success_rate=success_rate
    )
    
    if len(successful_fits) > 0:
        stats.mean_R2 = successful_fits['R2'].mean()
        stats.mean_RMSE = successful_fits['RMSE'].mean()
        stats.mean_A = successful_fits['A'].mean()
        stats.mean_um = successful_fits['um'].mean()
        stats.mean_t10 = successful_fits['t10'].mean()
        stats.mean_t50 = successful_fits['t50'].mean()
        stats.mean_t90 = successful_fits['t90'].mean()
        stats.mean_t_window = successful_fits['t_window'].mean()
        stats.mean_t_inflection = successful_fits['t_inflection'].mean()
        stats.mean_y_inflection = successful_fits['y_inflection'].mean()
        stats.mean_auc = successful_fits['auc'].mean()
        stats.mean_AIC = successful_fits['AIC'].mean()
        stats.mean_BIC = successful_fits['BIC'].mean()
    
    return stats


@logger.catch
def display_summary_table(stats: SummaryStatistics) -> None:
    """Display summary statistics in formatted table."""
    logger.info("=" * 50)
    logger.info("CURVE FITTING SUMMARY STATISTICS")
    logger.info("=" * 50)
    
    for key, value in stats.model_dump().items():
        if value is not None:
            if isinstance(value, float):
                logger.info(f"{key.replace('_', ' ').title():<25}: {value:.3f}")
            else:
                logger.info(f"{key.replace('_', ' ').title():<25}: {value}")
    
    logger.info("=" * 50)


# =============================== Main Function ===============================
def parse_arguments():
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


@logger.catch
def main():
    """Main entry point of the script."""
    start_time = time.time()
    
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args.verbose)
    
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
        sys.exit(1)
    
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


if __name__ == "__main__":
    main()