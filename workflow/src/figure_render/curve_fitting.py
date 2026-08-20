"""Curve fitting figure rendering.

Input
-----
- Fitting statistics TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex,
  tab-separated, with columns ``Status``, ``A``, ``DR``, ``DL``, ``R2``,
  ``RMSE``, ``time_points``.
- LFC TSV: same row MultiIndex, tab-separated, with one column per timepoint.

Output
------
- ``<stem>.pdf`` / ``<stem>.review.png`` — grid of observed vs fitted sigmoid
  curves for a random sample of genes.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns
import numpy as np
import pandas as pd
from loguru import logger

from depletion.curve_model import sigmoid_function
from figures import apply_house_style, save_dual

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_sample_data(
    fitting_stats_path: Path, lfc_path: Path, n_curves: int, random_seed: int
) -> tuple[pd.DataFrame, pd.DataFrame, list[float]]:
    """Load fitting stats and LFC data, sample successful fits, return aligned data and time points."""
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    fitting_stats = pd.read_csv(fitting_stats_path, sep='\t', index_col=[0, 1, 2, 3])

    logger.info(f"Loaded {len(fitting_stats)} rows")

    # Filter to successful fits
    successful = fitting_stats[fitting_stats['Status'] == 'Success'].copy()
    logger.info(f"Found {len(successful)} successful fits")

    if successful.empty:
        logger.warning("No successful fits found!")
        return successful, pd.DataFrame(), []

    # Sample n_curves
    n_sample = min(n_curves, len(successful))
    sampled = successful.sample(n=n_sample, random_state=random_seed)
    logger.info(f"Sampled {len(sampled)} curves for plotting")

    # Load LFC data
    logger.info(f"Loading LFC data from {lfc_path}...")
    lfc_data = pd.read_csv(lfc_path, sep='\t', index_col=[0, 1, 2, 3])

    # Align with sampled indices
    lfc_sampled = lfc_data.loc[sampled.index]

    # Extract time points from time_points column
    time_points_str = sampled['time_points'].iloc[0]
    time_points = [float(t) for t in time_points_str.split(',')]
    logger.info(f"Time points: {time_points}")

    return sampled, lfc_sampled, time_points


@logger.catch
def render_curve_fitting_figure(sampled_stats: pd.DataFrame, lfc_sampled: pd.DataFrame,
                                  time_points: list[float], output_stem: Path) -> None:
    """Render multipanel figure with observed + fitted curves."""
    logger.info("Rendering curve fitting figure...")

    if sampled_stats.empty:
        logger.warning("No data to plot!")
        return

    # Apply house style before creating figure
    apply_house_style()

    n_curves = len(sampled_stats)
    n_cols = 4
    n_rows = (n_curves + n_cols - 1) // n_cols

    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_curves} panels)...")

    # Create multipanel layout
    panel_width = 180
    panel_height = 150
    multipanel = cns.multipanel(max_width=panel_width * n_cols)

    # Time points as numpy array
    x_values = np.array(time_points)
    timepoint_cols = [col for col in lfc_sampled.columns if col in ['YES0', 'YES1', 'YES2', 'YES3', 'YES4']]

    for idx, (gene_idx, row) in enumerate(sampled_stats.iterrows()):
        # Get observed values
        y_observed = lfc_sampled.loc[gene_idx, timepoint_cols].values

        # Get fitted parameters
        A, DR, DL = row['A'], row['DR'], row['DL']
        R2, RMSE = row['R2'], row['RMSE']

        # Format gene ID for title
        gene_id = "=".join(map(str, gene_idx))
        panel_label = chr(65 + idx) if idx < 26 else f"A{idx-25}"

        # Create panel
        ax = multipanel.panel(label=panel_label, width=panel_width, height=panel_height)

        # Plot observed points with ax.scatter directly
        ax.scatter(
            x_values,
            y_observed,
            s=30,
            color='#1f77b4',
            alpha=0.8,
            edgecolor='white',
            linewidths=0.5
        )

        # Plot fitted curve
        x_smooth = np.linspace(x_values.min(), x_values.max(), 100)
        y_fit = sigmoid_function(x_smooth, A, DR, DL)
        ax.plot(x_smooth, y_fit, color='#ff7f0e', linewidth=2, label='Fitted')

        # Add constraint lines
        ax.axhline(y=A, color='gray', linestyle='--', alpha=0.3, linewidth=1)
        ax.axvline(x=DL, color='gray', linestyle='--', alpha=0.3, linewidth=1)

        # Add parameter text
        param_text = f'A={A:.2f} R²={R2:.3f}\nDR={DR:.2f} RMSE={RMSE:.3f}\nDL={DL:.2f}'
        ax.text(0.05, 0.95, param_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=8)

        ax.set_title(gene_id.replace("=", " "), fontsize=9)
        ax.set_ylim(-1.5, 8.5)
        ax.set_xlabel('Time', fontsize=9)
        ax.set_ylabel('LFC', fontsize=9)
        ax.tick_params(labelsize=8)

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
