"""Distribution of curve fitting results figure rendering.

Input
-----
- Fitting statistics TSV with index columns and metric columns such as
  ``A``, ``DR``, ``DL``, ``t50``, ``R2``, ``RMSE``.

Output
------
- ``<stem>.pdf`` — journal-quality vector multipanel histogram figure.
- ``<stem>.review.png`` — screen-review raster copy.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_fitting_stats(fitting_stats_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Load fitting statistics and identify numeric columns for histogramming."""
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    df = pd.read_csv(fitting_stats_path, sep='\t', index_col=[0, 1, 2, 3])

    logger.info(f"Loaded {len(df)} rows")

    # Filter to successful fits
    successful = df[df['Status'] == 'Success'].copy()
    logger.info(f"Found {len(successful)} successful fits")

    # Select numeric columns of interest (exclude fitted/residual columns)
    metric_cols = ['A', 'DR', 'DL', 't10', 't50', 't90', 't_window', 't_inflection', 'y_inflection', 'auc', 'R2', 'RMSE', 'normalized_RMSE', 'AIC', 'BIC']
    available_cols = [col for col in metric_cols if col in successful.columns]

    logger.info(f"Found {len(available_cols)} metric columns: {available_cols}")

    return successful[available_cols], available_cols


@logger.catch
def render_distribution_figure(df: pd.DataFrame, metric_cols: list[str], output_stem: Path, bins: int) -> None:
    """Render multipanel histogram figure with one panel per metric."""
    logger.info("Rendering distribution figure...")

    if df.empty or not metric_cols:
        logger.warning("No data to plot!")
        return

    # Apply house style before creating figure
    apply_house_style()

    n_metrics = len(metric_cols)
    n_cols = 4
    n_rows = (n_metrics + n_cols - 1) // n_cols

    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_metrics} panels)...")

    # Create multipanel layout
    panel_width = 180
    panel_height = 150
    multipanel = cns.multipanel(max_width=panel_width * n_cols)

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_metrics)]

    for metric, label in zip(metric_cols, panel_labels):
        data = df[metric].dropna()

        logger.info(f"  Panel {label}: {metric} (n={len(data)})")

        # Create panel
        ax = multipanel.panel(label=label, width=panel_width, height=panel_height)

        if data.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(metric)
            continue

        # Histogram - use ax.hist directly since cns.histplot requires DataFrame
        ax.hist(
            data,
            bins=bins,
            alpha=0.8,
            edgecolor='white',
            linewidth=0.5
        )

        # Add statistics text
        stats_text = f'n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=8)

        ax.set_title(metric, fontsize=10)
        ax.set_xlabel('Value', fontsize=9)
        ax.set_ylabel('Frequency', fontsize=9)
        ax.tick_params(labelsize=8)

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
