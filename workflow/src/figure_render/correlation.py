"""PBL-PBR correlation figure rendering.

Input
-----
- Pairs TSV with columns ``sample``, ``timepoint``, ``condition``, ``pbl``,
  ``pbr``.

Output
------
- ``<stem>.pdf`` — journal-quality vector multipanel figure.
- ``<stem>.review.png`` — screen-review raster copy.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import matplotlib.pyplot as plt
import cnsplots as cns
import numpy as np
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load pairs TSV, validate schema, filter positive values, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep='\t')

    # Validate required columns
    required_cols = ['sample', 'timepoint', 'condition', 'pbl', 'pbr']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")

    # Filter to positive values
    df_filtered = df[(df['pbl'] > 0) & (df['pbr'] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(df_filtered)} rows")

    if df_filtered.empty:
        logger.warning("No valid data points after filtering!")
        return df_filtered

    # Add log10 columns
    df_filtered['log10_pbl'] = np.log10(df_filtered['pbl'])
    df_filtered['log10_pbr'] = np.log10(df_filtered['pbr'])

    return df_filtered


@logger.catch
def render_correlation_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render multipanel correlation figure with one panel per sample-timepoint group."""
    logger.info("Rendering correlation figure...")

    # Apply house style before creating figure
    apply_house_style()

    # Group by sample and timepoint
    grouped = df.groupby(['sample', 'timepoint'], sort=True)
    n_panels = len(grouped)

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    # Create figure
    fig = cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    # Panel labels A, B, C...
    panel_labels = [chr(65 + i) for i in range(n_panels)]

    # Verified scatter kwargs for 131k points
    scatter_kws = dict(
        s=3,
        facecolor="none",
        edgecolor="gray",
        alpha=0.15,
        linewidths=0.25,
        rasterized=True
    )

    for (sample, timepoint), group_df in grouped:
        label = panel_labels.pop(0)
        logger.info(f"  Panel {label}: {sample} {timepoint} (n={len(group_df)})")

        if group_df.empty:
            # Handle empty data case
            ax = multipanel.panel(label=label)
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
            ax.set_xlabel('log$_{10}$ PBL')
            ax.set_ylabel('log$_{10}$ PBR')
            ax.set_title(f'{sample} {timepoint}')
            continue

        # Create panel and regplot
        ax = multipanel.panel(label=label)

        # regplot on log10 columns
        cns.regplot(
            data=group_df,
            x='log10_pbl',
            y='log10_pbr',
            scatter_kws=scatter_kws,
            ax=ax
        )

        # Set labels with proper subscript and title
        ax.set_xlabel('log$_{10}$ PBL')
        ax.set_ylabel('log$_{10}$ PBR')
        ax.set_title(f'{sample} {timepoint}')

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    plt.tight_layout()
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
