"""DESeq2 dispersion estimates figure rendering.

Input
-----
- Dispersion figure-data TSV written by ``write_dispersion_data_tsv()`` in
  ``insertion_level_depletion_analysis_has_replicates.py``: tab-separated, with
  index columns ``Chr, Coordinate, Strand, Target`` plus value columns
  ``normed_mean``, ``genewise_dispersion``, ``MAP_dispersion``,
  ``fitted_dispersion``.

Output
------
- ``<stem>.pdf`` — journal-quality vector figure (scatter rasterized).
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

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

# =============================================================================
# CONSTANTS
# =============================================================================
X_COLUMN = "normed_mean"

# Series order, legend labels and colours mirror pydeseq2's plot_dispersions():
# it passes [genewise, MAP, fitted] with labels ["Estimated", "Final", "Fitted"]
# and the matplotlib colour string "kbr" mapped positionally.
DISPERSION_SERIES: tuple[tuple[str, str, str], ...] = (
    ("genewise_dispersion", "Estimated", "k"),
    ("MAP_dispersion", "Final", "b"),
    ("fitted_dispersion", "Fitted", "r"),
)

REQUIRED_COLUMNS = (X_COLUMN, *(column for column, _, _ in DISPERSION_SERIES))

X_LABEL = "mean of normalized counts"
Y_LABEL = "dispersion"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_dispersion_data(dispersion_data_path: Path) -> pd.DataFrame:
    """Load dispersion figure-data TSV and validate schema."""
    logger.info(f"Loading dispersion data from {dispersion_data_path}...")
    df = pd.read_csv(dispersion_data_path, sep='\t', index_col=[0, 1, 2, 3])

    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"Loaded {len(df)} rows")
    return df


@logger.catch
def render_dispersion_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render the three dispersion series against normed mean on log-log axes."""
    logger.info("Rendering dispersion figure...")

    # Apply house style before creating figure
    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    ax = multipanel.panel()

    if df.empty:
        logger.warning("No data to plot!")
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel(X_LABEL)
        ax.set_ylabel(Y_LABEL)
        ax.set_title('DESeq2 dispersion estimates')
        save_dual(output_stem)
        return

    # Rasterize: ~93.6k points per series (3 series) would bloat the vector PDF.
    # Alpha is kept low because the three clouds overlap heavily — at pydeseq2's
    # own 0.5 the blue "Final" cloud paints over the black "Estimated" one.
    scatter_kws = dict(
        s=1.0,
        alpha=0.12,
        linewidths=0,
        rasterized=True,
    )

    x_values = df[X_COLUMN]
    for column, legend_label, color in DISPERSION_SERIES:
        logger.info(f"  Series {legend_label}: {column} (n={df[column].notna().sum()})")
        ax.scatter(x_values, df[column], c=color, label=legend_label, **scatter_kws)

    # Log both axes: dispersions floor at pydeseq2's min_disp of 1e-8
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_title('DESeq2 dispersion estimates')

    # Legend markers inherit the scatter's tiny size and low alpha, so scale them up
    legend = ax.legend(loc='best', frameon=False, markerscale=6)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)

    # Save dual artifacts
    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
