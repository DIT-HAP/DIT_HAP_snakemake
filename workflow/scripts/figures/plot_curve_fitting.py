#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Curve Fitting Figure Renderer
==============================

Render observed vs fitted sigmoid curves for genes and insertions.
Automatically detects data type (gene or insertion) from the index columns,
samples N curves (default 32, deterministic seed 42), and draws observed
points + fitted sigmoid overlay in a grid layout.

For insertion-level data, panel titles include both insertion position and
the gene/interval context.

Input
-----
- Fitting statistics TSV (``-s/--stats``): curve fitting results with Status,
  model parameters (A, DR, DL), and goodness-of-fit metrics (R2, RMSE).
- LFC TSV (``-l/--lfc``): log fold change values, same index as stats file,
  one column per timepoint.
- Optional annotation TSV (``-a/--annotation``): for insertion data, maps
  insertions to gene/interval context.

Output
------
- ``<output>.pdf`` — journal-quality vector figure (scatter rasterized).
- ``<output>.review.png`` — screen-review raster copy.

Usage
-----
    python plot_curve_fitting.py -s fitting_stats.tsv -l LFC.tsv -o figures/curves
    python plot_curve_fitting.py -s fitting_stats.tsv -l LFC.tsv -o figures/curves -a annotations.tsv --verbose

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-09-02
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from loguru import logger
from matplotlib import use

use("Agg")

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from logging_setup import setup_logger  # noqa: E402
from depletion.curve_model import sigmoid_function  # noqa: E402
from figure_render.curves import render_fitted_curves_figure  # noqa: E402

import pandas as pd  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated paths and sampling parameters."""
    fitting_stats_path: Path
    lfc_path: Path
    output_stem: Path
    n_curves: int = 32
    random_seed: int = 42
    annotation_path: Path | None = None

    def __post_init__(self) -> None:
        """Validate that inputs exist and output directory can be created."""
        if not self.fitting_stats_path.exists():
            raise ValueError(f"Fitting stats file does not exist: {self.fitting_stats_path}")
        if not self.lfc_path.exists():
            raise ValueError(f"LFC file does not exist: {self.lfc_path}")
        if self.annotation_path and not self.annotation_path.exists():
            raise ValueError(f"Annotation file does not exist: {self.annotation_path}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)
        if self.n_curves < 1:
            raise ValueError(f"n_curves must be positive, got {self.n_curves}")


# =============================================================================
# CONSTANTS
# =============================================================================
# The sigmoid's parameters, in the order sigmoid_function() accepts them.
# Both gene-level and insertion-level now use (A, DR, DL) consistently
MODEL_PARAM_COLUMNS = ["A", "DR", "DL"]

# Printed in each panel's text box.
ANNOTATION_COLUMNS = ["R2", "RMSE"]

SUCCESS_STATUS = "Success"
X_LABEL = "Time"
Y_LABEL = "LFC"


# =============================================================================
# CORE LOGIC
# =============================================================================
def detect_data_type(stats_df: pd.DataFrame) -> tuple[str, list[str]]:
    """Detect whether data is gene-level or insertion-level from index names.

    Returns:
        tuple: (data_type, index_column_names)
            data_type is either "gene" or "insertion"
    """
    index_names = list(stats_df.index.names)

    # Gene-level typically has: Systematic ID, Name, FYPOviability, DeletionLibrary_essentiality
    # Insertion-level typically has: Chr, Coordinate, Strand, Target
    if "Chr" in index_names and "Coordinate" in index_names:
        logger.info("Detected insertion-level data")
        return "insertion", index_names
    elif "Systematic ID" in index_names or "Name" in index_names:
        logger.info("Detected gene-level data")
        return "gene", index_names
    else:
        logger.warning(f"Unknown data type with index names: {index_names}, defaulting to generic")
        return "generic", index_names


def load_gene_annotations(annotation_path: Path) -> pd.DataFrame:
    """Load gene annotation file mapping insertions to genes/intervals.

    Expected format: Chr, Coordinate, Strand, Target, Gene/Interval columns
    """
    logger.info(f"Loading gene annotations from {annotation_path}...")
    annotations = pd.read_csv(annotation_path, sep="\t")

    # Set index to match insertion index
    index_cols = ["Chr", "Coordinate", "Strand", "Target"]
    if all(col in annotations.columns for col in index_cols):
        annotations = annotations.set_index(index_cols)

    logger.info(f"Loaded {len(annotations)} annotation entries")
    return annotations


def enhance_insertion_titles(sampled: pd.DataFrame, annotations: pd.DataFrame | None) -> pd.DataFrame:
    """Add gene/interval context to insertion titles.

    Creates a 'title' column with format: "Chr:Coord Strand | Gene/Interval"
    """
    if annotations is None:
        # Simple title without gene context
        sampled = sampled.copy()
        sampled["title"] = sampled.apply(
            lambda row: f"{row.name[0]}:{row.name[1]} {row.name[2]}",
            axis=1
        )
        return sampled

    # Join with annotations to get gene/interval info
    sampled = sampled.copy()

    # Find gene/interval columns with priority:
    # 1. Name (gene common name)
    # 2. Systematic ID (gene systematic name)
    # 3. Feature (genomic feature type)
    # 4. Chr_Interval (interval identifier)
    exclude_cols = {"time_points", "Status", "A", "DR", "DL", "t10", "t50", "t90",
                    "t_window", "t_inflection", "y_inflection", "auc", "R2", "RMSE",
                    "normalized_RMSE", "AIC", "BIC", "Chr", "Coordinate", "Strand",
                    "Target", "Start_Interval", "End_Interval", "Transcript", "Length",
                    "Strand_Interval", "Type", "Accumulated_CDS_bases", "FYPOviability",
                    "DeletionLibrary_essentiality"}

    # Priority order for annotation columns
    priority_cols = ["Name", "Systematic ID", "Feature", "Chr_Interval"]

    gene_col = None
    for col in priority_cols:
        if col in annotations.columns:
            gene_col = col
            break

    if gene_col is None:
        # Fallback: use first non-excluded column
        available_cols = [col for col in annotations.columns if col not in exclude_cols]
        if available_cols:
            gene_col = available_cols[0]

    if gene_col is None:
        logger.warning("No suitable annotation column found")
        sampled["title"] = sampled.apply(
            lambda row: f"{row.name[0]}:{row.name[1]} {row.name[2]}",
            axis=1
        )
        return sampled

    logger.info(f"Using '{gene_col}' column for gene context")

    # Join annotations
    sampled = sampled.join(annotations[[gene_col]], how="left")

    # Create enhanced titles
    # Filter out empty, NaN, or chromosome-only values
    def format_title(row):
        pos_str = f"{row.name[0]}:{row.name[1]} {row.name[2]}"
        gene_val = row[gene_col]

        # Check if gene_val is meaningful (not empty, not NaN, not just chromosome name)
        if pd.isna(gene_val):
            return pos_str

        gene_str = str(gene_val).strip()
        if not gene_str or gene_str == row.name[0]:  # Empty or just chromosome name
            return pos_str

        return f"{pos_str} | {gene_str}"

    sampled["title"] = sampled.apply(format_title, axis=1)

    # Drop the joined column to avoid issues downstream
    sampled = sampled.drop(columns=[gene_col])

    return sampled


@logger.catch(reraise=True)
def load_and_sample_data(
    fitting_stats_path: Path,
    lfc_path: Path,
    n_curves: int,
    random_seed: int,
    annotation_path: Path | None = None,
) -> tuple[pd.DataFrame, list[float], list[str], str]:
    """Sample successful fits, join their observed LFC values, and return the x values.

    Timepoint columns are taken from the LFC table, not a name whitelist: the old
    renderer tested `col in ['YES0'..'YES4']`, which selected the wrong subset for
    any project with a different timepoint count or naming scheme.

    Returns:
        tuple: (joined_data, time_points, timepoint_columns, data_type)
    """
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")

    # Read without assuming index structure
    stats_raw = pd.read_csv(fitting_stats_path, sep="\t")

    # Detect data type and set appropriate index
    # First row tells us the column names
    non_data_cols = {"time_points", "Status", "A", "DR", "DL", "t10", "t50", "t90",
                     "t_window", "t_inflection", "y_inflection", "auc", "R2", "RMSE",
                     "normalized_RMSE", "AIC", "BIC"}

    # Index columns are those before the data columns
    potential_index_cols = []
    for col in stats_raw.columns:
        if col in non_data_cols:
            break
        potential_index_cols.append(col)

    if not potential_index_cols:
        raise ValueError("Could not identify index columns in fitting statistics")

    logger.info(f"Using index columns: {potential_index_cols}")
    fitting_stats = stats_raw.set_index(potential_index_cols)

    logger.info(f"Loaded {len(fitting_stats)} rows")

    # Detect data type
    data_type, index_names = detect_data_type(fitting_stats)

    successful = fitting_stats[fitting_stats["Status"] == SUCCESS_STATUS].copy()
    logger.info(f"Found {len(successful)} successful fits")

    if successful.empty:
        logger.warning("No successful fits found!")
        return successful, [], [], data_type

    sampled = successful.sample(n=min(n_curves, len(successful)), random_state=random_seed)
    logger.info(f"Sampled {len(sampled)} curves for plotting")

    # Enhance titles for insertion data
    if data_type == "insertion" and annotation_path:
        annotations = load_gene_annotations(annotation_path)
        sampled = enhance_insertion_titles(sampled, annotations)
    elif data_type == "insertion":
        sampled = enhance_insertion_titles(sampled, None)

    logger.info(f"Loading LFC data from {lfc_path}...")
    lfc_raw = pd.read_csv(lfc_path, sep="\t")
    lfc_data = lfc_raw.set_index(potential_index_cols)

    timepoint_columns = list(lfc_data.columns)
    time_points = [float(t) for t in sampled["time_points"].iloc[0].split(",")]
    logger.info(f"Time points: {time_points} for columns {timepoint_columns}")

    # The stats table's own timepoint columns collide with the LFC ones, so the
    # LFC values are joined under a suffix and then renamed back.
    lfc_sampled = lfc_data.loc[sampled.index, timepoint_columns]
    joined = sampled.drop(columns=timepoint_columns, errors="ignore").join(lfc_sampled)

    return joined, time_points, timepoint_columns, data_type


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Render curve fitting figure from stats and LFC data. "
                    "Automatically detects gene-level or insertion-level data."
    )
    parser.add_argument("-s", "--stats", type=Path, required=True,
                       help="Input fitting statistics TSV file")
    parser.add_argument("-l", "--lfc", type=Path, required=True,
                       help="Input LFC TSV file")
    parser.add_argument("-o", "--output", type=Path, required=True,
                       help="Output file stem (extension will be added)")
    parser.add_argument("-a", "--annotation", type=Path, required=False,
                       help="Optional annotation file for insertion-to-gene mapping")
    parser.add_argument("-n", "--n-curves", type=int, default=32,
                       help="Number of curves to sample (default: 32)")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for sampling (default: 42)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load data, sample curves, render figure, and save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    # Strip extension from output if provided
    output_stem = args.output.with_suffix('')

    # Validate paths
    try:
        config = PlotConfig(
            fitting_stats_path=args.stats,
            lfc_path=args.lfc,
            output_stem=output_stem,
            n_curves=args.n_curves,
            random_seed=args.seed,
            annotation_path=args.annotation,
        )
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1

    logger.info("=== Curve Fitting Figure Rendering ===")

    try:
        joined, time_points, timepoint_columns, data_type = load_and_sample_data(
            config.fitting_stats_path,
            config.lfc_path,
            config.n_curves,
            config.random_seed,
            config.annotation_path,
        )

        if joined.empty:
            logger.warning("No data to plot")
            return 0

        # Use custom title column for insertion data if available
        title_col = "title" if "title" in joined.columns else None

        render_fitted_curves_figure(
            joined,
            config.output_stem,
            x_values=time_points,
            value_columns=timepoint_columns,
            model=sigmoid_function,
            model_params=MODEL_PARAM_COLUMNS,
            annotations=ANNOTATION_COLUMNS,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            title_column=title_col,
        )

    except Exception as e:
        logger.error(f"Error during rendering: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
