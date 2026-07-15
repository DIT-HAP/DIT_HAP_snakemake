#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
#     "matplotlib",
#     "numpy",
#     "pandas",
# ]
# ///
"""
Insertion Density Analysis for Transposon Insertion Sequencing
==============================================================

Analyzes insertion density patterns in transposon sequencing data, comparing
the initial and final timepoints to reveal the classic piggyBac phenotyping
signal: insertion sites well-covered at the initial timepoint that drop out
by the final timepoint mark depleted (typically essential) genes.

Loads raw insertion read counts (summed across biological replicates) and
genomic annotations, filters for in-gene insertions using established
criteria (non-intergenic and distance to stop codon > 4), and computes
per-gene density metrics at both timepoints plus their change.

For each gene it derives insertion-site density per kilobase at the initial
and final timepoint (and their difference / log2 fold-change), gap-distribution
metrics (including a Gini coefficient of insertion location), read-depth
inequality (Gini coefficient of depth) at both timepoints, and strand-preference
measures. A multi-page PDF of histogram distributions and initial-vs-final
scatter plots is produced alongside the density statistics table.

Input
-----
- Insertion data (``-i``): tab-separated file with a 4-level row index and a
  2-level column header (Sample, Timepoint). Read counts are summed across
  samples (replicates) at the initial (``-t``) and final (``-f``) timepoints.
- Annotation data (``-a``): tab-separated file with a 4-level index and gene
  annotation columns (Type, Distance_to_stop_codon, Systematic ID, Name,
  FYPOviability, Chr_Interval, Strand_Interval, ParentalRegion_start,
  ParentalRegion_end, ParentalRegion_length, Insertion_direction).

Output
------
- Density statistics table (``-o``): tab-separated, one row per gene.
- Histogram/scatter distribution report: multi-page PDF written next to the
  table as ``<output_stem>_histograms.pdf``.

Usage
-----
    python insertion_density_analysis.py -i raw_reads.filtered.tsv -a annotations.tsv -o density_stats.tsv -t YES0 -f YES4

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
import time
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
# The following is for plotting
SCRIPT_DIR = Path(__file__).parent.resolve()
plt.style.use(SCRIPT_DIR / "../../../config/DIT_HAP.mplstyle")
AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths and the initial/final timepoint column names."""

    insertion_data_path: Path
    annotations_path: Path
    output_path: Path
    initial_timepoint: str
    final_timepoint: str

    def __post_init__(self) -> None:
        for path in (self.insertion_data_path, self.annotations_path):
            if not path.exists():
                raise ValueError(f"Input file does not exist: {path}")
            if not path.is_file():
                raise ValueError(f"Input path is not a file: {path}")
            if path.suffix.lower() not in ['.csv', '.tsv', '.txt']:
                raise ValueError(f"Input file must be CSV or TSV format. Got: {path.suffix}")
            if path.stat().st_size == 0:
                raise ValueError(f"Input file is empty: {path}")

        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        if self.output_path.suffix.lower() not in ['.csv', '.tsv', '.txt']:
            raise ValueError(f"Output file must be CSV or TSV format. Got: {self.output_path.suffix}")

        if not self.initial_timepoint:
            raise ValueError("initial_timepoint must be a non-empty string")
        if not self.final_timepoint:
            raise ValueError("final_timepoint must be a non-empty string")


@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Aggregate summary metrics for the analyzed gene set."""

    total_genes_analyzed: int
    total_insertions_analyzed: int
    mean_insertion_density_per_kb_initial: float
    mean_insertion_density_per_kb_final: float
    mean_insertion_density_log2fc: float
    mean_gini_coefficient_of_depth_initial: float
    mean_strand_bias: float


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
def load_insertion_data(
    insertion_data_path: Path, initial_timepoint: str, final_timepoint: str
) -> pd.DataFrame:
    """Load raw insertion read counts and sum replicate samples at the initial/final timepoints."""
    logger.info(f"Loading insertion data from {insertion_data_path}")

    try:
        raw_counts = pd.read_csv(
            insertion_data_path, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t"
        )

        if not isinstance(raw_counts.columns, pd.MultiIndex):
            raise ValueError("Expected a 2-level (Sample, Timepoint) column header")

        available_timepoints = raw_counts.columns.get_level_values(1).unique().tolist()
        for timepoint in (initial_timepoint, final_timepoint):
            if timepoint not in available_timepoints:
                raise ValueError(
                    f"Timepoint '{timepoint}' not found. Available timepoints: {available_timepoints}"
                )

        # Sum read counts across replicate samples at each timepoint
        initial_counts = raw_counts.xs(initial_timepoint, level=1, axis=1).sum(axis=1, skipna=True)
        final_counts = raw_counts.xs(final_timepoint, level=1, axis=1).sum(axis=1, skipna=True)

        insertion_data = pd.DataFrame({
            "reads_initial": initial_counts,
            "reads_final": final_counts,
        })

        logger.info(f"Loaded {len(insertion_data)} insertions with read count data")
        logger.info(
            f"Initial ({initial_timepoint}) reads - Mean: {insertion_data['reads_initial'].mean():.2f}, "
            f"Median: {insertion_data['reads_initial'].median():.2f}"
        )
        logger.info(
            f"Final ({final_timepoint}) reads - Mean: {insertion_data['reads_final'].mean():.2f}, "
            f"Median: {insertion_data['reads_final'].median():.2f}"
        )

        return insertion_data

    except Exception as e:
        raise ValueError(f"Error loading insertion data: {e}")


@logger.catch
def load_annotation_data(annotations_path: Path) -> pd.DataFrame:
    """Load genomic annotations for insertions."""
    logger.info(f"Loading annotation data from {annotations_path}")

    try:
        required_cols = [
            'Type', 'Distance_to_stop_codon', 'Systematic ID', 'Name', 'FYPOviability',
            'Chr_Interval', 'Strand_Interval', 'ParentalRegion_start',
            'ParentalRegion_end', 'ParentalRegion_length', 'Insertion_direction',
        ]

        annotations = pd.read_csv(annotations_path, index_col=[0, 1, 2, 3], sep="\t")

        missing_cols = [col for col in required_cols if col not in annotations.columns]
        if missing_cols:
            raise ValueError(f"Missing required annotation columns: {missing_cols}")

        empty_systematic_ids = annotations['Systematic ID'].isna().sum()
        if empty_systematic_ids > 0:
            logger.warning(f"Found {empty_systematic_ids} annotations with empty Systematic ID")

        invalid_lengths = (annotations['ParentalRegion_length'] <= 0).sum()
        if invalid_lengths > 0:
            logger.warning(f"Found {invalid_lengths} annotations with invalid gene length (<= 0)")

        logger.info(f"Loaded annotations for {len(annotations)} insertions")
        logger.info(f"Unique gene types: {annotations['Type'].value_counts().to_dict()}")

        return annotations

    except Exception as e:
        raise ValueError(f"Error loading annotation data: {e}")


@logger.catch
def filter_in_gene_insertions(insertion_data: pd.DataFrame,
                             annotations: pd.DataFrame) -> pd.DataFrame:
    """Filter insertions to include only those within genes using established criteria."""
    logger.info("Filtering for in-gene insertions")

    merged_data = pd.merge(
        insertion_data, annotations,
        left_index=True, right_index=True,
        how='inner'
    )

    # Apply in-gene filtering criteria (same as gene_level_depletion_analysis.py)
    in_gene_mask = (
        (merged_data['Type'] != 'Intergenic region') &
        (merged_data['Distance_to_stop_codon'] > 4)
    )

    in_gene_insertions = merged_data[in_gene_mask].copy()

    logger.info(f"Found {len(in_gene_insertions)} in-gene insertions")
    logger.info(f"Filtered out {len(merged_data) - len(in_gene_insertions)} intergenic/near-stop insertions")

    return in_gene_insertions


@logger.catch
def calculate_insertion_statistics(gene_insertions: pd.DataFrame) -> dict[str, int | float]:
    """Calculate insertion-site density at the initial and final timepoint for a gene.

    A site counts as "detected" at a timepoint when its summed replicate read
    count is > 0. Since the input is already hard-filtered on the initial
    timepoint, every retained row is detected at the initial timepoint by
    construction; sites that lose signal by the final timepoint (density
    drop) are the classic piggyBac phenotyping readout for essential genes.
    """
    coordinates = gene_insertions.index.get_level_values(1)
    total_insertions = len(gene_insertions)
    unique_sites_initial = len(coordinates[gene_insertions["reads_initial"] > 0].unique())
    unique_sites_final = len(coordinates[gene_insertions["reads_final"] > 0].unique())

    gene_length = gene_insertions["ParentalRegion_length"].iloc[0]
    density_initial = (unique_sites_initial / gene_length) * 1000 if gene_length > 0 else 0
    density_final = (unique_sites_final / gene_length) * 1000 if gene_length > 0 else 0

    # Gene length cancels out in the ratio, so the fold-change of per-kb
    # density equals the fold-change of raw site counts (+1 pseudocount).
    density_log2fc = np.log2((unique_sites_final + 1) / (unique_sites_initial + 1))

    return {
        'total_insertions': total_insertions,
        'unique_sites_initial': unique_sites_initial,
        'unique_sites_final': unique_sites_final,
        'gene_length': gene_length,
        'insertion_density_per_kb_initial': round(density_initial, 3),
        'insertion_density_per_kb_final': round(density_final, 3),
        'insertion_density_change': round(density_final - density_initial, 3),
        'insertion_density_log2fc': round(density_log2fc, 3),
    }


def calculate_gap_statistics(gene_insertions: pd.DataFrame) -> dict[str, int | float | str]:
    """Calculate statistics about gaps between insertions within a gene."""
    coordinates = sorted(gene_insertions.index.get_level_values(1).unique())
    start_coordinate = gene_insertions["ParentalRegion_start"].iloc[0]
    end_coordinate = gene_insertions["ParentalRegion_end"].iloc[0]
    gene_length = gene_insertions["ParentalRegion_length"].iloc[0]

    coordinates_with_start_and_end = sorted(set([start_coordinate] + coordinates + [end_coordinate]))

    # Calculate gaps between consecutive insertions
    gaps = [coordinates_with_start_and_end[i+1] - coordinates_with_start_and_end[i] - 1 for i in range(len(coordinates_with_start_and_end)-1)]
    gaps = [gap for gap in gaps if gap > 0]  # Only count actual gaps

    normalized_gaps = [round(gap / gene_length, 3) for gap in gaps]
    gini_coefficient_of_location = calculate_gini_coefficient(normalized_gaps)

    if not gaps:
        return {
            'num_gaps': 0,
            'largest_gap': 0,
            'largest_gap_fraction': 0,
            'smallest_gap': 0,
            'smallest_gap_fraction': 0,
            'mean_gap_length': 0,
            'mean_gap_length_fraction': 0,
            'median_gap_length': 0,
            'median_gap_length_fraction': 0,
            'gap_length_sd': 0,
            'gap_length_sd_fraction': 0,
            'all_gap_lengths': "",
            'all_gap_lengths_fraction': "",
            'gini_coefficient_of_location': np.nan
        }

    return {
        'num_gaps': len(gaps),
        'largest_gap': max(gaps),
        'largest_gap_fraction': max(normalized_gaps),
        'smallest_gap': min(gaps),
        'smallest_gap_fraction': min(normalized_gaps),
        'mean_gap_length': round(np.mean(gaps), 2),
        'mean_gap_length_fraction': round(np.mean(normalized_gaps), 2),
        'median_gap_length': round(np.median(gaps), 2),
        'median_gap_length_fraction': round(np.median(normalized_gaps), 2),
        'gap_length_sd': round(np.std(gaps), 2),
        'gap_length_sd_fraction': round(np.std(normalized_gaps), 2),
        'all_gap_lengths': ",".join(map(str, gaps)),
        'all_gap_lengths_fraction': ",".join(map(str, normalized_gaps)),
        'gini_coefficient_of_location': round(gini_coefficient_of_location, 3)
    }


def calculate_gini_coefficient(values: np.ndarray) -> float:
    """Calculate Gini coefficient to measure inequality in read distribution."""
    if len(values) == 0:
        return 0.0

    # Sort values
    sorted_values = np.sort(values)
    n = len(sorted_values)

    # A gene fully depleted by the final timepoint has all-zero read counts;
    # the Gini ratio is 0/0 there, not an inequality signal, so short-circuit.
    cumsum = np.cumsum(sorted_values)
    if cumsum[-1] == 0:
        return 0.0

    # Calculate Gini coefficient
    gini = (2 * np.sum((np.arange(1, n+1) * sorted_values))) / (n * cumsum[-1]) - (n + 1) / n

    return max(0.0, min(1.0, gini))  # Ensure result is between 0 and 1


def calculate_read_statistics(gene_insertions: pd.DataFrame) -> dict[str, int | float]:
    """Calculate read distribution statistics at the initial and final timepoint for a gene."""
    initial_counts = gene_insertions["reads_initial"].values
    final_counts = gene_insertions["reads_final"].values

    if len(initial_counts) == 0:
        return {
            'total_reads_initial': 0,
            'total_reads_final': 0,
            'mean_reads_per_insertion_initial': 0,
            'mean_reads_per_insertion_final': 0,
            'gini_coefficient_of_depth_initial': 0,
            'gini_coefficient_of_depth_final': 0,
            'total_reads_log2fc': 0,
        }

    total_reads_initial = initial_counts.sum()
    total_reads_final = final_counts.sum()

    return {
        'total_reads_initial': int(total_reads_initial),
        'total_reads_final': int(total_reads_final),
        'mean_reads_per_insertion_initial': round(np.mean(initial_counts), 2),
        'mean_reads_per_insertion_final': round(np.mean(final_counts), 2),
        'gini_coefficient_of_depth_initial': round(calculate_gini_coefficient(initial_counts), 3),
        'gini_coefficient_of_depth_final': round(calculate_gini_coefficient(final_counts), 3),
        'total_reads_log2fc': round(np.log2((total_reads_final + 1) / (total_reads_initial + 1)), 3),
    }


def calculate_strand_statistics(gene_insertions: pd.DataFrame) -> dict[str, int | float]:
    """Calculate strand preference and pairing statistics."""
    strands = gene_insertions["Insertion_direction"].values
    coordinates = gene_insertions.index.get_level_values(1)

    # Count forward and reverse insertions
    forward_count = (strands == 'Forward').sum()
    reverse_count = (strands == 'Reverse').sum()
    total_insertions = len(strands)
    total_sites = len(coordinates.unique())

    # Calculate strand preference
    forward_preference = forward_count / total_insertions if total_insertions > 0 else 0
    reverse_preference = reverse_count / total_insertions if total_insertions > 0 else 0

    # Calculate strand bias (absolute difference from 50:50)
    strand_bias = abs(forward_preference - 0.5)

    # Count paired insertions (same coordinate, different strands)
    coord_strand_pairs = list(zip(coordinates, strands))
    coord_counts = defaultdict(lambda: {'forward': 0, 'reverse': 0})

    for coord, strand in coord_strand_pairs:
        if strand == 'Forward':
            coord_counts[coord]['forward'] += 1
        else:
            coord_counts[coord]['reverse'] += 1

    # Count sites with both forward and reverse insertions
    paired_sites = sum(1 for counts in coord_counts.values()
                      if counts['forward'] > 0 and counts['reverse'] > 0)

    paired_sites_fraction = paired_sites / total_sites if total_sites > 0 else 0

    return {
        'forward_insertions': forward_count,
        'reverse_insertions': reverse_count,
        'forward_preference': round(forward_preference, 3),
        'reverse_preference': round(reverse_preference, 3),
        'strand_bias': round(strand_bias, 3),
        'paired_sites': paired_sites,
        'paired_sites_fraction': round(paired_sites_fraction, 3)
    }


@logger.catch
def analyze_gene_insertions(gene_id: str, gene_insertions: pd.DataFrame) -> dict[str, str | int | float]:
    """Perform comprehensive analysis of insertions within a single gene."""
    # Calculate all statistics
    insertion_stats = calculate_insertion_statistics(gene_insertions)
    gap_stats = calculate_gap_statistics(gene_insertions)
    read_stats = calculate_read_statistics(gene_insertions)
    strand_stats = calculate_strand_statistics(gene_insertions)

    # Combine all statistics
    gene_analysis = {
        'Systematic ID': gene_id,
        'Name': gene_insertions['Name'].iloc[0],
        'Chr': gene_insertions['Chr_Interval'].iloc[0],
        'Start': gene_insertions['ParentalRegion_start'].iloc[0],
        'End': gene_insertions['ParentalRegion_end'].iloc[0],
        'Length': gene_insertions['ParentalRegion_length'].iloc[0],
        'Strand': gene_insertions['Strand_Interval'].iloc[0],
        'FYPOviability': gene_insertions['FYPOviability'].iloc[0],
    }

    gene_analysis.update(insertion_stats)
    gene_analysis.update(gap_stats)
    gene_analysis.update(read_stats)
    gene_analysis.update(strand_stats)

    return gene_analysis


def generate_summary_statistics(results_df: pd.DataFrame) -> dict[str, int | float]:
    """Generate summary statistics across all analyzed genes."""
    stats = {
        'total_genes_analyzed': len(results_df),
        'total_insertions_analyzed': results_df['total_insertions'].sum(),
        'total_unique_sites_initial': results_df['unique_sites_initial'].sum(),
        'total_unique_sites_final': results_df['unique_sites_final'].sum(),
        'mean_insertions_per_gene': results_df['total_insertions'].mean(),
        'median_insertions_per_gene': results_df['total_insertions'].median(),
        'mean_insertion_density_per_kb_initial': results_df['insertion_density_per_kb_initial'].mean(),
        'mean_insertion_density_per_kb_final': results_df['insertion_density_per_kb_final'].mean(),
        'median_insertion_density_per_kb_initial': results_df['insertion_density_per_kb_initial'].median(),
        'median_insertion_density_per_kb_final': results_df['insertion_density_per_kb_final'].median(),
        'mean_insertion_density_log2fc': results_df['insertion_density_log2fc'].mean(),
        'genes_with_density_drop': len(results_df[results_df['insertion_density_log2fc'] < -0.5]),
        'mean_gini_coefficient_of_location': results_df['gini_coefficient_of_location'].mean(),
        'genes_with_high_inequality_of_location': len(results_df[results_df['gini_coefficient_of_location'] > 0.5]),
        'mean_gini_coefficient_of_depth_initial': results_df['gini_coefficient_of_depth_initial'].mean(),
        'genes_with_high_inequality_of_depth': len(results_df[results_df['gini_coefficient_of_depth_initial'] > 0.5]),
        'mean_strand_bias': results_df['strand_bias'].mean(),
        'genes_with_strong_strand_bias': len(results_df[results_df['strand_bias'] > 0.2]),
        'mean_paired_sites_fraction': results_df['paired_sites_fraction'].mean(),
        'genes_with_high_paired_sites_fraction': len(results_df[results_df['paired_sites_fraction'] > 0.5])
    }

    return stats


def plot_initial_vs_final_scatter(
    results_df: pd.DataFrame, initial_timepoint: str, final_timepoint: str, pdf: PdfPages
) -> None:
    """Plot initial-vs-final insertion density scatter panels — the classic phenotyping view."""
    fig, axes = plt.subplots(2, 2, figsize=(AX_WIDTH * 2, AX_HEIGHT * 2))

    essentiality = None
    if 'FYPOviability' in results_df.columns:
        essentiality = results_df['FYPOviability']

    def scatter_by_essentiality(ax, x, y, xlabel, ylabel, title):
        """Draw a scatter plot, colored by FYPOviability when available."""
        if essentiality is not None:
            for label, color in [('viable', COLORS[0]), ('inviable', COLORS[3])]:
                mask = essentiality == label
                ax.scatter(x[mask], y[mask], s=6, alpha=0.4, color=color,
                           edgecolors='none', label=label, rasterized=True)
            ax.legend(fontsize=8, loc='best')
        else:
            ax.scatter(x, y, s=6, alpha=0.4, color=COLORS[0], edgecolors='none', rasterized=True)

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='semibold')
        ax.grid(True, alpha=0.3)

    # Panel 1: initial vs. final density, with the y = x reference line
    density_initial = results_df['insertion_density_per_kb_initial']
    density_final = results_df['insertion_density_per_kb_final']
    scatter_by_essentiality(
        axes[0, 0], density_initial, density_final,
        f'Insertion density per kb ({initial_timepoint})',
        f'Insertion density per kb ({final_timepoint})',
        'Initial vs. Final Insertion Density',
    )
    max_val = max(density_initial.max(), density_final.max(), 1)
    axes[0, 0].plot([0, max_val], [0, max_val], color='red', linestyle='--', alpha=0.6, linewidth=1)

    # Panel 2: initial density vs. log2 fold-change, with the y = 0 reference line
    log2fc = results_df['insertion_density_log2fc']
    scatter_by_essentiality(
        axes[0, 1], density_initial, log2fc,
        f'Insertion density per kb ({initial_timepoint})',
        f'log2FC density ({final_timepoint} / {initial_timepoint})',
        'Density Depletion vs. Initial Coverage',
    )
    axes[0, 1].axhline(0, color='red', linestyle='--', alpha=0.6, linewidth=1)

    # Panel 3: initial vs. final total reads (log scale)
    reads_initial = results_df['total_reads_initial']
    reads_final = results_df['total_reads_final']
    scatter_by_essentiality(
        axes[1, 0], reads_initial, reads_final,
        f'Total reads ({initial_timepoint})',
        f'Total reads ({final_timepoint})',
        'Initial vs. Final Read Depth',
    )
    axes[1, 0].set_xscale('symlog')
    axes[1, 0].set_yscale('symlog')

    # Panel 4: initial vs. final Gini coefficient of read depth
    gini_initial = results_df['gini_coefficient_of_depth_initial']
    gini_final = results_df['gini_coefficient_of_depth_final']
    scatter_by_essentiality(
        axes[1, 1], gini_initial, gini_final,
        f'Gini coefficient of depth ({initial_timepoint})',
        f'Gini coefficient of depth ({final_timepoint})',
        'Initial vs. Final Depth Inequality',
    )
    axes[1, 1].plot([0, 1], [0, 1], color='red', linestyle='--', alpha=0.6, linewidth=1)

    fig.suptitle('Initial vs. Final Timepoint Comparison', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def plot_numeric_distributions_to_pdf(
    results_df: pd.DataFrame, output_path: Path, initial_timepoint: str, final_timepoint: str
) -> None:
    """Generate histograms and scatter comparisons and save to a multi-page PDF."""
    logger.info("Generating histograms for numeric columns in PDF format")

    # Identify numeric columns (excluding string columns like gene names)
    numeric_columns = results_df.select_dtypes(include=[np.number]).columns.tolist()

    # Remove columns that are identifiers or coordinates
    exclude_patterns = ['Start', 'End', 'Length', 'Chr']
    numeric_columns = [col for col in numeric_columns
                      if not any(pattern in col for pattern in exclude_patterns)]

    if not numeric_columns:
        logger.warning("No numeric columns found for plotting")
        return

    # Group columns by category for better organization
    column_groups = {
        'Density Metrics': [col for col in numeric_columns if 'density' in col.lower()],
        'Gap Statistics': [col for col in numeric_columns if 'gap' in col.lower()],
        'Read Statistics': [col for col in numeric_columns if any(term in col.lower()
                           for term in ['reads', 'gini_coefficient_of_depth'])],
        'Strand Statistics': [col for col in numeric_columns if any(term in col.lower()
                             for term in ['forward', 'reverse', 'strand', 'paired'])],
        'Location Statistics': [col for col in numeric_columns if 'gini_coefficient_of_location' in col.lower()],
        'Count Statistics': [col for col in numeric_columns if any(term in col.lower()
                            for term in ['total_insertions', 'unique_sites', 'num_gaps'])]
    }

    # Add any remaining columns to a general category
    all_categorized = set()
    for group_cols in column_groups.values():
        all_categorized.update(group_cols)

    remaining_cols = [col for col in numeric_columns if col not in all_categorized]
    if remaining_cols:
        column_groups['Other Statistics'] = remaining_cols

    # Create PDF file
    pdf_path = output_path.parent / f"{output_path.stem}_histograms.pdf"

    COLOR_PALETTE = COLORS

    with PdfPages(pdf_path) as pdf:
        # Create title page
        create_title_page(pdf, results_df)

        # Create summary plot with key metrics first
        create_summary_histogram_plot_pdf(results_df, pdf)

        # Initial-vs-final scatter comparison — the classic phenotyping view
        plot_initial_vs_final_scatter(results_df, initial_timepoint, final_timepoint, pdf)

        # Generate plots for each group
        for group_name, group_columns in column_groups.items():
            if not group_columns:
                continue

            logger.info(f"Plotting {group_name}: {len(group_columns)} columns")

            # Calculate subplot layout
            n_cols = min(3, len(group_columns))
            n_rows = int(np.ceil(len(group_columns) / n_cols))

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
            if n_rows == 1 and n_cols == 1:
                axes = [axes]
            elif n_rows == 1 or n_cols == 1:
                axes = axes.flatten()
            else:
                axes = axes.flatten()

            for idx, column in enumerate(group_columns):
                ax = axes[idx] if len(group_columns) > 1 else axes[0]

                # Get data and remove NaN values
                data = results_df[column].dropna()

                if len(data) == 0:
                    ax.text(0.5, 0.5, 'No Data', transform=ax.transAxes,
                           ha='center', va='center', fontsize=12)
                    ax.set_title(column, fontsize=11, fontweight='semibold')
                    continue

                # Choose color based on column type
                if 'density' in column.lower():
                    color = COLOR_PALETTE[0]
                elif 'gap' in column.lower():
                    color = COLOR_PALETTE[1]
                elif any(term in column.lower() for term in ['reads', 'gini']):
                    color = COLOR_PALETTE[2]
                else:
                    color = COLOR_PALETTE[3]

                # Check if this is a count/density metric that should be log transformed
                # Include: reads, counts, insertion sites, insertion density
                # Exclude: gini coefficients (bounded [0,1]), percentages, ratios
                is_count_metric = any(term in column.lower() for term in
                                     ['reads', 'basemean', 'count', 'sites', 'density'])
                is_bounded_metric = any(term in column.lower() for term in
                                       ['gini', 'percent', 'fraction', 'bias'])
                # log2fc / change columns are already on a difference scale — never log-transform those
                is_diff_metric = any(term in column.lower() for term in ['log2fc', '_change'])

                # Transform data if it's count-like but not bounded or a difference metric
                # Always use log10(data + 1) to handle zeros consistently across initial/final pairs
                if is_count_metric and not is_bounded_metric and not is_diff_metric:
                    # Log transform with +1 pseudo-count (handles zeros consistently)
                    plot_data = np.log10(data + 1)
                    xlabel = 'log10(Value + 1)'

                    # Calculate statistics on original data for annotation
                    mean_val = data.mean()
                    median_val = data.median()

                    # Calculate mean and median of transformed data for reference lines
                    plot_mean = plot_data.mean()
                    plot_median = plot_data.median()
                else:
                    # Use original data for bounded metrics, percentages, and difference metrics
                    plot_data = data
                    xlabel = 'Value'
                    mean_val = data.mean()
                    median_val = data.median()
                    plot_mean = mean_val
                    plot_median = median_val

                # Create histogram
                n_bins = min(30, max(10, int(np.sqrt(len(plot_data)))))
                ax.hist(plot_data, bins=n_bins, color=color, alpha=0.7,
                       edgecolor='white', linewidth=0.5)

                # Add statistics text (always show original data statistics)
                stats_text = f'Mean: {mean_val:.3f}\nMedian: {median_val:.3f}\nN: {len(data)}'

                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                       fontsize=8, verticalalignment='top',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

                # Formatting
                ax.set_title(column.replace('_', ' ').title(), fontsize=11, fontweight='semibold')
                ax.set_xlabel(xlabel, fontsize=10)
                ax.set_ylabel('Frequency', fontsize=10)
                ax.grid(True, alpha=0.3)

                # Add vertical lines for mean and median (on transformed scale if applicable)
                ax.axvline(plot_mean, color='red', linestyle='--', alpha=0.7, linewidth=1)
                ax.axvline(plot_median, color='orange', linestyle=':', alpha=0.7, linewidth=1)

            # Hide unused subplots
            for idx in range(len(group_columns), len(axes)):
                axes[idx].set_visible(False)

            # Add group title
            fig.suptitle(f'{group_name} - Distribution Analysis',
                        fontsize=16, fontweight='bold', y=0.98)
            plt.tight_layout()

            # Save page to PDF
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    logger.info(f"Multi-page histogram PDF saved to {pdf_path}")


def create_title_page(pdf: PdfPages, results_df: pd.DataFrame) -> None:
    """Create a title page for the PDF with analysis summary."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')

    # Title
    ax.text(0.5, 0.9, 'Insertion Density Analysis',
           transform=ax.transAxes, fontsize=24, fontweight='bold',
           ha='center', va='center')

    ax.text(0.5, 0.85, 'Initial vs. Final Timepoint Comparison Report',
           transform=ax.transAxes, fontsize=16,
           ha='center', va='center', style='italic')

    # Analysis summary
    summary_text = f"""
Analysis Summary:
• Total genes analyzed: {len(results_df):,}
• Numeric metrics calculated: {len(results_df.select_dtypes(include=[np.number]).columns)}
• Analysis includes: initial/final insertion density, gap statistics, read
  distribution, and strand preferences

Report Contents:
1. Key Metrics Summary (6 most important measures)
2. Initial vs. Final Scatter Comparison (density, reads, depth inequality)
3. Density Metrics (initial/final density per kb and log2 fold-change)
4. Gap Statistics (gap lengths, counts, and distributions)
5. Read Statistics (read counts and depth inequality, initial/final)
6. Strand Statistics (forward/reverse preferences and pairing)
7. Location Statistics (spatial inequality measures)
8. Count Statistics (total insertions, unique sites, gap counts)

Statistical Annotations:
• Red dashed line: Mean value / y=x or y=0 reference (scatter panels)
• Orange dotted line: Median value
• Text box: Mean, Median, and sample size (N)
• Bin count optimized using square root rule
• Read depth metrics: log10 transformed for better visualization
• Statistics shown: Original data values (not transformed)
"""

    ax.text(0.1, 0.7, summary_text, transform=ax.transAxes, fontsize=12,
           va='top', ha='left', linespacing=1.5)

    # Footer
    ax.text(0.5, 0.1, f'Generated on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
           transform=ax.transAxes, fontsize=10, ha='center', va='center',
           style='italic', color='gray')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def create_summary_histogram_plot_pdf(results_df: pd.DataFrame, pdf: PdfPages) -> None:
    """Create a summary plot with the most important metrics for PDF."""
    # Select key metrics for summary plot
    key_metrics = [
        'insertion_density_per_kb_initial',
        'insertion_density_per_kb_final',
        'insertion_density_log2fc',
        'gini_coefficient_of_depth_initial',
        'gini_coefficient_of_location',
        'strand_bias',
    ]

    COLOR_PALETTE = COLORS

    # Filter to only include columns that exist in the data
    available_metrics = [col for col in key_metrics if col in results_df.columns]

    if not available_metrics:
        logger.warning("No key metrics available for summary plot")
        return

    # Create summary plot
    n_cols = 3
    n_rows = int(np.ceil(len(available_metrics) / n_cols))

    plot_width, plot_height = plt.rcParams['figure.figsize']
    fig_width = plot_width * n_cols
    fig_height = plot_height * n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
    if n_rows == 1:
        axes = axes.flatten() if n_cols > 1 else [axes]
    else:
        axes = axes.flatten()

    colors = [COLOR_PALETTE[0], COLOR_PALETTE[1],
              COLOR_PALETTE[2], COLOR_PALETTE[3]] * 2

    for idx, metric in enumerate(available_metrics):
        ax = axes[idx]
        data = results_df[metric].dropna()

        if len(data) == 0:
            ax.text(0.5, 0.5, 'No Data', transform=ax.transAxes,
                   ha='center', va='center')
            ax.set_title(metric.replace('_', ' ').title())
            continue

        # Check if this is a read depth related metric that should be log transformed
        is_read_depth = any(term in metric.lower() for term in
                          ['reads', 'basemean', 'count', 'depth', 'gini_coefficient_of_depth'])
        is_diff_metric = any(term in metric.lower() for term in ['log2fc', '_change'])

        # Transform data if it's read depth related
        if is_read_depth and not is_diff_metric and data.min() > 0:
            # Log transform the data values (add small constant to handle zeros)
            plot_data = np.log10(data + 1)
            xlabel = 'log10(Value + 1)'

            # Calculate statistics on original data for annotation
            mean_val = data.mean()
            median_val = data.median()
            std_val = data.std()

            # Calculate mean and median of transformed data for reference lines
            plot_mean = plot_data.mean()
            plot_median = plot_data.median()
        else:
            # Use original data for other metrics
            plot_data = data
            xlabel = 'Value'
            mean_val = data.mean()
            median_val = data.median()
            std_val = data.std()
            plot_mean = mean_val
            plot_median = median_val

        # Create histogram
        n_bins = min(25, max(10, int(np.sqrt(len(plot_data)))))
        ax.hist(plot_data, bins=n_bins, color=colors[idx], alpha=0.7,
               edgecolor='white', linewidth=0.5)

        # Add statistics (always show original data statistics)
        stats_text = f'Mean: {mean_val:.3f}\nMedian: {median_val:.3f}\nStd: {std_val:.3f}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top')

        # Formatting
        title = metric.replace('_', ' ').replace('per kb', '(per kb)').title()
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Number of Genes')
        ax.grid(True, alpha=0.3)

        # Add reference lines (on transformed scale if applicable)
        ax.axvline(plot_mean, color='red', linestyle='--', alpha=0.8, linewidth=1.5, label='Mean')
        ax.axvline(plot_median, color='orange', linestyle=':', alpha=0.8, linewidth=1.5, label='Median')

        # Add legend for first plot
        if idx == 0:
            ax.legend()

    # Hide unused subplots
    for idx in range(len(available_metrics), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('Key Insertion Density Metrics Distribution',
                y=0.98)

    # Save to PDF
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Insertion density analysis script")
    parser.add_argument("-i", "--insertion_data_path", type=Path, required=True, help="Input TSV file with raw insertion read counts (4-level index, 2-level (Sample, Timepoint) header)")
    parser.add_argument("-a", "--annotations_path", type=Path, required=True, help="Input TSV file with annotations")
    parser.add_argument("-o", "--output_path", type=Path, required=True, help="Output TSV file with density statistics")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint column name")
    parser.add_argument("-f", "--final_timepoint", type=str, required=True, help="Final timepoint column name")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point of the script."""
    args = parse_args()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logger(log_level)

    # Validate input and output paths using the config dataclass
    start_time = time.time()

    try:
        config = InputOutputConfig(
            insertion_data_path=args.insertion_data_path,
            annotations_path=args.annotations_path,
            output_path=args.output_path,
            initial_timepoint=args.initial_timepoint,
            final_timepoint=args.final_timepoint,
        )

        # Get file sizes for result tracking
        insertion_file_size = config.insertion_data_path.stat().st_size
        annotation_file_size = config.annotations_path.stat().st_size

        logger.info("Starting insertion density analysis")
        logger.info(f"Insertion data file: {config.insertion_data_path} ({insertion_file_size:,} bytes)")
        logger.info(f"Annotations file: {config.annotations_path} ({annotation_file_size:,} bytes)")
        logger.info(f"Initial timepoint: {config.initial_timepoint}")
        logger.info(f"Final timepoint: {config.final_timepoint}")

        # Load data
        insertion_data = load_insertion_data(
            config.insertion_data_path, config.initial_timepoint, config.final_timepoint
        )
        annotations = load_annotation_data(config.annotations_path)

        # Filter for in-gene insertions
        in_gene_insertions = filter_in_gene_insertions(insertion_data, annotations)

        # valid genes
        valid_genes = in_gene_insertions["Systematic ID"].unique().tolist()
        logger.info(f"Analyzing {len(valid_genes)} genes")

        # Analyze each gene
        gene_results = []

        for gene_id in valid_genes:
            gene_insertions = in_gene_insertions[
                in_gene_insertions['Systematic ID'] == gene_id
            ]

            if len(gene_insertions) > 0:
                gene_analysis = analyze_gene_insertions(gene_id, gene_insertions)
                gene_results.append(gene_analysis)

        # Create results DataFrame
        results_df = pd.DataFrame(gene_results)
        results_df = results_df.set_index('Systematic ID')

        # Generate summary statistics
        stats = generate_summary_statistics(results_df)

        # Update result object with final statistics
        end_time = time.time()
        analysis_duration = end_time - start_time

        res = AnalysisResult(
            total_genes_analyzed=stats['total_genes_analyzed'],
            total_insertions_analyzed=stats['total_insertions_analyzed'],
            mean_insertion_density_per_kb_initial=stats['mean_insertion_density_per_kb_initial'],
            mean_insertion_density_per_kb_final=stats['mean_insertion_density_per_kb_final'],
            mean_insertion_density_log2fc=stats['mean_insertion_density_log2fc'],
            mean_gini_coefficient_of_depth_initial=stats['mean_gini_coefficient_of_depth_initial'],
            mean_strand_bias=stats['mean_strand_bias'],
        )

        # Generate histogram plots in PDF format
        plot_numeric_distributions_to_pdf(
            results_df, config.output_path, config.initial_timepoint, config.final_timepoint
        )

        # Save results
        results_df.to_csv(config.output_path, index=True, sep="\t")

        # Final summary
        logger.success("Analysis completed successfully")
        logger.success(f"Results saved to {config.output_path}")
        logger.success(f"Histogram PDF saved to {config.output_path.parent / f'{config.output_path.stem}_histograms.pdf'}")
        logger.success(f"Analyzed {len(results_df)} genes with insertion data")

        # Log summary statistics
        summary = asdict(res)
        logger.info(f"Analysis summary: {summary}")
        logger.info(f"Performance: {analysis_duration:.2f} seconds for {res.total_genes_analyzed} genes")

    except ValueError as e:
        logger.error(f"Validation error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

