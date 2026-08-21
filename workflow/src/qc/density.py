"""Pure computation for insertion density analysis.

Moved out of ``workflow/scripts/quality_control/insertion_density_analysis.py``
so the density metrics (per-gene insertion density, gap statistics, Gini
coefficients, read-depth statistics, strand-preference statistics) can be
imported and unit tested independently of the CLI entrypoint.
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
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
