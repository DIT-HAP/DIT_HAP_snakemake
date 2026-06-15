"""
Gene Coverage Analysis Script

This script analyzes transposon insertion data to determine gene coverage across
protein coding genes in Schizosaccharomyces pombe. It generates donut charts
showing coverage statistics for different gene categories.

The script processes two input files:
1. Gene_level_statistics_fitted.tsv - Contains genes with transposon insertions
2. gene_IDs_names_products.tsv - Contains all genes in the genome with metadata

Key features:
- Focuses on protein coding genes only
- Categorizes genes by product type (dubious, retrotransposable elements, etc.)
- Generates donut charts for overall coverage and filtered coverage
- Outputs PNG files with transparent backgrounds
- Provides comprehensive statistics

Usage:
    python gene_coverage_analysis.py -c gene_level_statistics_fitted.tsv -a gene_IDs_names_products.tsv -o output_dir
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
import time
from loguru import logger
from pydantic import BaseModel, Field, field_validator

try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from tabulate import tabulate
except ImportError as e:
    logger.error(f"Error importing required packages: {e}")
    logger.error("Please install required packages: pip install pandas numpy matplotlib tabulate")
    sys.exit(1)

# =============================== Constants ===============================
# Color palette for plotting (following cursor rules)
COLORS = {
    'covered': '#962955',      # Deep pink-purple
    'not_covered': '#7fb775',  # Medium green
    'dubious': '#6479cc',      # Medium blue
    'transposon': '#ad933c',   # Golden brown
    'sp_specific': '#26b1fd',  # Bright blue
    'other': '#8c397b'         # Medium purple
}

# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    covered_genes_file: Path = Field(..., description="Path to gene level statistics file (TSV format)")
    all_genes_file: Path = Field(..., description="Path to all genes metadata file (TSV format)")
    output_dir: Path = Field(..., description="Output directory for plots and statistics")
    verbose: bool = Field(False, description="Enable verbose logging")

    @field_validator('covered_genes_file')
    def validate_covered_genes_file(cls, v):
        if not v.exists():
            raise ValueError(f"Covered genes file does not exist: {v}")
        return v
    
    @field_validator('all_genes_file')
    def validate_all_genes_file(cls, v):
        if not v.exists():
            raise ValueError(f"All genes file does not exist: {v}")
        return v
    
    @field_validator('output_dir')
    def validate_output_dir(cls, v):
        v.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True

class CoverageStats(BaseModel):
    """Pydantic model to hold and validate coverage statistics."""
    total: int = Field(..., ge=0, description="Total number of genes")
    covered: int = Field(..., ge=0, description="Number of covered genes")
    not_covered: int = Field(..., ge=0, description="Number of not covered genes")
    coverage_pct: float = Field(..., ge=0.0, le=100.0, description="Coverage percentage")

class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    overall_stats: CoverageStats = Field(..., description="Overall coverage statistics")
    filtered_stats: CoverageStats = Field(..., description="Filtered coverage statistics")
    category_stats: Dict[str, CoverageStats] = Field(..., description="Category-specific coverage statistics")
    total_processing_time: float = Field(..., ge=0.0, description="Total processing time in seconds")

# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )

# =============================== Core Functions ===============================
@logger.catch
def load_covered_genes(file_path: Path) -> Set[str]:
    """
    Load genes with transposon insertions from statistics file.
    
    Args:
        file_path: Path to gene level statistics file
        
    Returns:
        Set of systematic gene IDs that have coverage
    """
    logger.info(f"Loading covered genes from {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        covered_genes = set(df['Systematic ID'].astype(str))
        logger.info(f"Loaded {len(covered_genes)} covered genes")
        return covered_genes
    except Exception as e:
        logger.error(f"Error loading covered genes: {e}")
        raise

@logger.catch
def load_all_genes(file_path: Path) -> pd.DataFrame:
    """
    Load all genes metadata from gene information file.
    
    Args:
        file_path: Path to all genes metadata file
        
    Returns:
        DataFrame with all gene information
    """
    logger.info(f"Loading all genes metadata from {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        logger.info(f"Loaded {len(df)} total genes")
        return df
    except Exception as e:
        logger.error(f"Error loading all genes: {e}")
        raise

@logger.catch
def filter_protein_coding_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter genes to include only protein coding genes.
    
    Args:
        df: DataFrame with all genes
        
    Returns:
        DataFrame with only protein coding genes
    """
    logger.info("Filtering for protein coding genes")
    
    protein_coding = df[df['gene_type'] == 'protein coding gene'].copy()
    logger.info(f"Found {len(protein_coding)} protein coding genes")
    
    return protein_coding

@logger.catch
def categorize_genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Categorize genes by their product types.
    
    Args:
        df: DataFrame with protein coding genes
        
    Returns:
        Dictionary mapping category names to sets of gene IDs
    """
    logger.info("Categorizing genes by product type")
    
    categories = {
        'dubious': set(),
        'retrotransposable_element': set(),
        'sp_pombe_specific': set(),
        'sp_specific': set(),
        'other': set()
    }
    
    for _, row in df.iterrows():
        gene_id = str(row['gene_systematic_id'])
        product = str(row['gene_product']).lower()
        
        if 'dubious' in product:
            categories['dubious'].add(gene_id)
        elif 'retrotransposable element' in product or 'transposon tf2-type' in product:
            categories['retrotransposable_element'].add(gene_id)
        elif 'schizosaccharomyces pombe specific' in product or 's. pombe specific' in product:
            categories['sp_pombe_specific'].add(gene_id)
        elif 'schizosaccharomyces specific' in product:
            categories['sp_specific'].add(gene_id)
        else:
            categories['other'].add(gene_id)
    
    for cat, genes in categories.items():
        logger.info(f"Category '{cat}': {len(genes)} genes")
    
    return categories

@logger.catch
def calculate_coverage_stats(
    protein_coding_genes: Set[str],
    covered_genes: Set[str],
    categories: Dict[str, Set[str]]
) -> Dict[str, CoverageStats]:
    """
    Calculate coverage statistics for different gene categories.
    
    Args:
        protein_coding_genes: Set of all protein coding gene IDs
        covered_genes: Set of genes with transposon insertions
        categories: Dictionary of gene categories
        
    Returns:
        Dictionary with coverage statistics for each category
    """
    logger.info("Calculating coverage statistics")
    
    stats = {}
    
    # Overall protein coding gene coverage
    total_protein_coding = len(protein_coding_genes)
    covered_protein_coding = len(protein_coding_genes & covered_genes)
    
    stats['overall'] = CoverageStats(
        total=total_protein_coding,
        covered=covered_protein_coding,
        not_covered=total_protein_coding - covered_protein_coding,
        coverage_pct=(covered_protein_coding / total_protein_coding) * 100
    )
    
    # Coverage excluding special categories
    excluded_categories = {'dubious', 'retrotransposable_element', 'sp_pombe_specific', 'sp_specific'}
    excluded_genes = set()
    for cat in excluded_categories:
        excluded_genes.update(categories[cat])
    
    filtered_protein_coding = protein_coding_genes - excluded_genes
    filtered_covered = filtered_protein_coding & covered_genes
    
    stats['filtered'] = CoverageStats(
        total=len(filtered_protein_coding),
        covered=len(filtered_covered),
        not_covered=len(filtered_protein_coding) - len(filtered_covered),
        coverage_pct=(len(filtered_covered) / len(filtered_protein_coding)) * 100 if filtered_protein_coding else 0
    )
    
    # Individual category coverage
    for cat_name, cat_genes in categories.items():
        cat_protein_coding = cat_genes & protein_coding_genes
        cat_covered = cat_protein_coding & covered_genes
        
        if cat_protein_coding:
            stats[cat_name] = CoverageStats(
                total=len(cat_protein_coding),
                covered=len(cat_covered),
                not_covered=len(cat_protein_coding) - len(cat_covered),
                coverage_pct=(len(cat_covered) / len(cat_protein_coding)) * 100
            )
    
    return stats

@logger.catch
def create_donut_chart(
    stats: CoverageStats,
    title: str,
    output_path: Path,
    colors: Optional[List[str]] = None
) -> None:
    """
    Create a donut chart showing coverage statistics.
    
    Args:
        stats: CoverageStats object with coverage data
        title: Chart title
        output_path: Path to save the chart
        colors: List of colors to use for the chart
    """
    logger.info(f"Creating donut chart: {title}")
    
    if colors is None:
        colors = [COLORS['covered'], COLORS['not_covered']]
    
    # Data for the donut chart
    labels = ['Covered', 'Not Covered']
    sizes = [stats.covered, stats.not_covered]
    
    # Create figure with transparent background
    fig, ax = plt.subplots(figsize=(8, 8), facecolor='none')
    ax.set_facecolor('none')
    
    # Create donut chart
    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        pctdistance=0.85,
        textprops={'fontsize': 12, 'weight': 'bold'}
    )
    
    # Create donut hole
    centre_circle = plt.Circle((0, 0), 0.70, fc='white', linewidth=0, alpha=0)
    ax.add_artist(centre_circle)
    
    # Add title
    ax.set_title(title, fontsize=16, weight='bold', pad=20)
    
    # Add coverage percentage in center
    coverage_pct = stats.coverage_pct
    ax.text(0, 0, f'{coverage_pct:.1f}%\nCoverage', 
            horizontalalignment='center', verticalalignment='center',
            fontsize=20, weight='bold')
    
    # Add sample size
    total = stats.total
    ax.text(0, -0.3, f'n = {total}', 
            horizontalalignment='center', verticalalignment='center',
            fontsize=12, style='italic')
    
    # Remove axes
    ax.axis('equal')
    
    # Save with transparent background
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none', transparent=True)
    plt.close()
    
    logger.info(f"Saved donut chart to {output_path}")

@logger.catch
def display_summary_table(stats: Dict[str, CoverageStats]) -> None:
    """
    Display summary statistics in a formatted table.
    
    Args:
        stats: Dictionary with coverage statistics for each category
    """
    logger.info("Generating summary statistics table")
    
    table_data = []
    for category, data in stats.items():
        if category in ['overall', 'filtered']:
            category_display = category.replace('_', ' ').title()
        else:
            category_display = category.replace('_', ' ').title()
        
        table_data.append([
            category_display,
            data.total,
            data.covered,
            data.not_covered,
            f"{data.coverage_pct:.1f}%"
        ])
    
    headers = ['Category', 'Total Genes', 'Covered', 'Not Covered', 'Coverage %']
    
    print("\n" + "="*80)
    print("GENE COVERAGE ANALYSIS SUMMARY")
    print("="*80)
    print(tabulate(table_data, headers=headers, tablefmt='grid'))
    print("="*80)

@logger.catch
def save_statistics(stats: Dict[str, CoverageStats], output_path: Path) -> None:
    """
    Save detailed statistics to a TSV file.
    
    Args:
        stats: Dictionary with coverage statistics for each category
        output_path: Path to save the statistics file
    """
    logger.info(f"Saving detailed statistics to {output_path}")
    
    with open(output_path, 'w') as f:
        f.write('Category\tTotal_Genes\tCovered_Genes\tNot_Covered\tCoverage_Percentage\n')
        for category, data in stats.items():
            f.write(f'{category}\t{data.total}\t{data.covered}\t{data.not_covered}\t{data.coverage_pct:.2f}\n')
    
    logger.success(f"Saved detailed statistics to {output_path}")

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Analyze gene coverage from transposon insertion data")
    parser.add_argument("-c", "--covered-genes", type=Path, required=True, 
                       help="Path to gene level statistics file (TSV format)")
    parser.add_argument("-a", "--all-genes", type=Path, required=True,
                       help="Path to all genes metadata file (TSV format)")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path('.'),
                       help="Output directory for plots and statistics (default: current directory)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script."""
    start_time = time.time()
    
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)
    
    # Validate input and output paths using the Pydantic model
    try:
        config = InputOutputConfig(
            covered_genes_file=args.covered_genes,
            all_genes_file=args.all_genes,
            output_dir=args.output_dir,
            verbose=args.verbose
        )
        
        logger.info("Starting gene coverage analysis")
        logger.info(f"Covered genes file: {config.covered_genes_file}")
        logger.info(f"All genes file: {config.all_genes_file}")
        logger.info(f"Output directory: {config.output_dir}")
        
        # Load data
        covered_genes = load_covered_genes(config.covered_genes_file)
        all_genes_df = load_all_genes(config.all_genes_file)
        
        # Filter for protein coding genes
        protein_coding_df = filter_protein_coding_genes(all_genes_df)
        protein_coding_genes = set(protein_coding_df['gene_systematic_id'].astype(str))
        
        # Categorize genes
        categories = categorize_genes(protein_coding_df)
        
        # Calculate coverage statistics
        stats = calculate_coverage_stats(protein_coding_genes, covered_genes, categories)
        
        # Create donut charts
        # Overall coverage
        create_donut_chart(
            stats['overall'],
            'Overall Protein Coding Gene Coverage',
            config.output_dir / 'overall_gene_coverage.png'
        )
        
        # Filtered coverage (excluding special categories)
        create_donut_chart(
            stats['filtered'],
            'Protein Coding Gene Coverage\n(Excluding Dubious, Transposons, and Species-Specific)',
            config.output_dir / 'filtered_gene_coverage.png'
        )
        
        # Save detailed statistics
        stats_file = config.output_dir / 'coverage_statistics.tsv'
        save_statistics(stats, stats_file)
        
        # Display summary
        display_summary_table(stats)
        
        # Processing time
        end_time = time.time()
        processing_time = end_time - start_time
        logger.success(f"Analysis completed in {processing_time:.2f} seconds")
        
        # Create analysis result object
        result = AnalysisResult(
            overall_stats=stats['overall'],
            filtered_stats=stats['filtered'],
            category_stats={k: v for k, v in stats.items() if k not in ['overall', 'filtered']},
            total_processing_time=processing_time
        )
        
        logger.success(f"Analysis complete. Results saved to {config.output_dir}")
        
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
