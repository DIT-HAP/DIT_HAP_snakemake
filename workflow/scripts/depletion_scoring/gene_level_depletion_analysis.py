"""
Gene-Level Depletion Analysis for Transposon Insertion Sequencing

This script performs gene-level aggregation of transposon insertion depletion data,
combining multiple insertion sites within genes using weighted averages and statistical
meta-analysis methods. It processes log2 fold changes (LFC) and adjusted p-values
from insertion-level analysis to generate gene-level fitness scores.

Typical usage:
    python gene_level_depletion_analysis.py -l lfc_results.csv -a annotations.csv -o gene_results.csv

Input: CSV files with insertion-level LFC/p-values and genomic annotations
Output: CSV file with gene-level fitness scores and combined p-values
"""

import numpy as np
import pandas as pd
import sys
import time
from pathlib import Path
from typing import Tuple, Dict, Union
import argparse

from loguru import logger
from pydantic import BaseModel, field_validator, Field

# ======================== Configuration & Models ========================

class AnalysisConfig(BaseModel):
    lfc_path: Path = Field(..., description="Path to the LFC file")
    weights_path: Path = Field(..., description="Path to the weights file")
    annotations_path: Path = Field(..., description="Path to the annotations file")

    @field_validator("lfc_path", "weights_path", "annotations_path")
    def validate_input_exists(cls, v, field):
        if not v.exists():
            raise ValueError(f"Input file {v} does not exist")
        return v

    class Config:
        frozen = True

# ======================== Logging Setup ========================

def setup_logging() -> None:
    """Configure loguru for gene-level depletion analysis."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level="INFO",
        colorize=False
    )

# ======================== Data Loading ========================

@logger.catch
def load_data(config: AnalysisConfig) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Index]:
    """Load and validate all required data files."""
    logger.info(f"Loading LFC data from {config.lfc_path}")
    logger.info(f"Loading weights from {config.weights_path}")
    logger.info(f"Loading annotations from {config.annotations_path}")
    
    try:
        lfc_df = pd.read_csv(config.lfc_path, index_col=[0, 1, 2, 3], sep="\t")
        weights_df = pd.read_csv(config.weights_path, index_col=[0, 1, 2, 3], sep="\t")
        annotations_df = pd.read_csv(config.annotations_path, index_col=[0, 1, 2, 3], sep="\t")
        
        in_gene_insertions = annotations_df.query(
            "(Type != 'Intergenic region') and (Distance_to_stop_codon > 4)"
        ).index

        # transform weights
        transformed_weights_df = -np.log10(weights_df.fillna(1).clip(lower=1e-6, upper=1-1e-6))
        
        logger.info(f"Loaded {lfc_df.shape[0]} total insertions")
        logger.info(f"Found {len(in_gene_insertions)} in-gene insertions")
        
        return lfc_df, transformed_weights_df, annotations_df, in_gene_insertions
        
    except Exception as e:
        raise ValueError(f"Error loading data: {e}")

# ======================== Data Processing ========================

@logger.catch
def filter_in_gene_data(lfc_df: pd.DataFrame, transformed_weights_df: pd.DataFrame, 
                       in_gene_insertions: pd.Index) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Filter LFC and weights data for in-gene insertions."""
    in_gene_lfc = lfc_df[lfc_df.index.isin(in_gene_insertions)].copy()
    in_gene_weights = transformed_weights_df[transformed_weights_df.index.isin(in_gene_insertions)].copy()
    
    in_gene_lfc = in_gene_lfc.rename_axis("Timepoint", axis=1)
    in_gene_weights = in_gene_weights.rename_axis("Timepoint", axis=1)
    
    # Handle NaN values
    lfc_nan_count = in_gene_lfc.isna().any(axis=1).sum()
    weight_nan_count = in_gene_weights.isna().any(axis=1).sum()
    
    if lfc_nan_count > 0:
        logger.warning(f"Found {lfc_nan_count} rows with NaN LFC values")
    if weight_nan_count > 0:
        logger.warning(f"Found {weight_nan_count} rows with NaN weights - filling with 1")
        in_gene_weights = in_gene_weights.fillna(1)
    
    logger.info(f"Processed {in_gene_lfc.shape[0]} in-gene measurements")
    return in_gene_lfc, in_gene_weights

@logger.catch
def prepare_weighted_data(lfc_df: pd.DataFrame, transformed_weights_df: pd.DataFrame) -> pd.DataFrame:
    """Merge LFC and weights data into a single DataFrame."""
    lfc_series = lfc_df.stack().to_frame("LFC")
    weights_series = transformed_weights_df.stack().to_frame("Weights")

    merged = pd.merge(lfc_series, weights_series, left_index=True, right_index=True)
    
    weight_stats = merged["Weights"].describe()
    logger.info(f"Weight stats: mean={weight_stats['mean']:.2f}, max={weight_stats['max']:.2f}")
    
    return merged

@logger.catch
def annotate_and_normalize(data: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    """Annotate data with gene information and normalize weights."""
    annotated = pd.merge(data, annotations, left_index=True, right_index=True, how="left")
    
    # Remove missing weights
    initial_count = len(annotated)
    annotated = annotated[annotated["Weights"].notna()].copy()
    if len(annotated) != initial_count:
        logger.warning(f"Removed {initial_count - len(annotated)} rows with missing weights")
    
    # Normalize weights within gene-timepoint groups
    annotated["Normalized_weights"] = annotated.groupby(
        ["Systematic ID", "Timepoint"]
    )["Weights"].transform(lambda x: x / x.sum())
    
    logger.info(f"Annotated {len(annotated)} insertions")
    return annotated

# ======================== Gene-Level Analysis ========================

@logger.catch
def calculate_gene_lfc(gene_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate gene-level LFC for a single gene across timepoints."""
    gene_data = gene_data.reset_index().set_index(
        ["Chr", "Coordinate", "Strand", "Target", "Timepoint"]
    )
    
    results = pd.DataFrame()
    
    for timepoint, tp_data in gene_data.groupby("Timepoint"):
        lfcs = tp_data["LFC"].values
        weights = tp_data["Normalized_weights"].values
        
        gene_lfc = np.average(lfcs, weights=weights)
        results.loc[timepoint, "LFC"] = gene_lfc
    
    return results.sort_index().round(3)

@logger.catch
def analyze_all_genes(annotated_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate gene-level statistics for all genes."""
    logger.info("Calculating gene-level statistics")
    
    gene_groups = annotated_data.groupby([
        "Systematic ID", "Name", "FYPOviability", "DeletionLibrary_essentiality"
    ])
    
    results = []
    total_genes = len(gene_groups)
    start_time = time.time()
    
    for idx, ((sys_id, name, viability, essentiality), gene_data) in enumerate(gene_groups, 1):
        try:
            gene_results = calculate_gene_lfc(gene_data)
            
            for timepoint, lfc in gene_results["LFC"].items():
                results.append({
                    "Systematic ID": sys_id,
                    "Name": name,
                    "FYPOviability": viability,
                    "DeletionLibrary_essentiality": essentiality,
                    "Timepoint": timepoint,
                    "LFC": lfc
                })
                
        except Exception as e:
            logger.error(f"Error processing gene {sys_id}: {e}")
            continue
        
        if idx % 100 == 0 or idx == total_genes:
            elapsed = time.time() - start_time
            rate = idx / elapsed
            eta = (total_genes - idx) / rate if rate > 0 else 0
            logger.info(f"Processed {idx}/{total_genes} genes ({idx/total_genes*100:.1f}%) - ETA: {eta:.0f}s")
    
    # Pivot to wide format
    gene_df = pd.DataFrame(results)
    gene_wide = gene_df.pivot_table(
        index=["Systematic ID", "Name", "FYPOviability", "DeletionLibrary_essentiality"],
        columns="Timepoint",
        values="LFC",
        aggfunc="first"
    )
    
    # Remove genes with all NaN values
    gene_wide = gene_wide.dropna(how="all")
    
    logger.info(f"Completed analysis for {len(gene_wide)} genes")
    return gene_wide.reset_index()

# ======================== Summary Statistics ========================

def generate_summary(gene_df: pd.DataFrame) -> Dict[str, int]:
    """Generate summary statistics for the analysis."""
    return {
        'Total genes analyzed': len(gene_df),
        'FYPOviability: Essential genes': len(gene_df[gene_df['FYPOviability'] == 'inviable']),
        'FYPOviability: Non-essential genes': len(gene_df[gene_df['FYPOviability'] == 'viable']),
        'DeletionLibrary_essentiality: Essential genes': len(gene_df[gene_df['DeletionLibrary_essentiality'] == 'E']),
        'DeletionLibrary_essentiality: Non-essential genes': len(gene_df[gene_df['DeletionLibrary_essentiality'] == 'V'])
    }

def display_summary(stats: Dict[str, int]) -> None:
    """Display summary statistics."""
    logger.info("\n" + "="*60)
    logger.info("GENE-LEVEL DEPLETION ANALYSIS SUMMARY")
    logger.info("="*60)
    
    for key, value in stats.items():
        logger.info(f"{key:<40}: {value}")
    
    logger.info("="*60)

# ======================== Main ========================

@logger.catch
def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Gene-level depletion analysis for transposon insertion sequencing"
    )
    
    parser.add_argument('-l', '--lfc_path', type=Path, required=True,
                       help='Path to CSV file with LFC results')
    parser.add_argument('-a', '--annotations_path', type=Path, required=True,
                       help='Path to CSV file with annotations')
    parser.add_argument('-w', '--weights_path', type=Path, required=True,
                       help='Path to CSV file with weights')
    parser.add_argument('-o', '--output_path', type=Path, required=True,
                       help='Path for output CSV file')
    
    return parser.parse_args()

def main() -> None:
    """Execute the gene-level depletion analysis."""
    start_time = time.time()
    
    args = parse_arguments()
    setup_logging()
    
    logger.info("Starting gene-level depletion analysis")
    
    try:
        config = AnalysisConfig(
            lfc_path=args.lfc_path,
            weights_path=args.weights_path,
            annotations_path=args.annotations_path
        )
        
        args.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Load and process data
        lfc_df, transformed_weights_df, annotations_df, in_gene_insertions = load_data(config)
        in_gene_lfc, in_gene_weights = filter_in_gene_data(lfc_df, transformed_weights_df, in_gene_insertions)
        
        # Prepare weighted data and annotations
        weighted_data = prepare_weighted_data(in_gene_lfc, in_gene_weights)
        annotated_data = annotate_and_normalize(weighted_data, annotations_df)
        
        # Calculate gene-level statistics
        gene_results = analyze_all_genes(annotated_data)
        
        # Generate and display summary
        summary = generate_summary(gene_results)
        display_summary(summary)
        
        # Save results
        gene_results = gene_results.set_index("Systematic ID")
        gene_results.to_csv(args.output_path.parent / "LFC.tsv", sep="\t")
        transformed_weights_df.to_csv(args.output_path.parent / "transformed_weights.tsv", sep="\t")
        gene_results.to_csv(args.output_path, sep="\t")
        
        elapsed = time.time() - start_time
        logger.info(f"Analysis completed in {elapsed:.1f}s")
        logger.info(f"Results saved to: {args.output_path}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise

if __name__ == "__main__":
    main()
