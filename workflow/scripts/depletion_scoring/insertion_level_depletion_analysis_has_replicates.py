"""
Insertion-level depletion analysis script with replicates.

This script performs differential expression analysis on transposon insertion counts
using DESeq2 to identify genes that show depletion across time points. It handles
replicated samples and generates comprehensive statistical outputs including
log2 fold changes, p-values, and normalized counts.

The script processes count data with complex multi-level indexing and performs:
- Data preprocessing and quality control
- DESeq2 normalization and differential analysis
- Multiple testing correction
- Result visualization (MA plots)
- Comprehensive output files for downstream analysis

Typical Usage:
    python insertion_level_depletion_analysis_has_replicates.py \
        -i counts.tsv \
        -c control_insertions.tsv \
        -t 0h \
        -o output_directory/log2FoldChange.tsv

Input:
    - counts.tsv: Multi-index TSV file with insertion counts
    - control_insertions.tsv: TSV file with control insertion annotations

Output:
    - log2FoldChange.tsv: Primary output with log2 fold changes
    - Additional TSV files: baseMean, lfcSE, stat, pvalue, padj, normed_counts
    - PNG files: dispersions.png, MA_*.png plots
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import time
from loguru import logger
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# =============================== Configuration & Models ===============================

class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    counts_file: Path = Field(..., description="Path to the counts file")
    control_insertions_file: Path = Field(..., description="Path to the control insertions file")
    initial_timepoint: str = Field(..., description="Initial timepoint")
    output_file: Path = Field(..., description="Path to the output file")
    verbose: bool = Field(False, description="Enable verbose logging")

    @field_validator('counts_file', 'control_insertions_file')
    def validate_input_files(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True

class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    total_insertions_analyzed: int = Field(..., ge=0, description="Total number of insertions analyzed")
    timepoints_processed: int = Field(..., ge=0, description="Number of timepoints processed")
    control_insertions_count: int = Field(..., ge=0, description="Number of control insertions used")
    execution_time: float = Field(..., ge=0.0, description="Total execution time in seconds")

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
def load_and_preprocess_data(counts_file: Path, control_insertions_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame, List, List, pd.Index]:
    """Load and preprocess count data and control insertions."""
    logger.info(f"Loading counts data from {counts_file}")
    
    # Load counts data
    counts_df = pd.read_csv(counts_file, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t")
    counts_df_index_names = counts_df.index.names
    counts_df_columns_names = counts_df.columns.names

    counts_df.columns = ["#".join(col) for col in counts_df.columns]
    counts_df.index = ["=".join(map(str, index)) for index in counts_df.index]
    counts_df = counts_df.astype(int).T

    # Create metadata
    metadata = pd.DataFrame()
    metadata["sample"] = counts_df.index
    metadata["condition"] = [idx.split("#")[1] for idx in counts_df.index]
    metadata["group"] = [idx.split("#")[0] for idx in counts_df.index]
    metadata.set_index("sample", inplace=True)

    # Remove NA values
    counts_df = counts_df.loc[:, ~counts_df.isna().any(axis=0)].copy()

    # Load control insertions
    logger.info(f"Loading control insertions from {control_insertions_file}")
    control_insertion_annotations = pd.read_csv(control_insertions_file, index_col=[0, 1, 2, 3], sep="\t")
    control_insertion_annotations.index = ["=".join(map(str, index)) for index in control_insertion_annotations.index]

    logger.info(f"Loaded {len(counts_df.columns)} insertions and {len(control_insertion_annotations)} control insertions")
    
    return counts_df, metadata, counts_df_index_names, counts_df_columns_names, control_insertion_annotations.index

@logger.catch
def create_deseq_dataset(counts_df: pd.DataFrame, metadata: pd.DataFrame, control_insertions: pd.Index, initial_timepoint: str = "0h") -> DeseqDataSet:
    """Create and fit DESeq2 dataset for differential analysis."""
    logger.info("Creating DESeq2 dataset")
    
    inference = DefaultInference(n_cpus=36)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
        min_replicates=7,
    )
    
    logger.info("Fitting size factors using control insertions")
    dds.fit_size_factors(control_genes=control_insertions)
    logger.info("Fitting genewise dispersions")
    dds.fit_genewise_dispersions()
    logger.info("Fitting dispersion trend")
    dds.fit_dispersion_trend()
    logger.info("Fitting dispersion prior")
    dds.fit_dispersion_prior()
    logger.info("Fitting MAP dispersions")
    dds.fit_MAP_dispersions()
    logger.info("Fitting LFC")
    dds.fit_LFC()
    logger.info("Calculating Cook's distances")
    dds.calculate_cooks()
    if dds.refit_cooks:
        logger.info("Refitting after outlier removal")
        dds.refit()
    
    return dds

@logger.catch
def perform_differential_analysis(dds: DeseqDataSet, timepoints: List[str], initial_timepoint: str = "0h") -> Dict[str, DeseqStats]:
    """Perform differential expression analysis for all timepoints."""
    logger.info(f"Performing differential analysis for {len(timepoints)} timepoints")
    
    stat_res = {}
    inference = DefaultInference(n_cpus=36)
    
    for tp in timepoints:
        logger.info(f"Analyzing timepoint: {tp} vs {initial_timepoint}")
        stat_res[tp] = DeseqStats(
            dds, contrast=["condition", tp, initial_timepoint], inference=inference, 
            cooks_filter=True, independent_filter=True, quiet=True
        )
        stat_res[tp].summary()
        # Uncomment the following line if you want to perform LFC shrinkage
        # stat_res[tp].lfc_shrink(coeff=f"condition[T.{tp}]")
    
    return stat_res

@logger.catch
def plot_ma(stat_res: Dict[str, DeseqStats], output_dir: Path) -> None:
    """Generate MA plots for all timepoint comparisons."""
    logger.info("Generating MA plots")
    
    for tp, res in stat_res.items():
        output_path = output_dir / f"MA_{tp}.png"
        res.plot_MA(save_path=output_path)
        logger.debug(f"Saved MA plot for {tp} to {output_path}")

@logger.catch
def concatenate_results(stat_res: Dict[str, DeseqStats], timepoints: List[str]) -> pd.DataFrame:
    """Concatenate results from all timepoints into a single DataFrame."""
    logger.info("Concatenating results from all timepoints")
    
    result_df = {}
    for tp in timepoints:
        result_df[tp] = stat_res[tp].results_df
        result_df[tp]["log2FoldChange"] = -result_df[tp]["log2FoldChange"]
        result_df[tp]["stat"] = -result_df[tp]["stat"]
    concated_results = pd.concat(result_df, axis=1)
    concated_results.index = pd.MultiIndex.from_tuples(
        concated_results.index.str.split("=").tolist())
    # Convert string format numbers to numeric values in the MultiIndex
    new_index = []
    for idx in concated_results.index:
        chr_name = idx[0]
        # Convert coordinate from string to integer
        coordinate = int(idx[1]) if idx[1].isdigit() else idx[1]
        strand = idx[2]
        target = idx[3]
        new_index.append((chr_name, coordinate, strand, target))
    
    # Create a new MultiIndex with the converted values
    concated_results.index = pd.MultiIndex.from_tuples(
        new_index, names=concated_results.index.names)
    return concated_results

@logger.catch
def transform_index_to_multiindex(dds: DeseqDataSet, layer_name: str) -> pd.DataFrame:
    """Transform DESeq2 layer data to multi-index DataFrame."""
    logger.debug(f"Transforming {layer_name} layer to multi-index")
    
    df = pd.DataFrame(dds.layers[layer_name], index=dds.obs.index.tolist(), columns=dds.var.index.tolist()).T
    df.index = pd.MultiIndex.from_tuples(df.index.str.split("=").tolist())
    new_index = []
    for idx in df.index:
        chr_name = idx[0]
        # Convert coordinate from string to integer
        coordinate = int(idx[1]) if idx[1].isdigit() else idx[1]
        strand = idx[2]
        target = idx[3]
        new_index.append((chr_name, coordinate, strand, target))
    
    # Create a new MultiIndex with the converted values
    df.index = pd.MultiIndex.from_tuples(
        new_index, names=df.index.names)
    df.columns = pd.MultiIndex.from_tuples(df.columns.str.split("#").tolist())

    return df

# =============================== Main Function ===============================

def parse_arguments():
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perform differential expression analysis on insertion counts.")
    parser.add_argument("-i", "--counts_file", type=Path, required=True, help="Path to the counts file.")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint to analyze.")
    parser.add_argument("-c", "--control_insertions_file", type=Path, required=True, help="Path to the control insertions file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point for insertion-level depletion analysis."""
    start_time = time.time()
    
    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)
    
    try:
        # Validate input and output paths using the Pydantic model
        config = InputOutputConfig(
            counts_file=args.counts_file,
            control_insertions_file=args.control_insertions_file,
            initial_timepoint=args.initial_timepoint,
            output_file=args.output,
            verbose=args.verbose
        )
        
        logger.info("Starting insertion-level depletion analysis")
        logger.info(f"Counts file: {config.counts_file}")
        logger.info(f"Control insertions file: {config.control_insertions_file}")
        logger.info(f"Initial timepoint: {config.initial_timepoint}")
        logger.info(f"Output file: {config.output_file}")
        
        # Load and preprocess data
        counts_df, metadata, counts_df_index_names, counts_df_columns_names, control_insertions = load_and_preprocess_data(
            config.counts_file, config.control_insertions_file
        )

        logger.info("Metadata for analysis:")
        logger.info(f"\n{metadata}")
        logger.info(f"Number of control insertions: {len(control_insertions)}")

        timepoints = metadata["condition"].unique().tolist()
        timepoints.remove(config.initial_timepoint)
        logger.info(f"Control timepoint: {config.initial_timepoint}")
        logger.info(f"Timepoints for analysis: {timepoints}")
        
        # Create DESeq2 dataset
        dds = create_deseq_dataset(counts_df, metadata, control_insertions, config.initial_timepoint)
        logger.info("Plotting dispersions...")
        dds.plot_dispersions(save_path=config.output_file.parent / "dispersions.png")

        # Transform normalized counts
        logger.info("Transforming normalized counts to multi-index...")
        normalized_counts = transform_index_to_multiindex(dds, "normed_counts")
        normalized_counts = normalized_counts.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
        normalized_counts.to_csv(config.output_file.parent / "normed_counts.tsv", index=True, float_format="%.3f", sep="\t")

        # Transform count matrix
        logger.info("Transforming count matrix to multi-index...")
        count_X = pd.DataFrame(dds.X, index=dds.obs.index.tolist(), columns=dds.var.index.tolist()).T
        count_X.index = pd.MultiIndex.from_tuples(count_X.index.str.split("=").tolist())
        count_X.columns = pd.MultiIndex.from_tuples(count_X.columns.str.split("#").tolist())
        count_X = count_X.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
        count_X.to_csv(config.output_file.parent / "count_X.tsv", index=True, float_format="%.3f", sep="\t")

        # Transform Cook's distances
        cooks_df = transform_index_to_multiindex(dds, "cooks")
        cooks_df = cooks_df.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
        cooks_df.to_csv(config.output_file.parent / "cooks.tsv", index=True, float_format="%.3f", sep="\t")

        # Perform differential analysis
        logger.info("Performing differential analysis...")
        stat_res = perform_differential_analysis(dds, timepoints, config.initial_timepoint)
        plot_ma(stat_res, config.output_file.parent)

        # Concatenate results
        logger.info("Concatenating results...")
        concated_results = concatenate_results(stat_res, timepoints)
        concated_results = concated_results.rename_axis(counts_df_index_names, axis=0).rename_axis(["Timepoint", "Statistic"], axis=1)

        # Add the metrics for the initial timepoint
        logger.info("Adding metrics for the initial timepoint...")
        logger.info("All timepoints share the same baseMean")
        baseMean_initial = concated_results.xs("baseMean", axis=1, level="Statistic").iloc[:,0]
        logger.info("Set the initial log2FoldChange to 0")
        log2FoldChange_initial = 0
        logger.info("Set the initial lfcSE to NaN")
        lfcSE_initial = np.nan
        logger.info("Set the initial stat to NaN")
        stat_initial = np.nan
        logger.info("Set the initial pvalue to 1")
        pvalue_initial = 1
        logger.info("Set the initial padj to 1")
        padj_initial = 1
        logger.info("Insert the initial timepoint metrics to the first column...")
        concated_results.insert(0, (config.initial_timepoint,"padj"), padj_initial)
        concated_results.insert(0, (config.initial_timepoint,"pvalue"), pvalue_initial)
        concated_results.insert(0, (config.initial_timepoint,"stat"), stat_initial)
        concated_results.insert(0, (config.initial_timepoint,"lfcSE"), lfcSE_initial)
        concated_results.insert(0, (config.initial_timepoint,"log2FoldChange"), log2FoldChange_initial)
        concated_results.insert(0, (config.initial_timepoint,"baseMean"), baseMean_initial)
        
        numeric_columns = {"baseMean": 3, "log2FoldChange": 3, "lfcSE": 3, "stat": 3, "pvalue": 6, "padj": 6}
        logger.info("Rounding numeric columns...")
        for stat_name, decimal_places in numeric_columns.items():
            stat_columns = concated_results.xs(stat_name, axis=1, level="Statistic", drop_level=False)
            concated_results[stat_columns.columns] = stat_columns.round(decimal_places)
        logger.info("Saving results...")
        concated_results.to_csv(config.output_file.parent / "insertion_level_statistics.tsv", sep="\t")

        # Save individual statistic files
        baseMean_df = concated_results.xs("baseMean", axis=1, level="Statistic")
        baseMean_df.to_csv(config.output_file.parent/"baseMean.tsv", index=True, sep="\t")

        LFC_df = concated_results.xs("log2FoldChange", axis=1, level="Statistic")
        LFC_df.to_csv(config.output_file, index=True, sep="\t")

        lfcSE_df = concated_results.xs("lfcSE", axis=1, level="Statistic")
        lfcSE_df.to_csv(config.output_file.parent/"lfcSE.tsv", index=True, sep="\t")

        stat_df = concated_results.xs("stat", axis=1, level="Statistic")
        stat_df.to_csv(config.output_file.parent/"stat.tsv", index=True, sep="\t")

        pvalue_df = concated_results.xs("pvalue", axis=1, level="Statistic")
        pvalue_df.to_csv(config.output_file.parent/"pvalue.tsv", index=True, sep="\t")

        padj_df = concated_results.xs("padj", axis=1, level="Statistic")
        padj_df.to_csv(config.output_file.parent/"padj.tsv", index=True, sep="\t")
        
        # Create analysis result
        end_time = time.time()
        execution_time = end_time - start_time
        result = AnalysisResult(
            total_insertions_analyzed=len(counts_df.columns),
            timepoints_processed=len(timepoints),
            control_insertions_count=len(control_insertions),
            execution_time=execution_time
        )
        
        logger.success(f"Analysis completed in {execution_time:.2f} seconds")
        logger.info(f"Analyzed {result.total_insertions_analyzed} insertions across {result.timepoints_processed} timepoints")
        logger.info(f"Used {result.control_insertions_count} control insertions for normalization")
        logger.success(f"Results saved to {config.output_file.parent}")
        
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
