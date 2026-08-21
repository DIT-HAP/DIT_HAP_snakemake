"""Insertion-level depletion analysis for replicated samples (PyDESeq2).

Differential-abundance analysis of transposon insertion counts using PyDESeq2
(DESeq2). Identifies insertions that deplete across time points relative to an
initial time point, using replicated samples.

Counts are read with a four-level row MultiIndex (chromosome, coordinate,
strand, target) and a two-level column MultiIndex (group, condition). Size
factors are estimated from a supplied set of control insertions, dispersions
and log2 fold changes are fitted, Cook's-distance outliers are refit, and a
Wald test is run for each non-initial time point versus the initial time point.
Log2 fold changes and Wald statistics are negated so that depletion is reported
as a negative fold change.

Input
-----
- Counts TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex, ``header=[0, 1]`` column
  MultiIndex, tab-separated integer insertion counts.
- Control-insertions TSV: ``index_col=[0, 1, 2, 3]`` row MultiIndex, tab-separated;
  its index selects the control insertions used for size-factor estimation.

Output
------
- ``AnalysisResult`` dataclass, per-statistic DataFrames, and the dispersion
  data TSV writer, as used by the invoking script.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class AnalysisResult:
    """Summary statistics describing a completed analysis run."""
    total_insertions_analyzed: int
    timepoints_processed: int
    control_insertions_count: int
    execution_time: float


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def load_and_preprocess_data(counts_file: Path, control_insertions_file: Path) -> tuple[pd.DataFrame, pd.DataFrame, list, list, pd.Index]:
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
def perform_differential_analysis(dds: DeseqDataSet, timepoints: list[str], initial_timepoint: str = "0h") -> dict[str, DeseqStats]:
    """Perform differential expression analysis for all timepoints."""
    logger.info(f"Performing differential analysis for {len(timepoints)} timepoints")

    stat_res = {}
    inference = DefaultInference(n_cpus=36)

    for tp in timepoints:
        logger.info(f"Analyzing timepoint: {tp} vs {initial_timepoint}")
        stat_res[tp] = DeseqStats(
            dds, contrast=["condition", initial_timepoint, tp], inference=inference,
            cooks_filter=True, independent_filter=True, quiet=True
        )
        stat_res[tp].summary()
        # Uncomment the following line if you want to perform LFC shrinkage
        # stat_res[tp].lfc_shrink(coeff=f"condition[T.{tp}]")

    return stat_res


@logger.catch
def write_dispersion_data_tsv(dds: DeseqDataSet, output_path: Path) -> None:
    """Write per-insertion normed means and dispersion estimates for the rendering layer."""
    logger.info("Writing dispersion data")

    dispersion_df = pd.DataFrame(
        {
            "normed_mean": dds.var["_normed_means"].values,
            "genewise_dispersion": dds.var["genewise_dispersions"].values,
            "MAP_dispersion": dds.var["dispersions"].values,
            "fitted_dispersion": dds.var["fitted_dispersions"].values,
        },
        index=pd.MultiIndex.from_tuples(
            [tuple(idx.split("=")) for idx in dds.var.index],
            names=["Chr", "Coordinate", "Strand", "Target"],
        ),
    )
    dispersion_df.to_csv(output_path, sep="\t", index=True)


@logger.catch
def concatenate_results(stat_res: dict[str, DeseqStats], timepoints: list[str]) -> pd.DataFrame:
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
