#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (Optional) PEP 723 inline script metadata for self-contained execution with `uv`.
# Remove or adjust if managing dependencies via a traditional virtual environment.
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy",
#     "pandas",
#     "loguru",
#     "pybedtools",
# ]
# ///

"""
Extract Genome Regions from PomBase GFF3 Annotation
===================================================

Parse a PomBase GFF3 annotation together with supporting metadata to produce
BED-format interval files for downstream insertion-site annotation. Primary
coding transcripts are identified by matching accumulated CDS length to known
peptide lengths; intergenic regions are derived as the complement of the primary
transcript spans; non-coding RNA intervals are extracted from the annotation.

The pipeline additionally annotates each interval with gene names, FYPO
viability, and deletion-library essentiality, records parental-region spans,
finds overlapping coding regions via BedTools self-intersection, and derives a
non-coding RNA set that does not overlap expanded coding parental regions.

Input
-----
- ``--gff``            : PomBase GFF3 annotation (tab-separated, ``#`` comments).
- ``--fai``            : FASTA index (.fai) providing chromosome sizes.
- ``--peptide_stats``  : PomBase PeptideStats TSV (``Systematic_ID``, ``Residues``).
- ``--gene_ids``       : gene_IDs_names_products TSV (systematic id / name / synonyms).
- ``--fypo``           : FYPOviability TSV (headerless: systematic id, viability).
- ``--hayles``         : Hayles 2013 viability XLSX.

Output
------
- ``--out_primary``           : coding gene primary transcripts BED.
- ``--out_intergenic``        : intergenic regions BED.
- ``--out_ncrna``             : non-coding RNA BED.
- ``--out_genome_intervals``  : genome intervals BED (primary + intergenic).
- ``--out_overlapped``        : overlapping coding regions BED.
- A derived ``non_coding_rna_without_overlap_with_coding_gene.bed`` written
  alongside ``--out_ncrna``.

Usage
-----
    python extract_genome_region.py \
        --gff Schizosaccharomyces_pombe_all_chromosomes.gff3 \
        --fai Schizosaccharomyces_pombe_all_chromosomes.fa.fai \
        --peptide_stats PeptideStats.tsv \
        --gene_ids gene_IDs_names_products.tsv \
        --fypo FYPOviability.tsv \
        --hayles Hayles_2013_OB_merged_categories.xlsx \
        --out_primary coding_gene_primary_transcripts.bed \
        --out_intergenic intergenic_regions.bed \
        --out_ncrna non_coding_rna.bed \
        --out_genome_intervals Genome_intervals.bed \
        --out_overlapped overlapped_region.bed

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger
import pybedtools
from pybedtools import BedTool

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
CHR_ORDER = ["chr_II_telomeric_gap", "I", "II", "III", "mating_type_region", "mitochondrial"]
TRANSCRIPT_FEATURE_TYPES = frozenset(["mRNA", "tRNA", "rRNA", "snoRNA", "snRNA", "lncRNA"])
NON_CODING_RNA_FEATURES = ["tRNA", "rRNA", "snoRNA", "snRNA", "lncRNA"]
CODING_PARENTAL_EXPAND_BP = 200

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class Config:
    """Configuration for GFF3 genome region extraction."""
    gff_file: Path
    fai_file: Path
    peptide_stats_file: Path
    gene_ids_file: Path
    fypo_file: Path
    hayles_file: Path
    out_primary: Path
    out_intergenic: Path
    out_ncrna: Path
    out_genome_intervals: Path
    out_overlapped: Path

    def __post_init__(self) -> None:
        """Validate inputs and create output directories."""
        for field_name in ["gff_file", "fai_file", "peptide_stats_file", "gene_ids_file", "fypo_file", "hayles_file"]:
            p = getattr(self, field_name)
            if not Path(p).exists():
                raise ValueError(f"Input file does not exist: {p}")
        for field_name in ["out_primary", "out_intergenic", "out_ncrna", "out_genome_intervals", "out_overlapped"]:
            Path(getattr(self, field_name)).parent.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def update_sysID(genes: list, gene_IDs_names_products: pd.DataFrame) -> list:
    """Resolve gene identifiers to current systematic IDs via name/synonym lookup."""
    coding = gene_IDs_names_products.query("gene_type == 'protein coding gene'")
    synonyms2ID = (
        coding.set_index("gene_systematic_id")["synonyms"]
        .str.split(",").explode().str.strip().dropna()
        .reset_index().set_index("synonyms")
    )
    names2ID = (
        coding.set_index("gene_name")["gene_systematic_id"]
        .drop_duplicates().reset_index().set_index("gene_name")
    )
    sysIDs_now = coding["gene_systematic_id"].unique().tolist()

    updated = []
    for gene in genes:
        if pd.isna(gene):
            updated.append(gene)
            logger.debug(f"{gene} is NA")
        elif gene in sysIDs_now:
            updated.append(gene)
        elif gene in names2ID.index:
            val = names2ID.loc[gene, "gene_systematic_id"]
            if isinstance(val, str):
                updated.append(val)
                logger.debug(f"{gene} -> {val}")
            else:
                updated.append(np.nan)
                logger.warning(f"{gene} has multiple updates: {val.tolist()}")
        elif gene in synonyms2ID.index:
            val = synonyms2ID.loc[gene, "gene_systematic_id"]
            if isinstance(val, str):
                updated.append(val)
                logger.debug(f"{gene} -> {val}")
            else:
                updated.append(np.nan)
                logger.warning(f"{gene} has multiple synonym updates: {val.tolist()}")
        else:
            updated.append(gene)
            logger.debug(f"{gene} not found in gene metadata")
    return updated


@logger.catch
def get_gff_transcript_id(
    row: pd.Series,
    id_pattern: re.Pattern = re.compile(r"ID=([^:;]+)"),
    parent_pattern: re.Pattern = re.compile(r"Parent=([^;]+)"),
    transcript_feature_types: frozenset = TRANSCRIPT_FEATURE_TYPES,
) -> object:
    """Extract transcript identifier from a GFF row's attributes."""
    attributes = str(row["Attribute"])
    feature_type = row["Feature"]
    if feature_type in transcript_feature_types:
        match = id_pattern.search(attributes)
        if match:
            return match.group(1)
    elif "Parent=" in attributes:
        match = parent_pattern.search(attributes)
        if match:
            return match.group(1)
    return np.nan


@logger.catch
def parse_gff_data(gff_file_path: Path) -> pd.DataFrame:
    """Read and parse a GFF3 file, extracting Systematic ID and Transcript columns."""
    column_names = ["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]
    gff_df = pd.read_csv(
        gff_file_path,
        sep="\t",
        comment="#",
        names=column_names,
        dtype={"Chr": str, "Start": pd.Int64Dtype(), "End": pd.Int64Dtype()},
    )
    extract_systematic_ID_pattern = re.compile(r"ID=(\S+?)(?:$|(?:|\.\d(?::\S+|));))")
    gff_df["Systematic ID"] = gff_df["Attribute"].str.extract(extract_systematic_ID_pattern, expand=False)
    gff_df["Transcript"] = gff_df.apply(get_gff_transcript_id, axis=1)
    gff_df["Chr"] = pd.Categorical(gff_df["Chr"], categories=CHR_ORDER, ordered=True)
    logger.info(f"GFF parsed: {gff_df.shape[0]:,} rows, {gff_df['Systematic ID'].nunique():,} unique Systematic IDs")
    return gff_df


@logger.catch
def calculate_accumulated_cds_bases(transcript_features_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate accumulated CDS bases for features within a single transcript."""
    if "Length" not in transcript_features_df.columns or not pd.api.types.is_numeric_dtype(transcript_features_df["Length"]):
        raise ValueError("DataFrame must contain a numeric 'Length' column.")
    sorted_df = transcript_features_df.sort_values(["Start"]).copy()
    strand = sorted_df["Strand"].iloc[0]
    if strand == "+":
        iteration_indices = sorted_df.index
    elif strand == "-":
        iteration_indices = sorted_df.index[::-1]
    else:
        raise ValueError(f"Unknown strand: {strand}")
    current_cds_accumulation = 0
    sorted_df["Accumulated_CDS_bases"] = np.nan
    for idx in iteration_indices:
        sorted_df.loc[idx, "Accumulated_CDS_bases"] = current_cds_accumulation
        if sorted_df.loc[idx, "Feature"] == "CDS":
            current_cds_accumulation += sorted_df.loc[idx, "Length"]
    return sorted_df


@logger.catch
def gff_features_to_bed(
    transcript_features_group_df: pd.DataFrame,
    gene_type_label: str,
    peptide_length_map: dict,
) -> pd.DataFrame | None:
    """Convert GFF features for a transcript to BED-like format."""
    bed_df = transcript_features_group_df.copy()
    bed_df["Start"] = bed_df["Start"] - 1
    bed_df["Length"] = bed_df["End"] - bed_df["Start"]
    bed_columns_std = ["Chr", "Start", "End", "Transcript", "Length", "Strand"]
    other_info_cols = ["Feature", "Systematic ID", "Type"]
    current_gene_id = bed_df["Systematic ID"].iloc[0]
    current_transcript_id = bed_df["Transcript"].iloc[0]

    if gene_type_label == "Coding gene":
        boundary_feature_type = "CDS"
        cds_segments = bed_df[bed_df["Feature"] == boundary_feature_type]
        if cds_segments.empty:
            logger.warning(f"No CDS for transcript {current_transcript_id} of gene {current_gene_id}. Skipping.")
            return None
        total_cds_length = cds_segments["Length"].sum()
        bed_df = calculate_accumulated_cds_bases(bed_df)
        expected_peptide_len = peptide_length_map.get(current_gene_id, -1)
        if expected_peptide_len != -1 and total_cds_length > 0 and (total_cds_length % 3 == 0) and \
                (int(total_cds_length / 3) - 1 == expected_peptide_len):
            bed_df["Primary_transcript_flag"] = "Yes"
        elif expected_peptide_len != -1 and total_cds_length > 0 and (total_cds_length % 3 != 0) and \
                (int(total_cds_length // 3) - 1 == expected_peptide_len):
            bed_df["Primary_transcript_flag"] = "Yes"
            logger.warning(f"Gene:{current_gene_id} Transcript:{current_transcript_id}: CDS length not divisible by 3")
        else:
            bed_df["Primary_transcript_flag"] = "No"
        other_info_cols.extend(["Primary_transcript_flag", "Accumulated_CDS_bases"])
    elif gene_type_label == "Non-coding gene":
        boundary_feature_type = bed_df["Feature"].iloc[0]
    else:
        boundary_feature_type = "exon"

    boundary_defining_features = bed_df[bed_df["Feature"] == boundary_feature_type]
    if boundary_defining_features.empty:
        logger.warning(f"No '{boundary_feature_type}' to define boundaries for {current_transcript_id}. Using all features.")
        min_coord_start = bed_df["Start"].min()
        max_coord_end = bed_df["End"].max()
    else:
        min_coord_start = boundary_defining_features["Start"].min()
        max_coord_end = boundary_defining_features["End"].max()

    filtered_bed_df = bed_df[(bed_df["Start"] >= min_coord_start) & (bed_df["End"] <= max_coord_end)].copy()
    if filtered_bed_df.empty:
        logger.warning(f"No features after boundary filtering for {current_transcript_id}. Skipping.")
        return None
    filtered_bed_df.insert(3, "Type", gene_type_label)
    final_columns = bed_columns_std + other_info_cols
    for col in final_columns:
        if col not in filtered_bed_df.columns:
            filtered_bed_df[col] = np.nan
    output_df = (
        filtered_bed_df[final_columns]
        .rename(columns={"Chr": "#Chr"})
        .sort_values(["#Chr"] + bed_columns_std[1:])
    )
    return output_df


@logger.catch
def select_primary_transcripts(all_coding_features_bed_df: pd.DataFrame) -> pd.DataFrame:
    """Select one primary transcript per coding gene, preferring .1 when multiple candidates exist."""
    candidate_primary = all_coding_features_bed_df[
        all_coding_features_bed_df["Primary_transcript_flag"] == "Yes"
    ][["Systematic ID", "Transcript"]].drop_duplicates()

    primary_counts = candidate_primary.groupby("Systematic ID")["Transcript"].count()
    genes_with_multiple = primary_counts[primary_counts > 1].index.tolist()

    final_ids = []
    for gene_id, group in candidate_primary.groupby("Systematic ID"):
        if gene_id in genes_with_multiple:
            dot_one = [tid for tid in group["Transcript"] if tid.endswith(".1")]
            final_ids.append(dot_one[0] if dot_one else group["Transcript"].iloc[0])
        else:
            final_ids.append(group["Transcript"].iloc[0])

    primary_df = all_coding_features_bed_df[
        all_coding_features_bed_df["Transcript"].isin(final_ids)
    ].copy()
    primary_df["#Chr"] = pd.Categorical(primary_df["#Chr"], categories=CHR_ORDER, ordered=True)
    primary_df = primary_df.sort_values(["#Chr", "Start", "End"])
    logger.info(f"Selected {primary_df['Transcript'].nunique():,} primary transcripts")
    return primary_df


@logger.catch
def build_intergenic_bed(primary_transcripts_bed_df: pd.DataFrame, fai_file_path: Path) -> pd.DataFrame:
    """Compute intergenic regions as the BedTools complement of primary transcript spans."""
    primary_bt = BedTool.from_dataframe(primary_transcripts_bed_df)
    intergenic_bt = primary_bt.complement(g=str(fai_file_path))
    intergenic_df = intergenic_bt.to_dataframe(disable_auto_names=True, header=None)
    col_names = ["#Chr", "Start", "End"] + [f"col_{i+4}" for i in range(intergenic_df.shape[1] - 3)]
    intergenic_df.columns = col_names[: intergenic_df.shape[1]]
    logger.info(f"Identified {len(intergenic_df):,} intergenic regions")
    return intergenic_df


@logger.catch
def annotate_intergenic_region_flanks(
    intergenic_row: pd.Series,
    primary_transcripts_bed_df: pd.DataFrame,
) -> pd.Series:
    """Annotate an intergenic region with its flanking transcript information."""
    chrom = intergenic_row["#Chr"]
    intergenic_start = intergenic_row["Start"]
    intergenic_end = intergenic_row["End"]

    left = primary_transcripts_bed_df[
        (primary_transcripts_bed_df["#Chr"] == chrom) &
        (primary_transcripts_bed_df["End"] == intergenic_start)
    ]
    right = primary_transcripts_bed_df[
        (primary_transcripts_bed_df["#Chr"] == chrom) &
        (primary_transcripts_bed_df["Start"] == intergenic_end)
    ]

    def flank_vals(df: pd.DataFrame, side: str) -> dict:
        """Extract transcript, systematic ID, and strand for one flanking side."""
        if not df.empty:
            row = df.iloc[0]
            return {f"{side}_Transcript": row["Transcript"],
                    f"{side}_Systematic ID": row["Systematic ID"],
                    f"{side}_Strand": row["Strand"]}
        return {f"{side}_Transcript": "Boundary",
                f"{side}_Systematic ID": "Boundary",
                f"{side}_Strand": "Boundary"}

    left_info = flank_vals(left, "Left")
    right_info = flank_vals(right, "Right")

    ls, rs = left_info["Left_Strand"], right_info["Right_Strand"]
    if ls == "+" and rs == "-":
        orientation = "Convergent"
    elif ls == "-" and rs == "+":
        orientation = "Divergent"
    elif ls == "+" and rs == "+":
        orientation = "Tandem_Plus"
    elif ls == "-" and rs == "-":
        orientation = "Tandem_Minus"
    else:
        orientation = "Boundary_Adjacent"

    return pd.Series({
        "Transcript": left_info["Left_Transcript"] + "|" + right_info["Right_Transcript"],
        "Systematic ID": left_info["Left_Systematic ID"] + "|" + right_info["Right_Systematic ID"],
        "Strand": left_info["Left_Strand"] + "|" + right_info["Right_Strand"],
        "Feature": orientation + "_Region",
    })


@logger.catch
def find_overlapping_regions(feature) -> pybedtools.Interval:
    """Build a pybedtools Interval representing overlapping gene regions."""
    chr_a, chr_b = feature[0], feature[6]
    start_a, start_b = int(feature[1]), int(feature[7])
    end_a, end_b = int(feature[2]), int(feature[8])
    transcript_a, transcript_b = feature[3], feature[9]
    sysID_a, sysID_b = feature[4], feature[10]
    strand_a, strand_b = feature[5], feature[11]

    chrom = chr_a
    start = max(start_a, start_b)
    end = min(end_a, end_b)
    transcript = transcript_a
    sysID = sysID_a
    strand = strand_a
    score = ""
    if transcript_a != transcript_b:
        transcript = transcript_a + "," + transcript_b
        sysID = sysID_a + "," + sysID_b
        strand = strand_a + "," + strand_b
        score = "Overlapping genes"
    return pybedtools.create_interval_from_list([chrom, str(start), str(end), transcript, strand, score, sysID])


@logger.catch
def build_genome_intervals(
    primary_transcripts_bed_df: pd.DataFrame,
    intergenic_regions_df: pd.DataFrame,
) -> pd.DataFrame:
    """Concatenate primary transcripts and intergenic regions into a genome-wide interval file."""
    return pd.concat([
        primary_transcripts_bed_df.drop(columns=["Primary_transcript_flag"]),
        intergenic_regions_df,
    ])


@logger.catch
def run_pipeline(config: Config) -> None:
    """Execute the full GFF3 genome region extraction pipeline."""
    # --- Parse GFF ---
    gff_df = parse_gff_data(config.gff_file)

    # --- Load auxiliary data ---
    peptide_stats_df = pd.read_csv(config.peptide_stats_file, sep="\t")
    gene_to_peptide_length_map = dict(zip(peptide_stats_df["Systematic_ID"], peptide_stats_df["Residues"]))
    logger.info(f"Loaded peptide statistics for {len(gene_to_peptide_length_map):,} proteins")

    gene_IDs_names_products = pd.read_csv(config.gene_ids_file, sep="\t")
    gene_IDs_names_products["gene_name"] = gene_IDs_names_products["gene_name"].fillna(
        gene_IDs_names_products["gene_systematic_id"]
    )
    ID2name = dict(zip(gene_IDs_names_products["gene_systematic_id"], gene_IDs_names_products["gene_name"]))

    FYPOviability_df = pd.read_csv(config.fypo_file, sep="\t", header=None, names=["Systematic ID", "FYPOviability"])
    FYPOviability = FYPOviability_df.set_index("Systematic ID")["FYPOviability"].to_dict()

    Hayles_viability_df = pd.read_excel(config.hayles_file)
    Hayles_viability_df["Updated_Systematic_ID"] = update_sysID(
        Hayles_viability_df["Systematic ID"].tolist(), gene_IDs_names_products
    )
    DeletionLibrary_essentiality = dict(
        zip(Hayles_viability_df["Updated_Systematic_ID"],
            Hayles_viability_df["Gene dispensability. This study"].str.strip())
    )

    # --- Identify gene categories ---
    coding_gene_ids = gff_df[gff_df["Feature"] == "mRNA"]["Systematic ID"].unique().tolist()
    logger.info(f"Coding genes: {len(coding_gene_ids):,}")

    # --- Process coding genes ---
    coding_features_df = gff_df[
        gff_df["Systematic ID"].isin(coding_gene_ids) &
        gff_df["Feature"].isin(["CDS", "intron"])
    ].copy()
    logger.info(f"Processing {coding_features_df['Transcript'].nunique():,} coding transcript IDs")

    processed_beds = []
    for name, group in coding_features_df.groupby(["Systematic ID", "Transcript"]):
        result = gff_features_to_bed(group, "Coding gene", gene_to_peptide_length_map)
        if result is not None and not result.empty:
            processed_beds.append(result)

    if not processed_beds:
        logger.error("No coding gene BED features produced")
        return
    all_coding_features_bed_df = pd.concat(processed_beds).reset_index(drop=True)
    logger.info(f"Coding features BED: {all_coding_features_bed_df.shape[0]:,} rows")

    # --- Select primary transcripts ---
    primary_transcripts_bed_df = select_primary_transcripts(all_coding_features_bed_df)

    # --- Build intergenic regions ---
    intergenic_regions_df = build_intergenic_bed(primary_transcripts_bed_df, config.fai_file)
    intergenic_regions_df[["Transcript", "Systematic ID", "Strand", "Feature"]] = intergenic_regions_df.apply(
        annotate_intergenic_region_flanks,
        primary_transcripts_bed_df=primary_transcripts_bed_df,
        axis=1,
    )
    intergenic_regions_df["Length"] = intergenic_regions_df["End"] - intergenic_regions_df["Start"]
    intergenic_regions_df["Type"] = "Intergenic region"
    intergenic_regions_df = intergenic_regions_df[
        ["#Chr", "Start", "End", "Transcript", "Length", "Strand", "Feature", "Systematic ID", "Type"]
    ].copy()

    # --- Process non-coding RNAs ---
    non_coding_rna_df = gff_df[gff_df["Feature"].isin(NON_CODING_RNA_FEATURES)].copy().sort_values(
        ["Feature", "Chr", "Start", "End", "Systematic ID", "Transcript"]
    )
    non_coding_rna_bed_df = (
        non_coding_rna_df
        .groupby("Feature")
        .apply(gff_features_to_bed, "Non-coding gene", gene_to_peptide_length_map)
        .reset_index(drop=True)
    )
    logger.info(f"Non-coding RNA BED: {len(non_coding_rna_bed_df):,} rows")

    # --- Add gene name annotation ---
    for bed in [primary_transcripts_bed_df, intergenic_regions_df, non_coding_rna_bed_df]:
        if not bed.empty:
            if bed["Systematic ID"].astype(str).str.contains("|", regex=False).any():
                bed["Name"] = bed["Systematic ID"].apply(
                    lambda x: "|".join([ID2name.get(i, i) for i in str(x).split("|")])
                )
            else:
                bed["Name"] = bed["Systematic ID"].map(ID2name)

    # --- Add essentiality annotation ---
    for bed in [primary_transcripts_bed_df, intergenic_regions_df, non_coding_rna_bed_df]:
        if not bed.empty:
            if bed["Systematic ID"].astype(str).str.contains("|", regex=False).any():
                bed["FYPOviability"] = bed["Systematic ID"].apply(
                    lambda x: "|".join([FYPOviability.get(i, i) for i in str(x).split("|")])
                )
                bed["DeletionLibrary_essentiality"] = bed["Systematic ID"].apply(
                    lambda x: "|".join([DeletionLibrary_essentiality.get(i, "Not_determined") for i in str(x).split("|")])
                )
            else:
                bed["FYPOviability"] = bed["Systematic ID"].map(FYPOviability)
                bed["DeletionLibrary_essentiality"] = bed["Systematic ID"].map(DeletionLibrary_essentiality)

    # --- Add parental region info ---
    for bed in [primary_transcripts_bed_df, intergenic_regions_df, non_coding_rna_bed_df]:
        if not bed.empty:
            bed["ParentalRegion_start"] = bed.groupby("Systematic ID")["Start"].transform("min")
            bed["ParentalRegion_end"] = bed.groupby("Systematic ID")["End"].transform("max")
            bed["ParentalRegion_length"] = bed["ParentalRegion_end"] - bed["ParentalRegion_start"]

    # --- Find overlapping coding regions ---
    primary_spans = primary_transcripts_bed_df[
        ["#Chr", "ParentalRegion_start", "ParentalRegion_end", "Transcript", "Systematic ID", "Strand"]
    ].drop_duplicates()
    primary_spans_bt = BedTool.from_dataframe(primary_spans)
    overlaps = primary_spans_bt.intersect(primary_spans_bt, wa=True, wb=True, header=True)
    results_df = (
        overlaps.each(find_overlapping_regions)
        .saveas()
        .to_dataframe(names=["#Chr", "Start", "End", "Transcript", "Strand", "Feature", "Systematic ID"])
        .drop_duplicates(subset=["#Chr", "Start", "End"], keep="first")
    )
    overlapped_region_bed_df = results_df[results_df["Feature"] == "Overlapping genes"].copy()
    overlapped_region_bed_df["Length"] = overlapped_region_bed_df["End"] - overlapped_region_bed_df["Start"]
    overlapped_region_bed_df["Type"] = "Coding gene"
    overlapped_region_bed_df = overlapped_region_bed_df[
        ["#Chr", "Start", "End", "Transcript", "Length", "Strand", "Feature", "Systematic ID", "Type"]
    ].copy()
    overlapped_region_bed_df["Name"] = overlapped_region_bed_df["Systematic ID"].apply(
        lambda x: ",".join([ID2name.get(i, i) for i in str(x).split(",")])
    )
    overlapped_region_bed_df["FYPOviability"] = overlapped_region_bed_df["Systematic ID"].apply(
        lambda x: ",".join([FYPOviability.get(i, i) for i in str(x).split(",")])
    )
    overlapped_region_bed_df["DeletionLibrary_essentiality"] = overlapped_region_bed_df["Systematic ID"].apply(
        lambda x: ",".join([DeletionLibrary_essentiality.get(i, "Not_determined") for i in str(x).split(",")])
    )
    logger.info(f"Overlapping coding regions: {len(overlapped_region_bed_df):,}")

    # --- Save outputs ---
    primary_transcripts_bed_df.drop(columns=["Primary_transcript_flag"]).to_csv(
        config.out_primary, sep="\t", index=False
    )
    logger.success(f"Saved primary transcripts → {config.out_primary}")

    intergenic_regions_df.to_csv(config.out_intergenic, sep="\t", index=False)
    logger.success(f"Saved intergenic regions → {config.out_intergenic}")

    non_coding_rna_bed_df.to_csv(config.out_ncrna, sep="\t", index=False)
    logger.success(f"Saved non-coding RNAs → {config.out_ncrna}")

    build_genome_intervals(primary_transcripts_bed_df, intergenic_regions_df).to_csv(
        config.out_genome_intervals, sep="\t", index=False
    )
    logger.success(f"Saved genome intervals → {config.out_genome_intervals}")

    overlapped_region_bed_df.to_csv(config.out_overlapped, sep="\t", index=False)
    logger.success(f"Saved overlapped regions → {config.out_overlapped}")

    # --- Derive non-coding RNA without overlap with coding genes ---
    expanded = primary_transcripts_bed_df[
        ["#Chr", "ParentalRegion_start", "ParentalRegion_end", "Transcript", "ParentalRegion_length", "Strand", "Type", "Systematic ID"]
    ].copy()
    expanded["ParentalRegion_start"] = (expanded["ParentalRegion_start"] - CODING_PARENTAL_EXPAND_BP).clip(lower=0)
    expanded["ParentalRegion_end"] = expanded["ParentalRegion_end"] + CODING_PARENTAL_EXPAND_BP
    non_coding_nonoverlap_df = (
        BedTool.from_dataframe(non_coding_rna_bed_df)
        .subtract(BedTool.from_dataframe(expanded))
        .to_dataframe(disable_auto_names=True, header=None)
    )
    non_coding_nonoverlap_df.columns = non_coding_rna_bed_df.columns.tolist()
    out_ncrna_nooverlap = Path(config.out_ncrna).parent / "non_coding_rna_without_overlap_with_coding_gene.bed"
    non_coding_nonoverlap_df.to_csv(out_ncrna_nooverlap, sep="\t", index=False)
    logger.success(f"Saved non-coding RNA (no coding overlap) → {out_ncrna_nooverlap}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Set and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse GFF3 annotation to extract genomic region BED files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python extract_genome_region.py \\
      --gff Schizosaccharomyces_pombe_all_chromosomes.gff3 \\
      --fai Schizosaccharomyces_pombe_all_chromosomes.fa.fai \\
      --peptide_stats PeptideStats.tsv \\
      --gene_ids gene_IDs_names_products.tsv \\
      --fypo FYPOviability.tsv \\
      --hayles Hayles_2013_OB_merged_categories.xlsx \\
      --out_primary coding_gene_primary_transcripts.bed \\
      --out_intergenic intergenic_regions.bed \\
      --out_ncrna non_coding_rna.bed \\
      --out_genome_intervals Genome_intervals.bed \\
      --out_overlapped overlapped_region.bed
        """,
    )
    parser.add_argument("--gff", type=Path, required=True, help="GFF3 annotation file")
    parser.add_argument("--fai", type=Path, required=True, help="FASTA index (.fai) for chromosome sizes")
    parser.add_argument("--peptide_stats", type=Path, required=True, help="PomBase PeptideStats TSV")
    parser.add_argument("--gene_ids", type=Path, required=True, help="gene_IDs_names_products TSV")
    parser.add_argument("--fypo", type=Path, required=True, help="FYPOviability TSV")
    parser.add_argument("--hayles", type=Path, required=True, help="Hayles 2013 viability XLSX")
    parser.add_argument("--out_primary", type=Path, required=True, help="Output: coding gene primary transcripts BED")
    parser.add_argument("--out_intergenic", type=Path, required=True, help="Output: intergenic regions BED")
    parser.add_argument("--out_ncrna", type=Path, required=True, help="Output: non-coding RNA BED")
    parser.add_argument("--out_genome_intervals", type=Path, required=True, help="Output: genome intervals BED")
    parser.add_argument("--out_overlapped", type=Path, required=True, help="Output: overlapping coding regions BED")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main entry point of the script."""
    args = parse_args()
    setup_logging("DEBUG" if args.verbose else "INFO")

    try:
        config = Config(
            gff_file=args.gff,
            fai_file=args.fai,
            peptide_stats_file=args.peptide_stats,
            gene_ids_file=args.gene_ids,
            fypo_file=args.fypo,
            hayles_file=args.hayles,
            out_primary=args.out_primary,
            out_intergenic=args.out_intergenic,
            out_ncrna=args.out_ncrna,
            out_genome_intervals=args.out_genome_intervals,
            out_overlapped=args.out_overlapped,
        )

        logger.info(f"Extracting genome regions from {config.gff_file.name}")
        run_pipeline(config)
        logger.success("Script completed successfully!")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
