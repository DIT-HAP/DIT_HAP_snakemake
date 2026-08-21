"""Genome region extraction from PomBase GFF3 annotation.

Extracted from ``workflow/scripts/reference_data/extract_genome_region.py``.

Parse a PomBase GFF3 annotation to identify primary coding transcripts (by
matching accumulated CDS length to known peptide lengths), derive intergenic
regions as the complement of primary transcript spans, extract non-coding RNA
intervals, annotate flanking orientation for intergenic regions, find
overlapping coding regions via BedTools self-intersection, and build the
combined genome-interval BED.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pybedtools
from loguru import logger
from pybedtools import BedTool

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
CHR_ORDER = ["chr_II_telomeric_gap", "I", "II", "III", "mating_type_region", "mitochondrial"]
TRANSCRIPT_FEATURE_TYPES = frozenset(["mRNA", "tRNA", "rRNA", "snoRNA", "snRNA", "lncRNA"])

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
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
    extract_systematic_ID_pattern = re.compile(r"ID=(\S+?)(?:$|(?:|\.\d(?::\S+|));)")
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
