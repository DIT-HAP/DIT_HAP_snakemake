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
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from pybedtools import BedTool

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_table  # noqa: E402
from gene_metadata import resolve_gene_ids  # noqa: E402
from reference_data.genome_region import (  # noqa: E402
    build_genome_intervals,
    build_intergenic_bed,
    annotate_intergenic_region_flanks,
    find_overlapping_regions,
    gff_features_to_bed,
    parse_gff_data,
    select_primary_transcripts,
)

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
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
    setup_logger("DEBUG" if args.verbose else "INFO")

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

        # --- Parse GFF ---
        gff_df = parse_gff_data(config.gff_file)

        # --- Load auxiliary data ---
        peptide_stats_df = pd.read_csv(config.peptide_stats_file, sep="\t")
        gene_to_peptide_length_map = dict(zip(peptide_stats_df["Systematic_ID"], peptide_stats_df["Residues"]))
        logger.info(f"Loaded peptide statistics for {len(gene_to_peptide_length_map):,} proteins")

        gene_IDs_names_products = read_table(config.gene_ids_file)
        gene_IDs_names_products["gene_name"] = gene_IDs_names_products["gene_name"].fillna(
            gene_IDs_names_products["gene_systematic_id"]
        )
        ID2name = dict(zip(gene_IDs_names_products["gene_systematic_id"], gene_IDs_names_products["gene_name"]))

        FYPOviability_df = pd.read_csv(config.fypo_file, sep="\t", header=None, names=["Systematic ID", "FYPOviability"])
        FYPOviability = FYPOviability_df.set_index("Systematic ID")["FYPOviability"].to_dict()

        Hayles_viability_df = pd.read_excel(config.hayles_file)
        Hayles_viability_df["Updated_Systematic_ID"] = resolve_gene_ids(
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
            # Preserve original run_pipeline behavior: log and return early
            # without raising, so the outer try/except still reports success.
            logger.error("No coding gene BED features produced")
            logger.success("Script completed successfully!")
            return 0
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
        processed_ncrna_beds = []
        for name, group in non_coding_rna_df.groupby("Feature"):
            result = gff_features_to_bed(group, "Non-coding gene", gene_to_peptide_length_map)
            if result is not None and not result.empty:
                processed_ncrna_beds.append(result)
        non_coding_rna_bed_df = (
            pd.concat(processed_ncrna_beds).reset_index(drop=True)
            if processed_ncrna_beds
            else pd.DataFrame()
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

        logger.success("Script completed successfully!")
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
