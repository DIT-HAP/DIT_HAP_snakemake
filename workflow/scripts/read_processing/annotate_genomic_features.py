"""
Enhanced genomic feature annotation with Pydantic validation and Loguru logging.

Annotates transposon insertion sites with genomic features including genes,
intergenic regions, and coding sequences. Calculates distances to start/stop
codons and determines affected amino acid residues.

Typical Usage:
    python annotate_genomic_features.py --input [input_file] --genome-region [genome_file] --output [output_file]

Input: Input insertion file (TSV or CSV)
Output: Output annotation file with genomic features
Additional: Genome region BED file for feature annotation
"""

# =============================== Imports ===============================
# import necessary libraries, just keep the ones you need, the following are just examples
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd
from pybedtools import BedTool


# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_file: Path = Field(..., description="Path to the input insertion file")
    genome_region_file: Path = Field(..., description="Path to the genome region BED file")
    output_file: Path = Field(..., description="Path to the output annotation file")

    @field_validator('input_file', 'genome_region_file')
    def validate_input_file(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True) # Create dir if it doesn't exist
        return v
    
    class Config:
        frozen = True # Makes the model immutable after creation


class GenomicFeature(BaseModel):
    """Model for a genomic feature annotation."""
    
    chromosome: str = Field(..., min_length=1, description="Chromosome")
    coordinate: int = Field(..., ge=0, description="Insertion coordinate")
    strand: str = Field(..., pattern="^[+-]$", description="Strand")
    target: str = Field(..., description="Target sequence")
    feature_type: str = Field(..., description="Feature type")
    gene_name: Optional[str] = Field(None, description="Gene name if applicable")
    distance_to_start: Optional[float] = Field(None, description="Distance to start codon")
    distance_to_stop: Optional[float] = Field(None, description="Distance to stop codon")
    fraction_to_start: Optional[float] = Field(None, ge=0, le=1, description="Fraction to start")
    fraction_to_stop: Optional[float] = Field(None, ge=0, le=1, description="Fraction to stop")
    residue_affected: Optional[int] = Field(None, ge=1, description="Affected amino acid")
    residue_frame: Optional[int] = Field(None, ge=0, le=2, description="Reading frame")
    insertion_direction: Optional[str] = Field(None, description="Forward/Reverse")
    
    class Config:
        frozen = True


class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    total_insertions: int = Field(..., ge=0, description="Total number of insertions processed")
    annotated_insertions: int = Field(..., ge=0, description="Number of insertions successfully annotated")
    coding_insertions: int = Field(..., ge=0, description="Number of insertions in coding regions")
    intergenic_insertions: int = Field(..., ge=0, description="Number of insertions in intergenic regions")
    unique_genes: int = Field(..., ge=0, description="Number of unique genes affected")
    forward_insertions: int = Field(..., ge=0, description="Number of forward insertions")
    reverse_insertions: int = Field(..., ge=0, description="Number of reverse insertions")
    
    @property
    def coding_percentage(self) -> float:
        """Percentage of insertions in coding regions."""
        if self.total_insertions == 0:
            return 0.0
        return (self.coding_insertions / self.total_insertions) * 100.0


# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove() # Remove default logger
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Annotation Functions ========================

@logger.catch
def load_insertion_data(input_path: Path) -> pd.DataFrame:
    """
    Load insertion data from TSV or CSV file.
    
    Args:
        input_path: Path to insertion file
        
    Returns:
        DataFrame with insertion data
    """
    logger.info(f"Loading insertion data from {input_path}")
    
    file_ext = input_path.suffix.lower()
    
    try:
        if file_ext == ".tsv":
            df = pd.read_csv(input_path, header=0, usecols=[0, 1, 2, 3], sep="\t")
        elif file_ext == ".csv":
            df = pd.read_csv(input_path, header=0, usecols=[0, 1, 2, 3])
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        # Rename coordinate column for BED processing
        df = df.rename(columns={"Coordinate": "End"})
        df.insert(1, "Start", df["End"])
        
        logger.success(f"Loaded {len(df)} insertion sites")
        logger.debug("For BED processing, Start and End are the same (point insertions)")
        
        return df
        
    except Exception as e:
        logger.error(f"Error loading insertion data: {e}")
        raise


@logger.catch
def load_genome_regions(region_path: Path) -> pd.DataFrame:
    """
    Load genome region annotations.
    
    Args:
        region_path: Path to genome region BED file
        
    Returns:
        DataFrame with genome regions
    """
    logger.info(f"Loading genome regions from {region_path}")
    
    try:
        df = pd.read_csv(region_path, sep="\t", header=0)
        df = df.rename(columns={"#Chr": "Chr"})
        
        logger.success(f"Loaded {len(df)} genome regions")
        
        # Log region types if available
        if "Type" in df.columns:
            type_counts = df["Type"].value_counts()
            logger.debug("Region types:")
            for region_type, count in type_counts.items():
                logger.debug(f"  {region_type}: {count}")
        
        return df
        
    except Exception as e:
        logger.error(f"Error loading genome regions: {e}")
        raise


@logger.catch
def calculate_codon_distances(row: pd.Series) -> pd.Series:
    """
    Calculate distances to start and stop codons.
    
    Args:
        row: DataFrame row with annotation data
        
    Returns:
        Series with calculated distances
    """
    name_distance = [
        "Distance_to_start_codon",
        "Distance_to_stop_codon",
        "Fraction_to_start_codon",
        "Fraction_to_stop_codon"
    ]
    
    if row["Type"] != "Intergenic region":
        if row["Strand_Interval"] == "+":
            distance_values = [
                row["Distance_to_region_start"],
                row["Distance_to_region_end"],
                row["Fraction_to_region_start"],
                row["Fraction_to_region_end"]
            ]
        else:
            distance_values = [
                row["Distance_to_region_end"],
                row["Distance_to_region_start"],
                row["Fraction_to_region_end"],
                row["Fraction_to_region_start"]
            ]
    else:
        distance_values = [np.nan] * 4
    
    return pd.Series(distance_values, index=name_distance)


@logger.catch
def calculate_affected_residue(row: pd.Series) -> pd.Series:
    """
    Calculate affected amino acid residue and reading frame.
    
    Args:
        row: DataFrame row with annotation data
        
    Returns:
        Series with residue information
    """
    residue_stat = ["Residue_affected", "Residue_frame"]
    
    if row["Type"] in ["Intergenic region", "Non-coding gene"]:
        return pd.Series([np.nan, np.nan], index=residue_stat)
    
    cds_base = float(row.get("Accumulated_CDS_bases", 0))
    
    if row["Feature"] == "CDS":
        if row["Strand_Interval"] == "+":
            cds_base += int(row["Coordinate"]) - int(row["Start_Interval"])
        else:
            cds_base += int(row["End_Interval"] - row["Coordinate"])
    
    residue_affected = cds_base // 3 + 1
    residue_frame = cds_base % 3
    
    return pd.Series([residue_affected, residue_frame], index=residue_stat)


@logger.catch
def assign_insertion_direction(row: pd.Series) -> str:
    """
    Determine insertion direction relative to gene.
    
    Args:
        row: DataFrame row with annotation data
        
    Returns:
        Insertion direction (Forward/Reverse/NA)
    """
    if row["Type"] == "Intergenic region":
        return np.nan
    
    if row["Strand"] == row["Strand_Interval"]:
        return "Forward"
    else:
        return "Reverse"


@logger.catch
def drop_boundary_duplicates(sub_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove duplicated insertions at gene boundaries.
    
    Args:
        sub_df: DataFrame subset for a single insertion
        
    Returns:
        Deduplicated DataFrame
    """
    # Remove insertions exactly at boundaries
    boundary_mask = (
        (sub_df["Distance_to_start_codon"] == 0) |
        (sub_df["Distance_to_stop_codon"] == 0)
    )
    sub_df = sub_df[~boundary_mask]
    
    # If only one type remains, return it
    if sub_df["Type"].nunique() == 1:
        return sub_df
    
    # Otherwise, prefer coding regions over intergenic
    coding_mask = sub_df["Type"] != "Intergenic region"
    if coding_mask.any():
        return sub_df[coding_mask]
    
    return sub_df


@logger.catch
def annotate_insertions(
    insertions_df: pd.DataFrame,
    regions_df: pd.DataFrame
) -> Tuple[pd.DataFrame, AnalysisResult]:
    """
    Annotate insertions with genomic features.
    
    Args:
        insertions_df: DataFrame with insertion sites
        regions_df: DataFrame with genome regions
        
    Returns:
        Tuple of (annotated DataFrame, statistics)
    """
    logger.info("Annotating insertions with genomic features...")
    
    # Convert to BedTool objects
    insertions_bed = BedTool.from_dataframe(insertions_df)
    regions_bed = BedTool.from_dataframe(regions_df)
    
    # Get column names for merging
    insertion_cols = insertions_df.columns.tolist()
    region_cols = [
        f"{col}_Interval" if col in insertion_cols else col
        for col in regions_df.columns
    ]
    
    logger.debug(f"Insertion columns: {len(insertion_cols)}")
    logger.debug(f"Region columns: {len(region_cols)}")
    
    # Intersect insertions with regions
    intersected = insertions_bed.intersect(regions_bed, wa=True, wb=True)
    
    # Convert back to DataFrame
    annotated_df = intersected.to_dataframe(names=insertion_cols + region_cols)
    
    # Clean up DataFrame
    annotated_df.drop(columns=["Start"], inplace=True)
    annotated_df.rename(columns={"End": "Coordinate"}, inplace=True)
    
    # Replace "." with NaN
    annotated_df.replace(r"^\.$", np.nan, inplace=True, regex=True)
    
    logger.info("Calculating distances and annotations...")
    
    # Calculate distances to region boundaries
    annotated_df["Distance_to_region_start"] = (
        annotated_df["Coordinate"] - annotated_df["ParentalRegion_start"]
    )
    annotated_df["Distance_to_region_end"] = (
        annotated_df["ParentalRegion_end"] - annotated_df["Coordinate"]
    )
    annotated_df["Fraction_to_region_start"] = (
        annotated_df["Distance_to_region_start"] /
        annotated_df["ParentalRegion_length"]
    )
    annotated_df["Fraction_to_region_end"] = (
        annotated_df["Distance_to_region_end"] /
        annotated_df["ParentalRegion_length"]
    )
    
    # Calculate codon distances
    codon_distances = annotated_df.apply(calculate_codon_distances, axis=1)
    annotated_df = pd.concat([annotated_df, codon_distances], axis=1)
    
    # Calculate affected residues
    residue_info = annotated_df.apply(calculate_affected_residue, axis=1)
    annotated_df = pd.concat([annotated_df, residue_info], axis=1)
    
    # Assign insertion direction
    annotated_df["Insertion_direction"] = annotated_df.apply(
        assign_insertion_direction, axis=1
    )
    
    logger.info("Removing boundary duplicates...")
    
    # Drop duplicated insertions at boundaries
    annotated_df = (
        annotated_df.groupby(["Chr", "Coordinate", "Strand"])
        .apply(drop_boundary_duplicates)
        .reset_index(drop=True)
    )
    
    # Calculate statistics
    stats = AnalysisResult(
        total_insertions=len(insertions_df),
        annotated_insertions=len(annotated_df),
        coding_insertions=len(annotated_df[annotated_df["Type"] != "Intergenic region"]),
        intergenic_insertions=len(annotated_df[annotated_df["Type"] == "Intergenic region"]),
        unique_genes=annotated_df["GeneName"].nunique() if "GeneName" in annotated_df.columns else 0,
        forward_insertions=len(annotated_df[annotated_df["Insertion_direction"] == "Forward"]),
        reverse_insertions=len(annotated_df[annotated_df["Insertion_direction"] == "Reverse"])
    )
    
    return annotated_df, stats


@logger.catch
def save_annotations(
    annotated_df: pd.DataFrame,
    output_path: Path,
    stats: AnalysisResult
) -> None:
    """
    Save annotated insertions to file.
    
    Args:
        annotated_df: Annotated DataFrame
        output_path: Output file path
        stats: Annotation statistics
    """
    logger.info(f"Saving annotations to {output_path}")
    
    # Save with proper formatting
    annotated_df.to_csv(
        output_path,
        index=False,
        header=True,
        float_format="%.3f",
        sep="\t"
    )
    
    # Display statistics
    logger.info("=" * 60)
    logger.info("ANNOTATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total insertions: {stats.total_insertions:,}")
    logger.success(f"Annotated insertions: {stats.annotated_insertions:,}")
    
    logger.info("\nFeature distribution:")
    logger.info(f"  Coding regions: {stats.coding_insertions:,} ({stats.coding_percentage:.1f}%)")
    logger.info(f"  Intergenic regions: {stats.intergenic_insertions:,}")
    
    if stats.unique_genes > 0:
        logger.info(f"\nGene impact:")
        logger.info(f"  Unique genes affected: {stats.unique_genes:,}")
        logger.info(f"  Forward insertions: {stats.forward_insertions:,}")
        logger.info(f"  Reverse insertions: {stats.reverse_insertions:,}")
    
    logger.success(f"Annotations saved to {output_path}")


# =============================== Core Functions ===============================

@logger.catch
def main_processing_function(config: InputOutputConfig) -> AnalysisResult:
    """
    Main processing function that orchestrates the annotation workflow.
    
    Args:
        config: Input/output configuration object
        
    Returns:
        AnalysisResult containing annotation statistics
    """
    logger.info(f"Starting annotation workflow for {config.input_file}")
    logger.info(f"Genome regions: {config.genome_region_file}")
    logger.info(f"Output file: {config.output_file}")
    
    # Load data
    insertions_df = load_insertion_data(config.input_file)
    regions_df = load_genome_regions(config.genome_region_file)
    
    # Annotate insertions
    annotated_df, stats = annotate_insertions(insertions_df, regions_df)
    
    # Save results
    save_annotations(annotated_df, config.output_file, stats)
    
    logger.success("Annotation workflow completed successfully!")
    return stats

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="Annotate insertion sites with genomic features")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input insertion file")
    parser.add_argument("-g", "--genome-region", type=Path, required=True, help="Path to the genome region BED file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output annotation file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point for the script."""

    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    logger.info(f"Pandas version: {pd.__version__}")
    logger.info(f"NumPy version: {np.__version__}")
    
    try:
        # Create and validate configuration
        config = InputOutputConfig(
            input_file=args.input,
            genome_region_file=args.genome_region,
            output_file=args.output
        )
        
        # Run the core analysis/logic
        results = main_processing_function(config)
        
        logger.success(f"Annotation completed successfully! Processed {results.total_insertions} insertions.")
        return 0
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())