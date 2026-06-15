"""
Concatenate timepoint insertion count data with target sequence annotation.

Concatenates insertion count data across multiple timepoints for the same sample,
adding target sequence information from the reference genome. Generates separate
output files for PBL, PBR, and total read counts with target sequence annotation.

Typical Usage:
    python concatenate_timepoint_data.py -i file1.tsv file2.tsv -tp T1 T2 -g genome.fasta -o output_pbl.tsv --output_pbr output_pbr.tsv --output_reads output_reads.tsv

Input: Multiple timepoint insertion count files, reference genome FASTA
Output: Three output files (PBL, PBR, Reads) with concatenated timepoint data
Adds target sequence information from reference genome for each insertion site.
"""

# =============================== Imports ===============================
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple
from pydantic import BaseModel, Field, field_validator
import pandas as pd

# The following is for Bio sequence handling
from Bio import SeqIO


# =============================== Constants ===============================
# No plotting constants needed for this script

# =============================== Configuration & Models ===============================
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_files: List[Path] = Field(..., description="Input insertion count files", min_items=1)
    genome_file: Path = Field(..., description="Reference genome FASTA file")
    output_pbl: Path = Field(..., description="Output PBL file path")
    output_pbr: Path = Field(..., description="Output PBR file path") 
    output_reads: Path = Field(..., description="Output reads file path")
    sample_name: str = Field(..., description="Sample identifier", min_length=1)
    timepoints: List[str] = Field(..., description="Timepoint names", min_items=1)

    @field_validator('input_files')
    def validate_input_files(cls, v):
        missing = [f for f in v if not f.exists()]
        if missing:
            raise ValueError(f"Input files not found: {missing}")
        return v
    
    @field_validator('genome_file')
    def validate_genome_file(cls, v):
        if not v.exists():
            raise ValueError(f"Genome file not found: {v}")
        return v
    
    @field_validator('output_pbl', 'output_pbr', 'output_reads')
    def validate_output_dir(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True)
        return v
    
    @field_validator('timepoints')
    def validate_timepoint_count(cls, v, info):
        if hasattr(info, 'data') and 'input_files' in info.data:
            input_files = info.data['input_files']
            if len(v) != len(input_files):
                raise ValueError(
                    f"Number of timepoints ({len(v)}) must match number of input files "
                    f"({len(input_files)})"
                )
        return v
    
    class Config:
        frozen = True

# REPLACE THIS ENTIRE CLASS WITH A MODEL FOR YOUR SPECIFIC RESULTS.
class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    num_timepoints: int = Field(..., ge=1, description="Number of timepoints")
    num_insertions: int = Field(..., ge=0, description="Total unique insertions")
    num_chromosomes: int = Field(..., ge=0, description="Number of chromosomes")
    total_pbl_reads: int = Field(..., ge=0, description="Total PBL reads")
    total_pbr_reads: int = Field(..., ge=0, description="Total PBR reads")
    total_reads: int = Field(..., ge=0, description="Total combined reads")
    
    class Config:
        frozen = True


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

# =============================== Core Functions ===============================

@logger.catch
def load_reference_data(genome_path: Path) -> Dict:
    """Load reference genome sequences from FASTA file for target sequence extraction."""
    logger.info(f"Loading reference genome: {genome_path}")
    
    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        logger.success(f"Loaded {len(ref_dict)} sequences from genome")
        
        # Log chromosome names
        chroms = list(ref_dict.keys())[:5]
        if len(ref_dict) > 5:
            logger.debug(f"Chromosomes: {chroms} ... and {len(ref_dict) - 5} more")
        else:
            logger.debug(f"Chromosomes: {chroms}")
        
        return ref_dict
        
    except Exception as e:
        logger.error(f"Error loading genome: {e}")
        raise


@logger.catch
def extract_target_sequence(
    chrom: str,
    coordinate: int,
    ref_dict: Dict
) -> str:
    """Extract 4bp target sequence from reference genome at insertion coordinate."""
    try:
        if chrom not in ref_dict:
            logger.warning(f"Chromosome {chrom} not found in reference")
            return "NNNN"
        
        # Extract 4bp target sequence
        # Coordinate is 1-based, convert to 0-based
        start = coordinate - 4
        end = coordinate
        
        if start < 0:
            logger.warning(f"Coordinate {coordinate} too close to chromosome start")
            return "NNNN"
        
        seq = str(ref_dict[chrom].seq[start:end])
        
        if len(seq) != 4:
            logger.warning(f"Could not extract 4bp at {chrom}:{coordinate}")
            return "NNNN"
        
        return seq.upper()
        
    except Exception as e:
        logger.warning(f"Error extracting target at {chrom}:{coordinate}: {e}")
        return "NNNN"


@logger.catch
def process_concatenation_data(
    config: InputOutputConfig,
    ref_dict: Dict
) -> Tuple[pd.DataFrame, AnalysisResult]:
    """Concatenate insertion data across multiple timepoints with target sequence annotation."""
    logger.info(f"Concatenating {len(config.timepoints)} timepoints")
    
    tp_files = {}
    for tp in config.timepoints:
        for file in config.input_files:
            if f"_{tp}_" in file.name:
                tp_files[tp] = file

    # Load all timepoint files
    dfs = {}
    for tp, file in tp_files.items():
        logger.debug(f"Loading timepoint {tp} from {file}")
        
        try:
            df = pd.read_csv(file, header=0, index_col=[0, 1, 2], sep="\t")
            logger.debug(f"  Loaded {len(df)} insertions for {tp}")
            dfs[tp] = df
        except Exception as e:
            logger.error(f"Error loading {file}: {e}")
            raise
    
    # Concatenate all timepoints
    logger.info("Concatenating dataframes...")
    concatenated = pd.concat(dfs, axis=1, join="outer")
    
    # Sort by timepoint names and coordinates
    concatenated = concatenated.sort_index(
        level=0, axis=1, key=lambda x: x.str.lower()
    ).sort_index(axis=0)
    
    logger.success(f"Concatenated {len(concatenated)} unique insertion sites")
    
    # Add target sequence information
    logger.info("Adding target sequences...")
    target_sequences = []
    
    for idx in concatenated.index:
        chrom = idx[0]
        coordinate = idx[1]
        target = extract_target_sequence(chrom, coordinate, ref_dict)
        target_sequences.append(target)
    
    # Add target as new index level
    concatenated = concatenated.set_index(
        pd.Series(target_sequences, name="Target", index=concatenated.index),
        append=True
    )
    
    # Count unique targets
    unique_targets = concatenated.index.get_level_values("Target").unique()
    logger.info(f"Found {len(unique_targets)} unique target sequences")
    
    # Log target distribution if interesting
    target_counts = concatenated.index.get_level_values("Target").value_counts()
    if "TTAA" in target_counts.index:
        ttaa_fraction = target_counts["TTAA"] / len(concatenated) * 100
        logger.info(f"TTAA targets: {target_counts['TTAA']} ({ttaa_fraction:.1f}%)")
    
    # Calculate read totals before creating frozen stats object
    total_pbl_reads = 0
    total_pbr_reads = 0
    total_reads = 0
    
    if "PBL" in concatenated.columns.get_level_values(1):
        pbl_data = concatenated.xs("PBL", level=1, axis=1)
        total_pbl_reads = int(pbl_data.sum().sum())
    
    if "PBR" in concatenated.columns.get_level_values(1):
        pbr_data = concatenated.xs("PBR", level=1, axis=1)
        total_pbr_reads = int(pbr_data.sum().sum())
    
    if "Reads" in concatenated.columns.get_level_values(1):
        reads_data = concatenated.xs("Reads", level=1, axis=1)
        total_reads = int(reads_data.sum().sum())
    
    # Create statistics object with calculated values (can't modify after creation due to frozen=True)
    result = AnalysisResult(
        num_timepoints=len(config.timepoints),
        num_insertions=len(concatenated),
        num_chromosomes=concatenated.index.get_level_values("Chr").nunique(),
        total_pbl_reads=total_pbl_reads,
        total_pbr_reads=total_pbr_reads,
        total_reads=total_reads
    )
    
    return concatenated, result


@logger.catch
def save_processed_data(
    concatenated: pd.DataFrame,
    config: InputOutputConfig
) -> None:
    """Save concatenated timepoint data to separate output files for each read type."""
    logger.info("Saving concatenated data...")
    
    # Save PBL data
    if "PBL" in concatenated.columns.get_level_values(1):
        pbl_data = concatenated.xs("PBL", level=1, axis=1)
        pbl_data.fillna(0).astype(int).to_csv(
            config.output_pbl, index=True, sep="\t"
        )
        logger.success(f"Saved PBL data to {config.output_pbl}")
    else:
        logger.warning("No PBL data found in concatenated results")
    
    # Save PBR data
    if "PBR" in concatenated.columns.get_level_values(1):
        pbr_data = concatenated.xs("PBR", level=1, axis=1)
        pbr_data.fillna(0).astype(int).to_csv(
            config.output_pbr, index=True, sep="\t"
        )
        logger.success(f"Saved PBR data to {config.output_pbr}")
    else:
        logger.warning("No PBR data found in concatenated results")
    
    # Save Reads data
    if "Reads" in concatenated.columns.get_level_values(1):
        reads_data = concatenated.xs("Reads", level=1, axis=1)
        reads_data.fillna(0).astype(int).to_csv(
            config.output_reads, index=True, sep="\t"
        )
        logger.success(f"Saved Reads data to {config.output_reads}")
    else:
        logger.warning("No Reads data found in concatenated results")


# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments for concatenating timepoint insertion data."""
    parser = argparse.ArgumentParser(description="Concatenate insertion data across timepoints with target sequence annotation")
    parser.add_argument("-s", "--sample", type=str, required=True, help="Sample name")
    parser.add_argument("-i", "--input", type=Path, nargs="+", required=True, help="Path to the input insertion count files")
    parser.add_argument("-tp", "--timepoints", type=str, nargs="+", required=True, help="Timepoint names")
    parser.add_argument("-g", "--genome", type=Path, required=True, help="Reference genome FASTA file")
    parser.add_argument("-ol", "--output_pbl", type=Path, required=True, help="Output PBL file path")
    parser.add_argument("-or", "--output_pbr", type=Path, required=True, help="Output PBR file path")
    parser.add_argument("-o", "--output_reads", type=Path, required=True, help="Output reads file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point for concatenating timepoint insertion data with target sequence annotation."""

    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
    try:
        config = InputOutputConfig(
            input_files=args.input,
            timepoints=args.timepoints,
            genome_file=args.genome,
            output_pbl=args.output_pbl,
            output_pbr=args.output_pbr,
            output_reads=args.output_reads,
            sample_name=args.sample
        )

        logger.info(f"Starting processing of {config.sample_name}")

        # Load reference genome
        ref_dict = load_reference_data(config.genome_file)
        
        # Run the core analysis/logic
        concatenated, results = process_concatenation_data(config, ref_dict)
        
        # Save results
        save_processed_data(concatenated, config)
        
        logger.success(f"Analysis complete. Results saved to {config.output_pbl}, {config.output_pbr}, and {config.output_reads}")
        
        # Display summary statistics
        logger.info("=" * 60)
        logger.info("CONCATENATION SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Sample: {config.sample_name}")
        logger.info(f"Timepoints: {', '.join(config.timepoints)}")
        logger.success(f"Unique insertion sites: {results.num_insertions:,}")
        logger.info(f"Chromosomes: {results.num_chromosomes}")
        if results.total_reads > 0:
            logger.info("\nRead counts:")
            logger.info(f"  PBL reads: {results.total_pbl_reads:,}")
            logger.info(f"  PBR reads: {results.total_pbr_reads:,}")
            logger.info(f"  Total reads: {results.total_reads:,}")
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()