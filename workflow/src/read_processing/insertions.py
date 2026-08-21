"""Insertion-site extraction and counting from aligned-read tables.

Identifies and counts transposon insertion sites from the columnar read
alignment table produced by BAM parsing. Each aligned read carries a genomic
chromosome, reference start/end coordinates, and a strand orientation; the
insertion coordinate is derived from the strand-specific position of the TTAA
target motif.

Rows with a valid ``+``/``-`` strand and complete coordinates are kept, then
the insertion coordinate is computed per row: for the ``+`` strand the site is
``R1_Ref_Start + 4`` (the position immediately after ``TTAA``), and for the
``-`` strand it is ``R1_Ref_End``. Rows are grouped by
``(chromosome, coordinate, strand)`` and counted into a per-chunk
``InsertionCounts`` mapping.

``extract_insertion_sites.py``'s ``process_chunks``, ``write_empty_output``,
and ``count_insertion_sites`` are the chunk loop, reader, and writer and stay
in that script per Phase 4 convention 3.

Input
-----
- A chunk (DataFrame) of aligned reads with columns ``R1_Strand``,
  ``R1_Chrom``, ``R1_Ref_Start``, and ``R1_Ref_End``.

Output
------
- ``InsertionCounts`` mapping and ``ExtractionStats`` dataclass, and an output
  DataFrame with columns ``Chr``, ``Coordinate``, ``+``, ``-``.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-07-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================
type InsertionCounts = dict[tuple[str, int], dict[str, int]]

# =============================================================================
# DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class ExtractionStats:
    """Results of the insertion-site extraction."""
    total_rows: int
    valid_rows: int
    invalid_rows: int
    unique_sites: int
    total_plus_insertions: int
    total_minus_insertions: int

    @property
    def total_insertions(self) -> int:
        """Total number of insertions across both strands."""
        return self.total_plus_insertions + self.total_minus_insertions

    @property
    def validity_rate(self) -> float:
        """Percentage of valid rows."""
        if self.total_rows == 0:
            return 0.0
        return (self.valid_rows / self.total_rows) * 100

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def calculate_insertion_coordinates_vectorized(valid_df: pd.DataFrame) -> pd.Series:
    """Calculate insertion coordinates for all rows using vectorized operations."""
    # For + strand: ref_start + 4; for - strand: ref_end
    plus_coords = valid_df['R1_Ref_Start'].astype('Int64') + 4
    minus_coords = valid_df['R1_Ref_End'].astype('Int64')

    # Select based on strand with np.where
    plus_mask = valid_df['R1_Strand'] == '+'
    coordinates = np.where(plus_mask, plus_coords, minus_coords)

    return pd.Series(coordinates, index=valid_df.index, dtype='Int64')


@logger.catch
def create_validation_mask(df: pd.DataFrame) -> pd.Series:
    """Create a validation mask for filtering valid rows."""
    return (
        df['R1_Strand'].notna() &
        df['R1_Chrom'].notna() &
        df['R1_Ref_Start'].notna() &
        df['R1_Ref_End'].notna() &
        df['R1_Strand'].isin(['+', '-'])
    )


@logger.catch
def count_insertions_vectorized(valid_df: pd.DataFrame) -> InsertionCounts:
    """Count insertions using vectorized operations."""
    # Calculate coordinates for all valid rows at once (vectorized)
    valid_df = valid_df.copy()
    valid_df['Insertion_Coordinate'] = calculate_insertion_coordinates_vectorized(valid_df)

    # Group and count using pandas operations
    grouped = valid_df.groupby(['R1_Chrom', 'Insertion_Coordinate', 'R1_Strand']).size()

    # Convert to our dictionary format
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    for (chrom, coord, strand), count in grouped.items():
        insertion_counts[(chrom, int(coord))][strand] = count

    return dict(insertion_counts)


@logger.catch
def extract_insertion_sites(chunk: pd.DataFrame, chunk_num: int) -> tuple[InsertionCounts, int, int]:
    """Process a single chunk of data to extract insertion sites."""
    chunk_rows = len(chunk)

    # Filter valid rows
    valid_mask = create_validation_mask(chunk)
    valid_chunk = chunk[valid_mask].copy()
    valid_rows = len(valid_chunk)
    invalid_rows = chunk_rows - valid_rows

    if chunk_num == 1 or chunk_num % 10 == 0:
        retention_rate = (valid_rows / chunk_rows * 100) if chunk_rows > 0 else 0
        logger.info(f"Chunk {chunk_num}: {valid_rows:,}/{chunk_rows:,} valid rows ({retention_rate:.1f}%)")

    # Count insertions using vectorized operations
    insertion_counts = count_insertions_vectorized(valid_chunk) if valid_rows > 0 else {}

    return insertion_counts, valid_rows, invalid_rows


@logger.catch
def create_output_dataframe(insertion_counts: InsertionCounts) -> tuple[pd.DataFrame, int, int]:
    """Create output DataFrame from insertion counts."""
    output_data = []
    total_plus = 0
    total_minus = 0

    for (chrom, coord), strand_counts in insertion_counts.items():
        plus_count = strand_counts['+']
        minus_count = strand_counts['-']

        output_data.append({
            'Chr': chrom,
            'Coordinate': coord,
            '+': plus_count,
            '-': minus_count
        })

        total_plus += plus_count
        total_minus += minus_count

    output_df = pd.DataFrame(output_data)
    output_df = output_df.sort_values(['Chr', 'Coordinate'])

    return output_df, total_plus, total_minus
