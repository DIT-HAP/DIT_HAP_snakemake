"""Gene ID resolution utilities for PomBase systematic IDs."""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
from loguru import logger

# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def resolve_gene_ids(
    genes: list[str],
    gene_metadata: pd.DataFrame,
    gene_filter: str = "gene_type == 'protein coding gene'",
) -> list[str | float]:
    """Resolve gene identifiers to current PomBase systematic IDs.

    Looks up each gene in three mappings: current systematic IDs, gene names,
    and synonyms. Returns updated IDs where found, original identifiers where
    not, and np.nan for ambiguous matches.

    Args:
        genes: List of gene identifiers to resolve (systematic IDs, names, or synonyms).
        gene_metadata: DataFrame with columns gene_systematic_id, gene_name, synonyms, gene_type.
        gene_filter: Query string to filter genes (default: protein-coding only).

    Returns:
        List of resolved systematic IDs (same length as input).
    """
    # Normalize gene_name column
    gene_metadata = gene_metadata.copy()
    gene_metadata["gene_name"] = gene_metadata["gene_name"].fillna(
        gene_metadata["gene_systematic_id"]
    )

    # Filter and create mappings
    filtered = gene_metadata.query(gene_filter)

    synonyms2ID = (
        filtered.set_index("gene_systematic_id")["synonyms"]
        .str.split(",")
        .explode()
        .str.strip()
        .dropna()
        .reset_index()
        .set_index("synonyms")
    )

    names2ID = (
        filtered.set_index("gene_name")["gene_systematic_id"]
        .drop_duplicates()
        .reset_index()
        .set_index("gene_name")
    )

    sysIDs_now = filtered["gene_systematic_id"].unique().tolist()

    # Resolve each gene
    updated = []
    for gene in genes:
        # Normalize case for split IDs (e.g., "SPAC1.02" → "SPAC1.02")
        if isinstance(gene, str):
            gene = gene.strip()
            if "." in gene:
                parts = gene.split(".")
                gene = parts[0].upper() + "." + parts[1].lower()

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
                logger.warning(f"{gene} has multiple name updates: {val.tolist()}")
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
