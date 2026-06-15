# ================================ Imports =================================
from pathlib import Path
import pandas as pd
import numpy as np

# ================================= Utility Functions =================================
def read_file(
    file: Path,
    **kwargs,
) -> pd.DataFrame:
    """Read a file into a pandas DataFrame based on file extension."""
    if "tsv" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, sep=",", **kwargs)
    elif "xlsx" in file.name:
        return pd.read_excel(file, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")

def update_sysIDs(genes: list[str], gene_meta_file: Path, gene_filter: str = "gene_type == 'protein coding gene'") -> list[str | float]:
    """Update gene systematic IDs based on synonyms."""

    # Load gene metadata
    gene_meta = read_file(gene_meta_file)
    gene_meta["gene_name"] = gene_meta["gene_name"].fillna(gene_meta["gene_systematic_id"])

    # Create mappings
    filtered_genes = gene_meta.query(gene_filter)
    synonyms2ID = filtered_genes.set_index("gene_systematic_id")["synonyms"].str.split(",").explode().str.strip().dropna().reset_index().set_index("synonyms")
    names2ID = filtered_genes.set_index("gene_name")["gene_systematic_id"].drop_duplicates().reset_index().set_index("gene_name")
    sysIDs_now = filtered_genes["gene_systematic_id"].unique().tolist()

    # Update systematic IDs
    updated_sysIDs = []
    for gene in genes:
        if isinstance(gene, str):
            gene = gene.strip()
            if "." in gene:
                gene = gene.split(".")[0].upper() + "." + gene.split(".")[1].lower()
        else:
            gene = gene
        if pd.isna(gene):
            updated_sysIDs.append(gene)
            print(f"{gene} is NA")
        elif gene in sysIDs_now:
            updated_sysIDs.append(gene)
        elif gene in names2ID.index:
            updated = names2ID.loc[gene, "gene_systematic_id"]
            if isinstance(updated, str):
                updated_sysIDs.append(updated)
                print(f"{gene} is updated to {updated}")
            else:
                updated_sysIDs.append(np.nan)
                print(f"{gene} has multiple updates:", updated)
        elif gene in synonyms2ID.index:
            updated = synonyms2ID.loc[gene, "gene_systematic_id"]
            if isinstance(updated, str):
                updated_sysIDs.append(updated)
                print(f"{gene} is updated to {updated}")
            else:
                updated_sysIDs.append(np.nan)
                print(f"{gene} has multiple updates:", updated)
        else:
            updated_sysIDs.append(gene)
            print(f"{gene} is not found in geneid2symbol or synonyms2ID")
    return updated_sysIDs