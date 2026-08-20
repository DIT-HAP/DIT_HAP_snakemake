"""Table I/O utilities for DIT-HAP insertion and annotation data."""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX: list[str] = ["Chr", "Coordinate", "Strand", "Target"]

# =============================================================================
# CORE LOGIC
# =============================================================================
def read_table(path: Path, **kwargs) -> pd.DataFrame:
    """Read a table into a pandas DataFrame based on file extension.

    Supports TSV, BED, CSV, XLSX. Passes **kwargs to underlying pandas reader.
    """
    suffix = path.suffix.lower()
    name_lower = path.name.lower()

    if suffix == ".tsv" or "tsv" in name_lower:
        return pd.read_csv(path, sep="\t", **kwargs)
    elif suffix == ".bed" or "bed" in name_lower:
        return pd.read_csv(path, sep="\t", **kwargs)
    elif suffix == ".csv" or "csv" in name_lower:
        return pd.read_csv(path, sep=",", **kwargs)
    elif suffix in (".xlsx", ".xls") or "xlsx" in name_lower:
        return pd.read_excel(path, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {path}")


def read_insertion_table(path: Path, **kwargs) -> pd.DataFrame:
    """Read an insertion table with the standard 4-level index.

    Index: Chr, Coordinate, Strand, Target.
    Passes **kwargs to read_table (e.g., usecols).
    """
    if "index_col" not in kwargs:
        kwargs["index_col"] = [0, 1, 2, 3]
    return read_table(path, **kwargs)


def write_table(df: pd.DataFrame, path: Path, **kwargs) -> None:
    """Write a DataFrame to a file based on extension.

    Supports TSV, CSV, XLSX. Passes **kwargs to underlying pandas writer.
    """
    suffix = path.suffix.lower()
    name_lower = path.name.lower()

    if suffix == ".tsv" or "tsv" in name_lower:
        df.to_csv(path, sep="\t", **kwargs)
    elif suffix == ".csv" or "csv" in name_lower:
        df.to_csv(path, sep=",", **kwargs)
    elif suffix in (".xlsx", ".xls"):
        df.to_excel(path, **kwargs)
    else:
        raise ValueError(f"Unsupported file type for write: {path}")
