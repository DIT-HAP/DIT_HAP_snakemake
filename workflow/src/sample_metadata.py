"""Sample metadata extraction from pipeline file naming conventions.

Some pipeline outputs carry no sample/timepoint/condition columns, so the
metadata survives only in the filename stem. Two conventions exist:

- ``{{sample}}_{{timepoint}}_{{condition}}.tsv`` — per-(sample, timepoint,
  condition) tables from ``merge_strand_insertions``;
- ``{{sample}}_{{condition}}.*.tsv`` — per-(sample, condition) concatenations
  (e.g. ``11_merged``), where the sample is everything before the first dot.

This module is the single place that knows these conventions.

If the upstream schema later gains explicit sample/timepoint/condition columns,
callers should prefer those and this module becomes obsolete.
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

from loguru import logger


# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def parse_filename(file_path: Path) -> tuple[str, str, str] | None:
    """Parse filename stem into sample, timepoint, condition. Return None if format is invalid."""
    stem = file_path.stem
    parts = stem.split('_')

    if len(parts) != 3:
        logger.warning(
            f"Filename {file_path.name} does not match expected pattern "
            f"{{sample}}_{{timepoint}}_{{condition}}.tsv (got {len(parts)} parts)"
        )
        return None

    sample, timepoint, condition = parts
    return sample, timepoint, condition


@logger.catch
def parse_sample_name(file_path: Path) -> str:
    """Derive the sample label from a filename by dropping every dotted suffix."""
    return file_path.name.split(".")[0]

