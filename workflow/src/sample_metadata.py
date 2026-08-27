"""Sample metadata extraction from pipeline file naming conventions.

The ``merge_strand_insertions`` output carries only a Chr/Coordinate/Strand
index; sample/timepoint/condition metadata survives only in the filename stem,
following the pattern ``{{sample}}_{{timepoint}}_{{condition}}.tsv``. This
module is the single place that knows that convention.

If the upstream schema later gains explicit Sample/Timepoint/Condition columns,
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
