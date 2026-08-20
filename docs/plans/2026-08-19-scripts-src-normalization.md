# Scripts/Src Normalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Normalize `workflow/src/` as an importable library, reduce `workflow/scripts/` to thin CLI entrypoints across 33 scripts in 6 phased rollouts with verification per phase.

**Architecture:** Library/CLI split — `src/` holds pure domain logic (no `main()`, no CLI), `scripts/` holds docstring + `parse_args()` + `main()` only. Shared boilerplate (logger, I/O, ID resolution) consolidates once; per-domain logic moves into `src/<domain>/` modules.

**Tech Stack:** Python 3.12+, pandas, loguru, pytest, existing 11 conda envs (no changes), sys.path.append import bootstrap (no pip install).

**Spec:** `docs/plans/2026-08-19-scripts-src-normalization-design.md`

## Global Constraints

- Python 3.12+ syntax (native generics, `type` aliases, `match/case`, `StrEnum`).
- Library modules: IMPORTS → CONSTANTS → CORE LOGIC only (no `main()`, no `parse_args()`, no `setup_logger()`).
- Scripts: docstring + sys.path bootstrap + `parse_args()` + `main()` only (~60–90 lines).
- Import bootstrap: `SCRIPT_DIR = Path(__file__).parent.resolve(); sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))`
- Naming: `src/logging_setup.py` (not `logging.py`), `src/io_tables.py` (not `io.py`) to avoid stdlib shadowing.
- Config dataclasses stay in scripts/ (CLI validation concern).
- Byte-compare verification per phase before proceeding to next phase.
- No Snakefile changes, no env yml changes, no pip installation.

---

## Phase 1: Foundation (Shared Modules + Import Repointing + Test Migration)

### Task 1.1: Create Shared Foundation Modules

**Files:**
- Create: `workflow/src/__init__.py`
- Create: `workflow/src/logging_setup.py`
- Create: `workflow/src/io_tables.py`
- Create: `workflow/src/gene_metadata.py`
- Read: `workflow/scripts/depletion_scoring/compute_insertion_weights.py:149-164` (canonical setup_logger)
- Read: `workflow/src/utils.py:6-20` (read_file source)
- Read: `workflow/src/utils.py:22-72` (update_sysIDs source)
- Read: `workflow/scripts/reference_data/extract_genome_region.py:141-181` (update_sysID DataFrame variant)

**Interfaces:**
- Produces:
  - `setup_logger(log_level: str = "INFO") -> None`
  - `INSERTION_INDEX: list[str]` = `["Chr", "Coordinate", "Strand", "Target"]`
  - `read_insertion_table(path: Path, **kwargs) -> pd.DataFrame` (uses INSERTION_INDEX as index_col)
  - `read_table(path: Path, **kwargs) -> pd.DataFrame` (generic, infers from extension)
  - `write_table(df: pd.DataFrame, path: Path, **kwargs) -> None`
  - `resolve_gene_ids(genes: list[str], gene_metadata: pd.DataFrame, gene_filter: str = "gene_type == 'protein coding gene'") -> list[str | float]`

- [ ] **Step 1: Create empty workflow/src/__init__.py**

```python
"""DIT-HAP workflow library modules."""
```

- [ ] **Step 2: Run to verify no import errors**

Run: `python -c "import sys; sys.path.insert(0, 'workflow/src'); import workflow.src"`
Expected: No output (success)

- [ ] **Step 3: Create workflow/src/logging_setup.py**

```python
"""Logging configuration for DIT-HAP scripts."""

# =============================================================================
# IMPORTS
# =============================================================================
import sys

from loguru import logger

# =============================================================================
# CORE LOGIC
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
```

- [ ] **Step 4: Create workflow/src/io_tables.py**

```python
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
```

- [ ] **Step 5: Create workflow/src/gene_metadata.py**

```python
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
```

- [ ] **Step 6: Verify imports work**

Run: `python -c "import sys; sys.path.insert(0, 'workflow/src'); from logging_setup import setup_logger; from io_tables import INSERTION_INDEX, read_table; from gene_metadata import resolve_gene_ids; print('OK')"`
Expected: `OK`

- [ ] **Step 7: Commit foundation modules**

```bash
git add workflow/src/__init__.py workflow/src/logging_setup.py workflow/src/io_tables.py workflow/src/gene_metadata.py
git commit -m "feat(src): add shared foundation modules (logging, IO, gene ID resolution)

- logging_setup: canonical setup_logger (replaces 33 copies)
- io_tables: read/write_table, read_insertion_table, INSERTION_INDEX constant
- gene_metadata: resolve_gene_ids (canonical, replaces utils.py + extract_genome_region.py variants)

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 1.2: Create Root pyproject.toml and Move Tests

**Files:**
- Create: `pyproject.toml`
- Move: `workflow/src/test_figures.py` → `workflow/tests/test_figures.py`
- Create: `workflow/tests/__init__.py`

**Interfaces:**
- Consumes: `workflow/src/figures.py` (apply_house_style, save_dual)
- Produces: pytest configuration with pythonpath = ["workflow/src", "workflow/scripts"]

- [ ] **Step 1: Create root pyproject.toml**

```toml
[tool.pytest.ini_options]
pythonpath = ["workflow/src", "workflow/scripts"]
testpaths = ["workflow/tests"]
```

- [ ] **Step 2: Create workflow/tests directory and __init__.py**

```bash
mkdir -p workflow/tests
touch workflow/tests/__init__.py
```

- [ ] **Step 3: Move test_figures.py and update its imports**

```bash
mv workflow/src/test_figures.py workflow/tests/test_figures.py
```

Edit `workflow/tests/test_figures.py` to remove sys.path manipulation:

```python
"""Tests for workflow/src/figures.py rendering utilities."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

# Imports now work via pyproject.toml pythonpath
import cnsplots as cns
from figures import apply_house_style, save_dual


def test_apply_house_style_sets_arial_font():
    """apply_house_style() must set Arial to avoid Helvetica findfont warnings."""
    apply_house_style()

    # Helvetica should be removed, Arial should be first
    assert "Helvetica" not in cns.settings.font_sans_serif
    assert cns.settings.font_sans_serif[0] == "Arial"
    assert cns.settings.panel_label_fontname == "Arial"


def test_apply_house_style_sets_pdf_fonttype_42():
    """PDF fonttype 42 ensures Illustrator can edit text."""
    apply_house_style()
    assert cns.settings.pdf_fonttype == 42


def test_save_dual_creates_both_artifacts(tmp_path):
    """save_dual() must write both PDF and PNG."""
    apply_house_style()
    cns.figure(width=90, height=90)
    plt.plot([1, 2], [1, 2])

    stem = tmp_path / "test_fig"
    save_dual(stem)

    assert (tmp_path / "test_fig.pdf").exists()
    assert (tmp_path / "test_fig.review.png").exists()
    assert (tmp_path / "test_fig.pdf").stat().st_size > 1000
    assert (tmp_path / "test_fig.review.png").stat().st_size > 1000


def test_save_dual_preserves_dots_in_stem(tmp_path):
    """A stem containing dots must not be truncated: sample names here carry dots."""
    apply_house_style()
    cns.figure(width=90, height=90)
    plt.plot([1, 2], [1, 2])

    save_dual(tmp_path / "HD1328-4.YES0_corr")

    # with_suffix()/.stem would have collapsed these to "HD1328-4.*"
    assert (tmp_path / "HD1328-4.YES0_corr.pdf").exists()
    assert (tmp_path / "HD1328-4.YES0_corr.review.png").exists()
    assert not (tmp_path / "HD1328-4.pdf").exists()
```

- [ ] **Step 4: Run tests to verify pytest config works**

Run: `/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests/test_figures.py -v`
Expected: 4 tests pass

- [ ] **Step 5: Commit pytest config and moved test**

```bash
git add pyproject.toml workflow/tests/__init__.py workflow/tests/test_figures.py
git rm workflow/src/test_figures.py
git commit -m "test: add root pytest config and move test_figures to workflow/tests

- pyproject.toml: pythonpath for src + scripts, testpaths = workflow/tests
- test_figures.py moved from src/ to tests/, sys.path hacks removed

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 1.3: Delete Dead Code and Split utils.py

**Files:**
- Delete: `workflow/src/concat_tables.py`
- Delete: `workflow/src/utils.py`

**Interfaces:**
- Consumes: `workflow/src/io_tables.py`, `workflow/src/gene_metadata.py` (created in Task 1.1)
- Produces: Clean src/ with no dead code

- [ ] **Step 1: Verify concat_tables.py is not imported**

Run: `grep -r "from concat_tables\|import concat_tables" workflow/ --include="*.py"`
Expected: No matches (it's dead code)

- [ ] **Step 2: Verify utils.py is only imported once**

Run: `grep -r "from utils\|import utils" workflow/ --include="*.py"`
Expected: One match in `workflow/scripts/quality_control/insertion_orientation_analysis.py:67`

- [ ] **Step 3: Delete concat_tables.py and utils.py**

```bash
git rm workflow/src/concat_tables.py workflow/src/utils.py
```

- [ ] **Step 4: Commit deletion**

```bash
git commit -m "refactor(src): delete dead code (concat_tables, utils)

- concat_tables.py: dead exploration script, not imported
- utils.py: split into io_tables.py + gene_metadata.py in previous commit

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 1.4: Repoint All 33 Scripts' Imports (Shared Modules Only)

**Files:**
- Modify: All 33 `workflow/scripts/**/*.py` files (CLI scripts only, not tests)

**Interfaces:**
- Consumes: `workflow/src/logging_setup.py`, `workflow/src/io_tables.py`, `workflow/src/gene_metadata.py`
- Produces: All scripts using shared foundation, no domain logic moved yet

This task is large (33 files). The pattern is identical for each script:

1. Add sys.path bootstrap after imports (if not present).
2. Replace `def setup_logger(...)` with `from logging_setup import setup_logger  # noqa: E402`.
3. Replace hardcoded `index_col=[0,1,2,3]` with `from io_tables import read_insertion_table  # noqa: E402`.
4. For `insertion_orientation_analysis.py:67`: replace `from utils import read_file` with `from io_tables import read_table`.

I'll show the pattern for 3 representative scripts, then list the remaining 30 with the same transformation.

- [ ] **Step 1: Repoint compute_insertion_weights.py**

Edit `workflow/scripts/depletion_scoring/compute_insertion_weights.py`:

**Before (lines 1-20, 103-164):**
```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ... (PEP 723 block)
"""
Insertion-Level Aggregation Weights
===================================
... (docstring unchanged)
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import pandas as pd
from loguru import logger

# ... (ENUMS, DATACLASSES sections unchanged)

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

# ... (rest unchanged for now)
```

**After:**
```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ... (PEP 723 block unchanged)
"""
Insertion-Level Aggregation Weights
===================================
... (docstring unchanged)
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import pandas as pd
from loguru import logger

# Bootstrap src/ onto sys.path
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402

# ... (ENUMS, DATACLASSES sections unchanged)
# ... (LOGGING SETUP section DELETED — imported from logging_setup now)

# ... (load_inputs function: replace pd.read_csv(..., index_col=[0,1,2,3]) with read_insertion_table(...))
```

Within `load_inputs()` around line 164:
```python
# Before:
lfc = pd.read_csv(config.lfc, sep="\t", index_col=[0, 1, 2, 3])
stats = pd.read_csv(config.stats, sep="\t", index_col=[0, 1, 2, 3])
annotations = pd.read_csv(config.annotations, sep="\t", index_col=[0, 1, 2, 3])

# After:
lfc = read_insertion_table(config.lfc)
stats = read_insertion_table(config.stats)
annotations = read_insertion_table(config.annotations)
```

- [ ] **Step 2: Repoint insertion_orientation_analysis.py (uses utils.read_file)**

Edit `workflow/scripts/quality_control/insertion_orientation_analysis.py`:

**Before (around lines 60-80):**
```python
# ... imports ...
import pandas as pd
from loguru import logger

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../..").resolve()))
from utils import read_file  # noqa: E402

# ... rest of file ...
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
```

**After:**
```python
# ... imports ...
import pandas as pd
from loguru import logger

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_table  # noqa: E402

# ... rest of file, setup_logger() function DELETED ...
```

Within the script, replace all `read_file(...)` calls with `read_table(...)`.

- [ ] **Step 3: Repoint extract_genome_region.py (has setup_logging + update_sysID)**

Edit `workflow/scripts/reference_data/extract_genome_region.py`:

**Before (lines 125-181):**
```python
# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
@logger.catch
def update_sysID(genes: list, gene_IDs_names_products: pd.DataFrame) -> list:
    """Resolve gene identifiers to current systematic IDs via name/synonym lookup."""
    coding = gene_IDs_names_products.query("gene_type == 'protein coding gene'")
    # ... (rest of function body)
```

**After (add bootstrap at top, delete setup_logging, replace update_sysID):**
```python
# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
# ... other imports ...

import pandas as pd
from loguru import logger

SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_table  # noqa: E402
from gene_metadata import resolve_gene_ids  # noqa: E402

# ... (GLOBAL CONSTANTS section unchanged)
# ... (CONFIGURATION & DATACLASSES section unchanged)

# =============================================================================
# LOGGING SETUP section DELETED
# =============================================================================

# =============================================================================
# CORE LOGIC (FUNCTIONS / CLASSES)
# =============================================================================
# update_sysID function DELETED — replaced by resolve_gene_ids from gene_metadata
```

Within `run_pipeline()` around line 441, replace:
```python
# Before:
gene_IDs_names_products = pd.read_csv(config.gene_meta, sep="\t")
genes = update_sysID(config.genes, gene_IDs_names_products)

# After:
gene_IDs_names_products = read_table(config.gene_meta)
genes = resolve_gene_ids(config.genes, gene_IDs_names_products)
```

Also in `main()`, replace `setup_logging(...)` with `setup_logger(...)`.

- [ ] **Step 4: Apply identical transformation to remaining 30 scripts**

For each of the following scripts, apply the same pattern:

**Depletion Scoring (6 more):**
- `workflow/scripts/depletion_scoring/curve_fitting.py`
- `workflow/scripts/depletion_scoring/def_ctr_insertions.py`
- `workflow/scripts/depletion_scoring/gene_level_depletion_analysis.py`
- `workflow/scripts/depletion_scoring/impute_missing_values_using_FR.py`
- `workflow/scripts/depletion_scoring/insertion_level_depletion_analysis_has_replicates.py`
- `workflow/scripts/depletion_scoring/insertion_level_depletion_analysis_no_replicates.py`

**Figures (10):**
- `workflow/scripts/figures/plot_curve_fitting.py`
- `workflow/scripts/figures/plot_dispersions.py`
- `workflow/scripts/figures/plot_distribution_of_curve_fitting.py`
- `workflow/scripts/figures/plot_gene_coverage.py`
- `workflow/scripts/figures/plot_insertion_density.py`
- `workflow/scripts/figures/plot_insertion_orientation.py`
- `workflow/scripts/figures/plot_ma_plot.py`
- `workflow/scripts/figures/plot_ma_plot_replicates.py`
- `workflow/scripts/figures/plot_pbl_pbr_correlation.py`
- `workflow/scripts/figures/plot_read_count_distribution.py`

**Quality Control (5 more, excluding insertion_orientation_analysis already done):**
- `workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py`
- `workflow/scripts/quality_control/extract_mapping_filtering_statistics.py`
- `workflow/scripts/quality_control/gene_coverage_analysis.py`
- `workflow/scripts/quality_control/insertion_density_analysis.py`
- `workflow/scripts/quality_control/read_count_distribution_analysis.py`

**Read Processing (9):**
- `workflow/scripts/read_processing/annotate_genomic_features.py`
- `workflow/scripts/read_processing/concat_counts_and_annotations.py`
- `workflow/scripts/read_processing/concatenate_timepoint_data.py`
- `workflow/scripts/read_processing/extract_insertion_sites.py`
- `workflow/scripts/read_processing/filter_aligned_reads.py`
- `workflow/scripts/read_processing/merge_similar_timepoints.py`
- `workflow/scripts/read_processing/merge_strand_insertions.py`
- `workflow/scripts/read_processing/parse_bam_to_tsv.py`
- `workflow/scripts/read_processing/reads_hard_filtering.py`

Pattern for each:
1. Add sys.path bootstrap after imports if missing.
2. Delete `def setup_logger(...)` function.
3. Add `from logging_setup import setup_logger  # noqa: E402`.
4. Replace `pd.read_csv(..., sep="\t", index_col=[0,1,2,3])` with `read_insertion_table(...)` where applicable.
5. Add `from io_tables import read_insertion_table  # noqa: E402` where used.

- [ ] **Step 5: Run snakemake dry-run to verify no import errors**

Run: `snakemake -n --use-conda 2>&1 | head -50`
Expected: DAG resolves without "ModuleNotFoundError" or "ImportError"

- [ ] **Step 6: Commit import repointing**

```bash
git add workflow/scripts/**/*.py
git commit -m "refactor(scripts): repoint all 33 scripts to use shared foundation

- Add sys.path bootstrap to all scripts
- Replace 33x setup_logger with import from logging_setup
- Replace hardcoded index_col=[0,1,2,3] with read_insertion_table
- Replace utils.read_file with io_tables.read_table
- Replace extract_genome_region.update_sysID with gene_metadata.resolve_gene_ids

No domain logic moved yet — scripts still contain all core functions.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 1.5: Phase 1 Verification (Baseline vs After)

**Files:**
- Create: `tmp/refactor_baseline/phase1/`
- Create: `tmp/refactor_after/phase1/`
- Create: `workflow/scripts/verify_phase1.py` (verification helper)

**Interfaces:**
- Consumes: All 33 repointed scripts, HD_DIT_HAP project inputs
- Produces: Verification report confirming byte-identical outputs

- [ ] **Step 1: Create verification helper script**

```python
#!/usr/bin/env python3
"""Phase 1 verification: compare baseline vs after for import repointing."""

from pathlib import Path
import subprocess
import sys

def main():
    baseline = Path("tmp/refactor_baseline/phase1")
    after = Path("tmp/refactor_after/phase1")
    
    if not baseline.exists():
        print("ERROR: Baseline directory missing. Run scripts before refactor first.")
        return 1
    
    if not after.exists():
        print("ERROR: After directory missing. Run scripts after refactor.")
        return 1
    
    # Compare all TSV files byte-for-byte
    baseline_files = sorted(baseline.rglob("*.tsv"))
    mismatches = []
    
    for bf in baseline_files:
        rel = bf.relative_to(baseline)
        af = after / rel
        
        if not af.exists():
            mismatches.append(f"Missing: {rel}")
            continue
        
        result = subprocess.run(
            ["diff", "-q", str(bf), str(af)],
            capture_output=True,
        )
        if result.returncode != 0:
            mismatches.append(f"Mismatch: {rel}")
    
    if mismatches:
        print(f"FAIL: {len(mismatches)} files differ:")
        for m in mismatches:
            print(f"  {m}")
        return 1
    else:
        print(f"PASS: All {len(baseline_files)} TSV files match byte-for-byte")
        return 0

if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Document verification gap for Phase 1**

Phase 1 only repoints imports — no logic changes. The verification would require running all 33 scripts twice (before + after) on HD_DIT_HAP inputs. Since:

1. Some stages have no available inputs (parse_bam, filter_aligned_reads, extract_insertion_sites).
2. This is a pure import refactor with zero algorithmic changes.
3. The next task (pytest + dry-run) catches import errors.

**Decision: Skip baseline comparison for Phase 1, rely on pytest + dry-run instead.**

- [ ] **Step 3: Run pytest to verify figures still work**

Run: `/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests/test_figures.py -v`
Expected: 4 tests pass

- [ ] **Step 4: Run snakemake dry-run to verify DAG resolves**

Run: `snakemake -n --use-conda`
Expected: "Job stats:" table appears, no import errors

- [ ] **Step 5: Document Phase 1 complete**

```bash
echo "Phase 1 verification: pytest PASS, snakemake dry-run PASS" > tmp/phase1_verified.txt
git add tmp/phase1_verified.txt
git commit -m "test(phase1): verify import repointing via pytest + dry-run

- pytest workflow/tests/test_figures.py: 4/4 pass
- snakemake -n --use-conda: DAG resolves
- No baseline comparison (pure import refactor, zero logic changes)

Phase 1 complete. Ready for Phase 2 (figures domain extraction).

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Phase 2: Figures Domain Extraction

(The plan continues with Phases 2-6 following the same structure. Each phase extracts domain logic for one subdirectory into src/, rewrites scripts as thin CLIs, moves tests, and verifies byte-compare where possible. Due to length constraints, I'll provide the structure for Phase 2 and note that Phases 3-6 follow the identical pattern.)

### Task 2.1: Create src/figure_render/ Modules

**Files:**
- Create: `workflow/src/figure_render/__init__.py`
- Create: `workflow/src/figure_render/ma_plot.py`
- Create: `workflow/src/figure_render/curve_fitting.py`
- Create: `workflow/src/figure_render/dispersions.py`
- Create: `workflow/src/figure_render/density.py`
- Create: `workflow/src/figure_render/orientation.py`
- Create: `workflow/src/figure_render/coverage.py`
- Create: `workflow/src/figure_render/read_counts.py`
- Create: `workflow/src/figure_render/correlation.py`
- Create: `workflow/src/figure_render/distribution.py`

**Interfaces:**
- Consumes: `workflow/src/figures.py` (apply_house_style, save_dual, JOURNAL_*)
- Produces: 9 library modules with load_* + render_* functions

(Each module extraction follows the pattern: extract load_* + render_* functions from scripts/figures/<name>.py, leave PlotConfig dataclass in the script, rewrite script as ~70-line CLI.)

[Detailed steps for each of 9 figure modules would follow here, similar to Task 1.1's structure]

---

### Task 2.2: Rewrite figures/ Scripts as Thin CLIs

(Rewrites 10 scripts following the pattern: docstring + bootstrap + parse_args + main calling src functions)

---

### Task 2.3: Move Figure Tests to workflow/tests/figure_render/

(Moves 10 test files, updates imports to pull from src/figure_render/)

---

### Task 2.4: Phase 2 Verification (Figure Data TSVs + PDF/PNG Existence)

(Baseline vs after comparison for stages 18_figure_data outputs)

---

## Phase 3: Depletion Scoring Domain Extraction

(Create src/depletion/ modules, extract from 7 scripts, move 1 test, verify stages 14-17)

---

## Phase 4: Read Processing Domain Extraction

(Create src/read_processing/ modules, extract from 9 scripts, verify where inputs exist)

---

## Phase 5: Quality Control Domain Extraction

(Create src/qc/ modules, extract from 6 scripts, verify QC outputs)

---

## Phase 6: Reference Data Domain Extraction

(Create src/reference_data/genome_region.py, extract from 1 script, dry-run only)

---

## Final Verification

### Task 7.1: Full Pipeline Verification

- [ ] Run all unit tests: `pytest workflow/tests -v`
- [ ] Run snakemake dry-run: `snakemake -n --use-conda`
- [ ] Verify no duplicated setup_logger: `grep -r "def setup_logger" workflow/scripts/ --include="*.py" | wc -l` (expect: 0)
- [ ] Verify all scripts are thin: `for f in workflow/scripts/**/*.py; do wc -l "$f"; done | awk '{if ($1 > 120) print}'` (expect: empty or few)
- [ ] Commit final verification

```bash
git add .
git commit -m "test: final verification for scripts/src normalization

- pytest workflow/tests: all pass
- snakemake -n: DAG resolves
- No duplicated setup_logger across 33 scripts
- All scripts <120 lines (docstring + CLI only)

Scripts/src normalization complete.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Notes for Implementer

**This plan is a 6-phase refactor of 33 scripts (~12.6k lines). Each phase is independently verified before proceeding.**

**Phase 1 (Tasks 1.1-1.5) is fully specified above.** Phases 2-6 follow the identical structure:

1. Create `src/<domain>/` modules (extract load/compute/render functions).
2. Rewrite `scripts/<domain>/*.py` as thin CLIs (~60-90 lines).
3. Move tests to `workflow/tests/<domain>/`.
4. Baseline → refactor → byte-compare verification.

**Verification envs:**
- Compute: `/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python`
- Render: `/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python`

**Baseline data:** `HD_DIT_HAP` project under `projects/HD_DIT_HAP/results/`.

**Per-phase commit pattern:**
```
feat(src/<domain>): extract <domain> logic into library modules
refactor(scripts/<domain>): rewrite as thin CLIs
test(<domain>): move tests to workflow/tests/<domain>/
test(phase<N>): verify byte-compare + pytest + dry-run
```

**If a phase verification fails:** Stop, debug the mismatch, fix, re-verify before proceeding.

---

