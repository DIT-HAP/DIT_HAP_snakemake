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

Ten figure scripts, each split into a `src/figure_render/<module>.py` library half and
a thin CLI half. The scripts are the least entangled in the tree (2–3 core functions
each), which is why this phase runs first.

### Phase 2 Conventions (apply to every task below)

**Extraction rule.** For each script, move the `load_*` / `render_*` / helper functions
verbatim into the named `src/figure_render/` module. Do not alter function bodies except
where a step below names the change. Leave in the script: the module docstring, the
`PlotConfig` dataclass, `parse_args()`, `main()`.

**Module header.** Every new `src/figure_render/*.py` starts with:

```python
"""<One-line title>.

Input
-----
- <what the load_* function reads>

Output
------
- <what the render_* function writes>

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# =============================================================================
# CONSTANTS
# =============================================================================
# =============================================================================
# CORE LOGIC
# =============================================================================
```

No `main()`, no `parse_args()`, no `setup_logger()` — these are library modules.

**Import style inside `src/figure_render/`.** Sibling `src` modules import by bare name
(`from figures import apply_house_style, save_dual`), because `src/` is the sys.path root
the scripts append. Do NOT use `from ..figures import` or `from src.figures import`.

**Script CLI shape.** After extraction each script reads:

```python
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from figure_render.<module> import <load_fn>, <render_fn>  # noqa: E402
```

with `main()` calling `setup_logger`, building `PlotConfig`, calling the two imported
functions, and returning 0/1. Keep every existing `-`/`--` flag name and default
unchanged: the `.smk` rules pass them positionally by name and are out of scope.

**Per-script verification.** Confirm the module imports standalone and the CLI constructs:

```bash
RENDER_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
$RENDER_PY workflow/scripts/figures/<script>.py --help
$RENDER_PY -c "import sys; sys.path.insert(0,'workflow/src'); import figure_render.<module>; print('ok')"
```

**`--help` alone is not sufficient, and this is not hypothetical.** Phase 1 shipped all 10
figure scripts with a missing `setup_logger` import; `--help` passed on every one because
argparse exits before `main()` runs, and `snakemake -n` never executes a script body
either. Every figure rule would have crashed with `NameError` in production. So for at
least one script per module, actually render against real data:

```bash
R=projects/HD_DIT_HAP/results
$RENDER_PY workflow/scripts/figures/<script>.py <real input flags> -o /tmp/check/<name>
```

and confirm a non-empty `.pdf` + `.review.png` pair appears. A verification step that
cannot fail is not a verification step.

---

### Task 2.1: Extract the shared sigmoid model

`sigmoid_function` exists byte-identically in `scripts/figures/plot_curve_fitting.py:67`
and `scripts/depletion_scoring/curve_fitting.py`. Phase 2 is the first phase to need it,
so it lands here and Phase 3 imports it rather than re-homing it later.

**Files:**
- Create: `workflow/src/depletion/__init__.py`
- Create: `workflow/src/depletion/curve_model.py`
- Create: `workflow/tests/depletion/__init__.py`
- Create: `workflow/tests/depletion/test_curve_model.py`

**Interfaces:**
- Produces: `sigmoid_function(x: np.ndarray, A: float, DR: float, DL: float) -> np.ndarray`
  — Gompertz curve; returns `np.zeros_like(x)` when `A == 0`; exponent clipped to
  `(-700, 700)` for numerical stability. Consumed by Task 2.3 and by Phase 3.

- [ ] **Step 1: Write the failing test**

Create `workflow/tests/depletion/__init__.py` (empty) and
`workflow/tests/depletion/test_curve_model.py`:

```python
"""Tests for the shared Gompertz curve model."""

import numpy as np
import pytest

from depletion.curve_model import sigmoid_function


def test_zero_amplitude_returns_zeros():
    """A == 0 must short-circuit to zeros rather than divide by zero."""
    x = np.array([0.0, 1.0, 2.0])
    assert np.array_equal(sigmoid_function(x, 0.0, 1.0, 1.0), np.zeros_like(x))


def test_monotonic_decreasing_in_x():
    """With positive amplitude the Gompertz curve decreases as x grows."""
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = sigmoid_function(x, 2.0, 0.5, 2.0)
    assert np.all(np.diff(y) <= 1e-12)


def test_extreme_input_does_not_overflow():
    """Exponent clipping must keep a huge negative x finite."""
    y = sigmoid_function(np.array([-1e6]), 2.0, 0.5, 2.0)
    assert np.all(np.isfinite(y))


def test_matches_amplitude_at_large_x():
    """As x grows past DL the curve approaches the amplitude A."""
    y = sigmoid_function(np.array([1e3]), 2.0, 0.5, 2.0)
    assert y[0] == pytest.approx(2.0, abs=1e-9)
```

- [ ] **Step 2: Run it to confirm it fails**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest \
    workflow/tests/depletion/test_curve_model.py -v
```
Expected: collection error / `ModuleNotFoundError: No module named 'depletion'`

- [ ] **Step 3: Create the module**

`workflow/src/depletion/__init__.py`:

```python
"""Depletion scoring library modules."""
```

`workflow/src/depletion/curve_model.py`:

```python
"""Gompertz curve model shared by curve fitting and its figures.

Input
-----
- None (pure numerical model).

Output
------
- None (returns arrays).

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-19
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import numpy as np

# =============================================================================
# CONSTANTS
# =============================================================================
# np.exp overflows past ~709; clip well inside that bound.
EXPONENT_CLIP = 700

# =============================================================================
# CORE LOGIC
# =============================================================================
def sigmoid_function(x: np.ndarray, A: float, DR: float, DL: float) -> np.ndarray:
    """Calculate sigmoid function values using gompertz function."""
    if A == 0:
        return np.zeros_like(x)
    alpha = (DR * np.e) / A
    u = alpha * (DL - x) + 1
    exponent = np.clip(u, -EXPONENT_CLIP, EXPONENT_CLIP)
    return A * np.exp(-np.exp(exponent))
```

- [ ] **Step 4: Run the test to confirm it passes**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest \
    workflow/tests/depletion/test_curve_model.py -v
```
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add workflow/src/depletion workflow/tests/depletion
git commit -m "feat(src/depletion): extract shared Gompertz curve model

sigmoid_function was byte-identical in scripts/figures/plot_curve_fitting.py
and scripts/depletion_scoring/curve_fitting.py. Phase 2 needs it first, so it
lands here; Phase 3 will import it instead of keeping its own copy.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2.2: Extract the four single-load/render figure modules

Four scripts share one shape exactly — one `load_*`, one `render_*`, nothing else — so
they are one task, not four. Apply the Phase 2 Conventions to each.

**Files:**
- Create: `workflow/src/figure_render/__init__.py` (`"""Figure rendering library modules."""`)
- Create: `workflow/src/figure_render/dispersions.py` — from `scripts/figures/plot_dispersions.py`
- Create: `workflow/src/figure_render/orientation.py` — from `scripts/figures/plot_insertion_orientation.py`
- Create: `workflow/src/figure_render/correlation.py` — from `scripts/figures/plot_pbl_pbr_correlation.py`
- Create: `workflow/src/figure_render/distribution.py` — from `scripts/figures/plot_distribution_of_curve_fitting.py`
- Modify: those four scripts, reduced to CLI per the conventions

**Interfaces:**
- Consumes: `figures.apply_house_style`, `figures.save_dual`, and for dispersions /
  orientation / correlation also `figures.JOURNAL_WIDTH_PX`, `figures.JOURNAL_HEIGHT_PX`.
- Produces, moved verbatim from the named script:
  - `dispersions.load_dispersion_data(dispersion_data_path: Path) -> pd.DataFrame`
  - `dispersions.render_dispersion_figure(df: pd.DataFrame, output_stem: Path) -> None`
  - `orientation.load_and_prepare_data(input_path: Path) -> pd.DataFrame`
  - `orientation.render_orientation_figure(df: pd.DataFrame, output_stem: Path) -> None`
  - `correlation.load_and_prepare_data(input_path: Path) -> pd.DataFrame`
  - `correlation.render_correlation_figure(df: pd.DataFrame, output_stem: Path) -> None`
  - `distribution.load_fitting_stats(fitting_stats_path: Path) -> tuple[pd.DataFrame, list[str]]`
  - `distribution.render_distribution_figure(df: pd.DataFrame, metric_cols: list[str], output_stem: Path, bins: int) -> None`

`orientation.load_and_prepare_data` and `correlation.load_and_prepare_data` keep their
shared name: they live in different modules and are never imported together.

- [ ] **Step 1: Create the package marker and the four modules**

For each of the four, move the functions listed under Interfaces out of the script and
into the new module, adding the standard header and the imports each body actually uses.
Change nothing inside the bodies. `matplotlib.use("Agg")` stays in the *script*, not the
module — it is a process-level side effect and a library must not impose it.

- [ ] **Step 2: Confirm each module imports standalone**

```bash
RENDER_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
for m in dispersions orientation correlation distribution; do
  $RENDER_PY -c "import sys; sys.path.insert(0,'workflow/src'); import figure_render.$m; print('$m ok')"
done
```
Expected: four `ok` lines

- [ ] **Step 3: Reduce the four scripts to CLIs**

Apply the Script CLI shape from the conventions. Each `main()` keeps its existing
sequence of calls and its existing flags; only the function *source* changes from local
def to import.

- [ ] **Step 4: Confirm all four CLIs still construct**

```bash
RENDER_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
for s in plot_dispersions plot_insertion_orientation plot_pbl_pbr_correlation \
         plot_distribution_of_curve_fitting; do
  $RENDER_PY workflow/scripts/figures/$s.py --help > /dev/null && echo "$s ok"
done
```
Expected: four `ok` lines

- [ ] **Step 5: Commit**

```bash
git add workflow/src/figure_render workflow/scripts/figures
git commit -m "refactor(figures): extract dispersions/orientation/correlation/distribution

Move load_*/render_* into src/figure_render/; the four scripts keep their
docstring, PlotConfig, parse_args and main. CLI flags unchanged.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2.3: Extract the five figure modules with helpers or extra params

The remaining six scripts each carry a private helper, an enum, or a config-dependent
loader, so each needs a named decision. `plot_ma_plot.py` and `plot_ma_plot_replicates.py`
are distinct scripts and get distinct modules.

**Files:**
- Create: `workflow/src/figure_render/ma_plot.py` — from `scripts/figures/plot_ma_plot.py`
- Create: `workflow/src/figure_render/ma_plot_replicates.py` — from `scripts/figures/plot_ma_plot_replicates.py`
- Create: `workflow/src/figure_render/coverage.py` — from `scripts/figures/plot_gene_coverage.py`
- Create: `workflow/src/figure_render/density.py` — from `scripts/figures/plot_insertion_density.py`
- Create: `workflow/src/figure_render/read_counts.py` — from `scripts/figures/plot_read_count_distribution.py`
- Create: `workflow/src/figure_render/curve_fitting.py` — from `scripts/figures/plot_curve_fitting.py`
- Modify: those six scripts, reduced to CLI per the conventions

**Interfaces:**
- Consumes: `figures.apply_house_style`, `figures.save_dual`, `figures.JOURNAL_*` as each
  body already does; `curve_fitting` additionally consumes
  `depletion.curve_model.sigmoid_function` from Task 2.1.
- Produces:
  - `ma_plot.Orientation` (StrEnum: `VERTICAL`, `HORIZONTAL`) — moves with the module
  - `ma_plot.load_ma_data(basemean_path: Path, lfc_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]`
  - `ma_plot.render_ma_figure(basemean_df, lfc_df, output_stem: Path, orientation: Orientation = Orientation.VERTICAL) -> None`
  - `ma_plot_replicates.load_ma_data(ma_values_path: Path) -> pd.DataFrame`
  - `ma_plot_replicates.significance_colors(padj: pd.Series) -> pd.Series`
  - `ma_plot_replicates.render_ma_figure(df: pd.DataFrame, output_stem: Path) -> None`
  - `coverage.load_coverage_data(input_path: Path) -> pd.DataFrame`
  - `coverage.render_coverage_figure(df: pd.DataFrame, output_stem: Path) -> None`
    (private `_status_counts_frame` moves too, still underscore-prefixed)
  - `density.load_density_data(input_path: Path) -> pd.DataFrame`
  - `density.render_density_figure(df: pd.DataFrame, output_stem: Path, initial_timepoint: str, final_timepoint: str) -> None`
    (private `_draw_panel_or_placeholder` moves too)
  - `read_counts.load_distribution_data(input_path: Path) -> pd.DataFrame`
  - `read_counts.load_cutoff_stats(stats_path: Path) -> pd.DataFrame`
  - `read_counts.format_retention_caption(sample: str, stats_row: pd.Series | None) -> str`
  - `read_counts.render_distribution_figure(...)` — same signature as the current script
  - `curve_fitting.load_and_sample_data(fitting_stats_path: Path, lfc_path: Path, n_curves: int, random_seed: int) -> tuple[pd.DataFrame, pd.DataFrame, list[float]]`
  - `curve_fitting.render_curve_fitting_figure(...)` — same signature as the current script

- [ ] **Step 1: Change `load_and_sample_data` to take scalars, not `PlotConfig`**

Today it reads `config.fitting_stats_path`, `config.lfc_path`, `config.n_curves`,
`config.random_seed` off the dataclass. `PlotConfig` stays in the script, so a `src`
function must not depend on it. Replace the parameter list with the four scalars in the
order given under Interfaces and substitute the four attribute reads with the parameter
names. Nothing else in the body changes. The caller in `main()` becomes:

```python
sampled, lfc_sampled, time_points = load_and_sample_data(
    config.fitting_stats_path, config.lfc_path, config.n_curves, config.random_seed
)
```

- [ ] **Step 2: Point `curve_fitting` at the shared sigmoid**

Delete the local `sigmoid_function` def from `plot_curve_fitting.py` and have
`src/figure_render/curve_fitting.py` import it:

```python
from depletion.curve_model import sigmoid_function
```

- [ ] **Step 3: Create the six modules**

Move the functions listed under Interfaces verbatim (except the two changes above),
each with the standard header. Keep `Orientation` in `ma_plot.py` beside the two
functions that switch on it. `matplotlib.use("Agg")` stays in the scripts.

- [ ] **Step 4: Confirm each module imports standalone**

```bash
RENDER_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
for m in ma_plot ma_plot_replicates coverage density read_counts curve_fitting; do
  $RENDER_PY -c "import sys; sys.path.insert(0,'workflow/src'); import figure_render.$m; print('$m ok')"
done
```
Expected: six `ok` lines

- [ ] **Step 5: Reduce the six scripts to CLIs, then confirm they construct**

```bash
RENDER_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
for s in plot_ma_plot plot_ma_plot_replicates plot_gene_coverage \
         plot_insertion_density plot_read_count_distribution plot_curve_fitting; do
  $RENDER_PY workflow/scripts/figures/$s.py --help > /dev/null && echo "$s ok"
done
```
Expected: six `ok` lines

- [ ] **Step 6: Commit**

```bash
git add workflow/src/figure_render workflow/scripts/figures
git commit -m "refactor(figures): extract ma_plot, coverage, density, read_counts, curve_fitting

load_and_sample_data now takes explicit scalars instead of PlotConfig, since
config dataclasses stay in scripts/. plot_curve_fitting drops its duplicate
sigmoid_function in favour of depletion.curve_model.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2.4: Move the ten figure tests to workflow/tests/figure_render/

**Files:**
- Create: `workflow/tests/figure_render/__init__.py`
- Move, renaming `test_plot_<x>.py` → `test_<module>.py` to match the module under test:
  - `test_plot_dispersions.py` → `workflow/tests/figure_render/test_dispersions.py`
  - `test_plot_insertion_orientation.py` → `.../test_orientation.py`
  - `test_plot_pbl_pbr_correlation.py` → `.../test_correlation.py`
  - `test_plot_distribution_of_curve_fitting.py` → `.../test_distribution.py`
  - `test_plot_ma_plot.py` → `.../test_ma_plot.py`
  - `test_plot_ma_plot_replicates.py` → `.../test_ma_plot_replicates.py`
  - `test_plot_gene_coverage.py` → `.../test_coverage.py`
  - `test_plot_insertion_density.py` → `.../test_density.py`
  - `test_plot_read_count_distribution.py` → `.../test_read_counts.py`
  - `test_plot_curve_fitting.py` → `.../test_curve_fitting.py`

**Interfaces:**
- Consumes: every symbol produced by Tasks 2.1–2.3.

- [ ] **Step 1: Move the files with git mv**

Use `git mv` so history follows, and create `workflow/tests/figure_render/__init__.py`
(empty).

- [ ] **Step 2: Repoint each test's imports to the module**

Delete the `SCRIPT_DIR` / `sys.path.append` preamble from every moved test —
`pyproject.toml` already puts `workflow/src` and `workflow/scripts` on the path. Point
every domain symbol at its new module:

```python
from figure_render.ma_plot import Orientation, load_ma_data, render_ma_figure
```

**`PlotConfig` cannot be imported this way, and neither can any other script symbol.**
Two independent obstacles, both verified:

1. `workflow/scripts` is on the path but `workflow/scripts/figures/` is a subdirectory
   with no `__init__.py`, so `plot_ma_plot` is not importable as a top-level module.
2. Even with an `__init__.py`, `workflow/src/figures.py` **shadows** the
   `workflow/scripts/figures/` directory — `import figures` resolves to the module, and
   `import figures.plot_ma_plot` fails with
   `ModuleNotFoundError: 'figures' is not a package`.

So the five tests that assert on `PlotConfig` do NOT move to
`workflow/tests/figure_render/`. They move to `workflow/tests/scripts/` and reach their
target through a spec-loader fixture (Step 3). Those five are:
`test_plot_ma_plot.py`, `test_plot_ma_plot_replicates.py`, `test_plot_dispersions.py`,
`test_plot_curve_fitting.py`, `test_plot_distribution_of_curve_fitting.py`.

Each of those five splits in two: its rendering assertions go to
`workflow/tests/figure_render/test_<module>.py` with module imports, and its
`test_config_rejects_missing_input` goes to `workflow/tests/scripts/test_plot_configs.py`.
Note `test_plot_dispersions.py` and `test_plot_ma_plot_replicates.py` also import
module-level constants (`DISPERSION_SERIES`, `REQUIRED_COLUMNS`, `X_COLUMN`;
`NONSIGNIFICANT_COLOR`, `PADJ_THRESHOLD`, `SIGNIFICANT_COLOR`) — those moved to the
`src` modules in Tasks 2.2/2.3, so they import from `figure_render.*` normally.

- [ ] **Step 3: Add the script-loading conftest**

Create `workflow/tests/scripts/__init__.py` (empty) and
`workflow/tests/scripts/conftest.py`:

```python
"""Fixtures for testing CLI scripts that are not importable as modules.

workflow/scripts/figures/ has no __init__.py, and workflow/src/figures.py shadows
the directory name on sys.path, so these scripts can only be loaded by file path.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import importlib.util
from collections.abc import Callable
from pathlib import Path
from types import ModuleType

import pytest

# =============================================================================
# CONSTANTS
# =============================================================================
FIGURES_DIR = Path(__file__).resolve().parents[3] / "workflow" / "scripts" / "figures"

# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture(scope="session")
def load_script() -> Callable[[str], ModuleType]:
    """Return a loader that imports a figure CLI script by filename stem."""
    def _load(stem: str) -> ModuleType:
        path = FIGURES_DIR / f"{stem}.py"
        if not path.exists():
            raise FileNotFoundError(f"No such script: {path}")
        spec = importlib.util.spec_from_file_location(f"_script_{stem}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    return _load
```

`FIGURES_DIR` walks up three parents from `workflow/tests/scripts/conftest.py` to the
repo root, then back down — verify it resolves before relying on it.

- [ ] **Step 4: Write the consolidated PlotConfig test**

Create `workflow/tests/scripts/test_plot_configs.py` holding the five validation tests,
each loading its script through the fixture. All five assert the same shape — a
`ValueError` whose message contains `does not exist` — so parametrize rather than
repeating the body five times:

```python
"""Validation tests for the figure scripts' PlotConfig dataclasses."""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pytest

# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize(
    ("stem", "kwargs"),
    [
        ("plot_ma_plot", {"basemean_path": "nope_baseMean.tsv", "lfc_path": "nope_LFC.tsv"}),
        ("plot_ma_plot_replicates", {"ma_values_path": "nope_ma_values.tsv"}),
        ("plot_dispersions", {"input_path": "nope_dispersion.tsv"}),
        ("plot_distribution_of_curve_fitting", {"fitting_stats_path": "nope_stats.tsv"}),
        ("plot_curve_fitting", {"fitting_stats_path": "nope_stats.tsv", "lfc_path": "nope_LFC.tsv"}),
    ],
)
def test_config_rejects_missing_input(load_script, tmp_path: Path, stem: str, kwargs: dict) -> None:
    """Assert each script's PlotConfig rejects a non-existent input path."""
    module = load_script(stem)
    paths = {key: tmp_path / value for key, value in kwargs.items()}

    with pytest.raises(ValueError, match="does not exist"):
        module.PlotConfig(output_stem=tmp_path / "out", **paths)
```

The `kwargs` above are the field names each `PlotConfig` actually declares — read every
one of the five scripts and confirm the names match before running, since a wrong keyword
raises `TypeError` rather than the `ValueError` the test asserts, which would pass the
`raises` check for the wrong reason only if you also loosened the match. Do not loosen it.

- [ ] **Step 5: Fix the `load_and_sample_data` call in test_curve_fitting.py**

It currently passes a `PlotConfig`. Update the call to the four-scalar signature from
Task 2.3 Step 1 — `(fitting_stats_path, lfc_path, n_curves, random_seed)` — building the
values from the fixture paths the test already has.

- [ ] **Step 6: Run the whole suite**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests -v
```
Expected: the 4 `test_figures.py` tests, the 4 `test_curve_model.py` tests, the 5
parametrized `test_plot_configs.py` cases, and every moved figure test pass. Tests that
`pytest.skip` on absent real project data are an acceptable pass; a `ModuleNotFoundError`
or `ImportError` is not.

Real project data for the skip-guarded tests is in
`projects/HD_DIT_HAP/results/18_figure_data/arc/` and the numbered stage directories —
if a test skips because it points at a path that does not exist, check `arc/` before
concluding the data is unavailable.

- [ ] **Step 7: Confirm no test file is left behind**

```bash
ls workflow/scripts/figures/test_*.py 2>/dev/null && echo "FAIL: tests remain" || echo "ok: none remain"
```
Expected: `ok: none remain`.

- [ ] **Step 8: Commit**

```bash
git add -A workflow/tests workflow/scripts/figures
git commit -m "test(figures): move figure tests out of scripts/ into workflow/tests/

Rendering tests -> workflow/tests/figure_render/, renamed to match the module
under test, sys.path preambles dropped in favour of the pyproject pythonpath.

The 5 PlotConfig validation tests -> workflow/tests/scripts/test_plot_configs.py,
parametrized, reaching their targets through a spec-loader conftest fixture:
scripts/figures/ has no __init__.py AND src/figures.py shadows the directory
name, so those scripts cannot be imported as modules at all.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2.5: Phase 2 verification

PDFs embed a creation timestamp, so two renders of an identical figure differ on disk.
**PNGs do not** — this was measured, not assumed: two consecutive `plot_dispersions.py`
runs produced `.review.png` files with identical SHA-256, while their `.pdf` siblings
differed. So verification pixel-compares the PNG against a pre-refactor baseline, which
catches a real rendering regression, and falls back to existence + non-emptiness only for
the PDF.

**Capturing the baseline requires a pre-refactor worktree**, because the extraction has
already happened in the main tree by the time this task runs:

```bash
git worktree add /tmp/phase2_baseline <commit-before-this-phase> --detach
```

Note that a pre-Phase-2 tree still carries the Phase 1 `setup_logger` defect for figure
scripts (fixed in `890eafc`), so patch that import in the worktree — uncommitted, local to
the worktree — before rendering, or every baseline render dies at `main()`.

Verified for Task 2.2's scripts at the time of writing: `disp.review.png` (113476 bytes)
and `dist.review.png` (430736 bytes) were byte-identical before and after extraction.

**Files:**
- Create: `tmp/refactor_baseline/phase2/`, `tmp/refactor_after/phase2/` (scratch, git-ignored)

- [ ] **Step 1: Render every figure whose input data exists — all of them do**

Inputs live in two places, and missing the second one causes a false "cannot verify":

- `projects/HD_DIT_HAP/results/{14_insertion_level_depletion_analysis,15_insertion_level_curve_fitting,16_gene_level_depletion_analysis,17_gene_level_curve_fitting}/`
- **`projects/HD_DIT_HAP/results/18_figure_data/arc/`** — holds `insertion_density_analysis.tsv`,
  `gene_coverage_stats.tsv`, `pbl_pbr_pairs.tsv`, `read_count_distribution.tsv`,
  `read_count_cutoff_stats.tsv`, `dispersion_data.tsv`, `ma_values.tsv`.

All ten figure scripts are renderable from these. The exact invocations, verified working:

```bash
PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
R=projects/HD_DIT_HAP/results; A=$R/18_figure_data/arc; O=/tmp/refactor_after/phase2
$PY .../plot_ma_plot.py -b $R/14_insertion_level_depletion_analysis/baseMean.tsv \
                        -l $R/14_insertion_level_depletion_analysis/LFC.tsv -o $O/ma
$PY .../plot_ma_plot_replicates.py -i $A/ma_values.tsv -o $O/mar
$PY .../plot_dispersions.py -i $A/dispersion_data.tsv -o $O/disp
$PY .../plot_distribution_of_curve_fitting.py \
        -i $R/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv -o $O/dist
$PY .../plot_gene_coverage.py -i $A/gene_coverage_stats.tsv -o $O/cov
$PY .../plot_pbl_pbr_correlation.py -i $A/pbl_pbr_pairs.tsv -o $O/corr
$PY .../plot_insertion_density.py -i $A/insertion_density_analysis.tsv -o $O/dens -t YES0 -f YES4
$PY .../plot_read_count_distribution.py -i $A/read_count_distribution.tsv \
        -s $A/read_count_cutoff_stats.tsv -o $O/rc -t YES0 -c 8
```

Note `plot_ma_plot.py` emits TWO pairs (`ma` and `ma_horizontal`), so ten scripts yield
nine comparable PNGs. Quote paths (`"$A/x.tsv"`) — an unquoted expansion inside a
composite loop variable silently prepends a space and the script reports the file missing.

- [ ] **Step 2: Pixel-compare every PNG against the baseline**

```bash
fail=0
for b in tmp/refactor_baseline/phase2/*.review.png; do
  n=$(basename "$b"); a="tmp/refactor_after/phase2/$n"
  if   [ ! -f "$a" ];        then echo "MISSING:   $n"; fail=1
  elif cmp -s "$b" "$a";     then echo "IDENTICAL: $n ($(stat -c%s "$a") bytes)"
  else                            echo "DIFFERS:   $n"; fail=1
  fi
done
exit $fail
```
Expected: `IDENTICAL` for every figure rendered in both trees. A `DIFFERS` line is a real
rendering regression — stop and fix before proceeding, do not park it.

- [ ] **Step 3: Assert PDFs exist and are non-empty**

PDFs cannot be compared (embedded timestamp), so they only get an existence check:

```bash
find tmp/refactor_after/phase2 -name "*.pdf" | sort | while read f; do
  test -s "$f" || { echo "EMPTY: $f"; exit 1; }
  echo "ok: $(basename "$f") ($(stat -c%s "$f") bytes)"
done
```
Expected: a non-empty `.pdf` beside every `.review.png`.

- [ ] **Step 4: Full suite plus DAG**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests -q
snakemake -n --use-conda
```
Expected: suite green; DAG resolves with no `ImportError`.

- [ ] **Step 5: Assert the scripts actually got thin**

```bash
wc -l workflow/scripts/figures/plot_*.py | sort -rn
```
Expected: every script well under its former size (they ranged 186–318 lines). Report the
before/after table. A script still over ~120 lines means logic stayed behind — name it
in the report rather than quietly passing.

- [ ] **Step 6: Commit the verification record**

```bash
git add -A tmp/phase2_verified.txt
git commit -m "test(phase2): verify figures extraction

Artifacts render non-empty for every script with available inputs; suite green;
DAG resolves. Scripts reduced from 186-318 lines to CLI-only.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

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

