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

Seven scripts, 2,582 lines. Unlike Phase 2 these produce **TSV** outputs, so verification
byte-compares the data itself — the strongest check available in this refactor.

### Phase 3 Conventions (apply to every task below)

Identical to Phase 2's conventions, with these deltas:

**Extraction rule.** Move function bodies VERBATIM. Do not reformat, reorder, or fix
anything — including code that looks wrong. This is a scientific pipeline; silent
numerical drift corrupts published results.

**Module header.** Same 7-section library layout as Phase 2 (title/Input/Output/metadata
docstring → IMPORTS → CONSTANTS → CORE LOGIC). No `main()`, no `parse_args()`, no
`setup_logger()`.

**Import style.** Bare sibling names from `src/`: `from logging_setup import setup_logger`,
`from io_tables import read_insertion_table`, `from depletion.curve_model import
sigmoid_function`. Never `from ..x` or `from src.x`.

**What stays in the script.** The module docstring, the **config** dataclass, `parse_args()`,
`main()`. Config dataclasses are the ones that validate CLI paths in `__post_init__`:
`WeightsConfig`, `CurveFittingConfig`, `ImputationConfig`, `AnalysisConfig`,
`InputOutputConfig`. They stay.

**What moves.** Domain logic plus the **result/data** dataclasses that functions return or
pass between themselves: `WeightInputs`, `FittingResult`, `SummaryStatistics`,
`ControlSelectionResult`, `ImputationResult`, `AnalysisResult`. Enums are domain vocabulary
and move too: `Scheme`, `AnnotCol`, `StatsCol`.

**Config-coupled functions must take scalars.** A `src` module must never import a config
dataclass from a script. Two functions currently violate this and are named in their tasks.

**Verification is per-script and mandatory.** `--help` proves nothing (Phase 1 shipped 10
broken scripts that all passed `--help`). Every script must be **executed** against real
data and its TSV output byte-compared to a pre-refactor baseline. Capture the baseline
FIRST, from a worktree at the pre-Phase-3 commit:

```bash
git worktree add /tmp/phase3_baseline <pre-phase-3-sha> --detach
```

Real data: `projects/HD_DIT_HAP/results/{13_filtered,14_insertion_level_depletion_analysis,15_insertion_level_curve_fitting,16_gene_level_depletion_analysis,17_gene_level_curve_fitting}/`.
Compute env (pandas 3, pydeseq2, scipy, joblib):
`/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python`.
Read the invoking rule in `workflow/rules/depletion_scoring.smk` for each script's exact
flags rather than guessing them.

**Known dead code — leave it alone.** `curve_fitting.py` holds `create_fitted_plot` and
`generate_fitting_plots`; `insertion_level_depletion_analysis_no_replicates.py` holds
`generate_MA_plots`. All three are matplotlib figure code, called nowhere in the repo, and
superseded by Phase 2's `plot_curve_fitting.py` / `plot_ma_plot.py`. Per user decision they
stay in the scripts for now — do NOT move them into `src/`, do NOT delete them, do NOT
"fix" them. They are on the final cleanup list.

---

### Task 3.1: Extract curve fitting and retire the sigmoid duplicate

`curve_fitting.py` is the largest script (623 lines) and holds the last duplicate of
`sigmoid_function`, which Task 2.1 already canonicalised in `src/depletion/curve_model.py`.

**Files:**
- Create: `workflow/src/depletion/curve_fitting.py`
- Modify: `workflow/scripts/depletion_scoring/curve_fitting.py`

**Interfaces:**
- Consumes: `depletion.curve_model.sigmoid_function`
- Produces (moved verbatim unless noted):
  - `FittingResult`, `SummaryStatistics` (result dataclasses)
  - `sigmoid_derivative(x, A, DR, DL) -> np.ndarray`
  - `time_at_p_effect(p, A, DR, DL) -> float`
  - `objective_function(params, x, y, ...) -> float`
  - `constraint_function1(params, t_last) -> float`
  - `constraint_function2(params) -> float`
  - `fit_single_curve(x_values, y_values, ...) -> FittingResult`
  - `fit_and_augment(x_values, y_data, ...)`
  - `process_depletion_data(input_file: Path, time_points: list[float], ...)`
  - `generate_summary_statistics(results_df) -> SummaryStatistics`
  - `display_summary_table(stats: SummaryStatistics) -> None`

- [ ] **Step 1: Delete the duplicate sigmoid, import the canonical one**

Remove `def sigmoid_function` from the script (it is at `curve_fitting.py:174`). The new
`src/depletion/curve_fitting.py` imports it:

```python
from depletion.curve_model import sigmoid_function
```

Confirm byte-identity before deleting, so you know the swap is safe:

```bash
diff <(sed -n '/^def sigmoid_function/,/^$/p' workflow/scripts/depletion_scoring/curve_fitting.py) \
     <(sed -n '/^def sigmoid_function/,/^$/p' workflow/src/depletion/curve_model.py)
```
Expected: differs only in the docstring line, if at all. If the bodies differ, STOP and
report — that would mean the two copies had drifted and Task 2.1's premise was wrong.

- [ ] **Step 2: Move the listed functions and dataclasses**

Verbatim, into `src/depletion/curve_fitting.py`. Leave `CurveFittingConfig`, `parse_args`,
`main`, and the two dead plotting functions in the script. Because the dead plotting
functions stay, `import matplotlib.pyplot as plt` stays in the script too.

- [ ] **Step 3: Confirm the module imports and the script runs**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python
$PY -c "import sys; sys.path.insert(0,'workflow/src'); import depletion.curve_fitting; print('ok')"
$PY workflow/scripts/depletion_scoring/curve_fitting.py --help > /dev/null && echo "help ok"
```

- [ ] **Step 4: Execute against real data and byte-compare**

Read the `insertion_level_curve_fitting` rule in `workflow/rules/depletion_scoring.smk`
for the exact flags, then run the script in both the baseline worktree and the current
tree, writing to separate directories, and `diff` every TSV produced. Report the exact
command used and the diff result per file.

- [ ] **Step 5: Commit**

```bash
git add workflow/src/depletion/curve_fitting.py workflow/scripts/depletion_scoring/curve_fitting.py
git commit -m "refactor(depletion): extract curve fitting, retire sigmoid duplicate

The last duplicate of sigmoid_function is gone; curve fitting now imports the
canonical one from depletion.curve_model (Task 2.1). Dead plotting functions
left in the script per user decision, pending final cleanup.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3.2: Extract weights, gene level, and control insertions

Three pandas-only scripts with no heavy dependencies — batched because they share a shape.

**Files:**
- Create: `workflow/src/depletion/weights.py` — from `compute_insertion_weights.py`
- Create: `workflow/src/depletion/gene_level.py` — from `gene_level_depletion_analysis.py`
- Create: `workflow/src/depletion/ctr_insertions.py` — from `def_ctr_insertions.py`
- Modify: those three scripts

**Interfaces:**
- Produces:
  - `weights.Scheme`, `weights.AnnotCol`, `weights.StatsCol` (StrEnums), `weights.WeightInputs`
  - `weights.load_inputs(stats_file: Path, lfc_file: Path, annotations_file: Path, scheme: Scheme) -> WeightInputs` — **signature change, see Step 1**
  - `weights.filter_in_gene(inputs) -> WeightInputs`
  - `weights.neg_log10(values) -> pd.DataFrame`
  - `weights.r2_confidence_frame(stats) -> pd.DataFrame`
  - `weights.raw_weights(scheme, inputs) -> pd.DataFrame`
  - `weights.normalise_per_gene_timepoint(weights, inputs) -> pd.DataFrame`
  - `weights.display_summary(weights) -> None`
  - `gene_level.load_lfc_long(lfc_path) -> pd.DataFrame`
  - `gene_level.load_weights_long(weights_path) -> pd.DataFrame`
  - `gene_level.load_gene_metadata(annotations_path) -> pd.DataFrame`
  - `gene_level.aggregate_to_gene_level(lfc_long, weights, gene_metadata) -> pd.DataFrame`
  - `gene_level.generate_summary(gene_df) -> dict[str, int]`
  - `gene_level.display_summary(stats) -> None`
  - `ctr_insertions.ControlSelectionResult`
  - `ctr_insertions.load_and_preprocess_data(counts_file: Path, annotations_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]`
  - `ctr_insertions.get_control_insertions(counts_df, insertion_annotations) -> pd.DataFrame`

Note `weights.display_summary` and `gene_level.display_summary` share a name in different
modules and are never imported together — that is fine, do not rename either.

- [ ] **Step 1: Decouple `load_inputs` from `WeightsConfig`**

It currently reads `config.stats_file`, `config.lfc_file`, `config.annotations_file`,
`config.scheme`. Replace the parameter with those four as explicit scalars in that order
and substitute the attribute reads with the parameter names. Nothing else in the body
changes. `WeightsConfig` stays in the script; its `main()` call site becomes:

```python
inputs = load_inputs(config.stats_file, config.lfc_file, config.annotations_file, config.scheme)
```

- [ ] **Step 2: Dissolve `process_control_insertions` and `save_results` into `main()`**

`def_ctr_insertions.py`'s `process_control_insertions(config)` both computes and writes,
and `save_results` is a one-line `to_csv`. Neither belongs in a library. Move only
`load_and_preprocess_data` and `get_control_insertions` into
`src/depletion/ctr_insertions.py`; inline their orchestration and the write into `main()`,
preserving the exact same call order, the same `ControlSelectionResult` construction, and
the same output path and `to_csv` arguments. If the existing `to_csv` passes any non-default
argument, carry it over exactly — a changed separator or index flag silently corrupts the
output.

- [ ] **Step 3: Move the remaining functions verbatim**

Per the Interfaces list. `Scheme` moves to `weights.py`; the script's `parse_args()` imports
it for its `--scheme` choices.

- [ ] **Step 4: Confirm modules import and all three scripts run**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python
for m in weights gene_level ctr_insertions; do
  $PY -c "import sys; sys.path.insert(0,'workflow/src'); import depletion.$m; print('$m ok')"
done
```

- [ ] **Step 5: Execute all three against real data and byte-compare**

Flags come from `workflow/rules/depletion_scoring.smk`. Available inputs include
`13_filtered/raw_reads.filtered.tsv`, `13_filtered/control_insertions.tsv`,
`14_insertion_level_depletion_analysis/{LFC,padj}.tsv`,
`15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv`,
`16_gene_level_depletion_analysis/transformed_weights.tsv`. Run each script in the baseline
worktree and the current tree; `diff` every TSV. Report per-file results.

- [ ] **Step 6: Commit**

```bash
git add workflow/src/depletion workflow/scripts/depletion_scoring
git commit -m "refactor(depletion): extract weights, gene_level, ctr_insertions

load_inputs takes 4 scalars instead of WeightsConfig (config dataclasses stay
in scripts). def_ctr_insertions' process_control_insertions/save_results
dissolved into main(), since they mixed compute with I/O.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3.3: Extract imputation and both insertion-level branches

**Files:**
- Create: `workflow/src/depletion/imputation.py` — from `impute_missing_values_using_FR.py`
- Create: `workflow/src/depletion/insertion_level_replicates.py` — from `insertion_level_depletion_analysis_has_replicates.py`
- Create: `workflow/src/depletion/insertion_level_no_replicates.py` — from `insertion_level_depletion_analysis_no_replicates.py`
- Modify: those three scripts

The two insertion-level branches get **separate modules**, not one shared file: the
replicate branch depends on `pydeseq2` while the no-replicate branch is pure pandas, and
merging them would force every consumer to import pydeseq2.

**Interfaces:**
- Produces:
  - `imputation.ImputationResult`
  - `imputation.filter_insertions(insertion_annotations) -> tuple[pd.Index, pd.Index]`
  - `imputation.transfer_FR_index(idxs: tuple) -> tuple`
  - `imputation.impute_missing_values(in_gene_counts_df) -> tuple[pd.DataFrame, list[tuple]]`
  - `imputation.calculate_imputation_statistics(counts_df, in_gene_insertions, ...) -> ImputationResult`
  - `imputation.print_imputation_statistics(result) -> None`
  - `insertion_level_replicates.AnalysisResult`
  - `insertion_level_replicates.load_and_preprocess_data(counts_file, control_insertions_file) -> tuple[...]`
  - `insertion_level_replicates.create_deseq_dataset(counts_df, metadata, control_insertions, initial_timepoint="0h") -> DeseqDataSet`
  - `insertion_level_replicates.perform_differential_analysis(dds, timepoints, initial_timepoint="0h") -> dict[str, DeseqStats]`
  - `insertion_level_replicates.write_dispersion_data_tsv(dds, output_path) -> None`
  - `insertion_level_replicates.concatenate_results(stat_res, timepoints) -> pd.DataFrame`
  - `insertion_level_replicates.transform_index_to_multiindex(dds, layer_name) -> pd.DataFrame`
  - `insertion_level_no_replicates.AnalysisResult`
  - `insertion_level_no_replicates.load_and_preprocess_data(...)`
  - `insertion_level_no_replicates.perform_median_normalization(...)`
  - `insertion_level_no_replicates.calculate_MA_values(...)`

`write_dispersion_data_tsv` writes a file but is a **format serialiser** for a pydeseq2
object, not orchestration — it moves with its module. Both branches define
`load_and_preprocess_data` and `AnalysisResult`; separate modules, so no collision. Do not
rename them.

- [ ] **Step 1: Move the functions and result dataclasses verbatim**

Leave `ImputationConfig` / `InputOutputConfig`, `parse_args()`, `main()`, and
`generate_MA_plots` (dead, per conventions) in the scripts.

- [ ] **Step 2: Confirm modules import**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python
for m in imputation insertion_level_replicates insertion_level_no_replicates; do
  $PY -c "import sys; sys.path.insert(0,'workflow/src'); import depletion.$m; print('$m ok')"
done
```
If `pydeseq2` is missing from this env, try
`/data/a/yangyusheng/miniforge3/envs/pydeseq2/bin/python` for the replicate branch and say
which env you used for which module.

- [ ] **Step 3: Execute against real data and byte-compare**

The replicate branch is the pipeline's most numerically sensitive script (DESeq2 dispersion
estimation). Byte-compare **every** TSV it writes — `14_insertion_level_depletion_analysis/`
holds 11 files. Any single differing file is a STOP condition: report it rather than
proceeding.

- [ ] **Step 4: Commit**

```bash
git add workflow/src/depletion workflow/scripts/depletion_scoring
git commit -m "refactor(depletion): extract imputation and both insertion-level branches

The replicate and no-replicate branches get separate modules: the former needs
pydeseq2, the latter is pure pandas, and merging would force every consumer to
import pydeseq2.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3.4: Move the depletion test and add coverage for the extracted core

**Files:**
- Move: `workflow/scripts/depletion_scoring/test_insertion_level_depletion_analysis_has_replicates.py`
  → `workflow/tests/depletion/test_insertion_level_replicates.py`

- [ ] **Step 1: Move with `git mv` and repoint imports**

Delete the `SCRIPT_DIR` / `sys.path.append` preamble; `pyproject.toml` already provides the
path. Domain symbols import from `depletion.insertion_level_replicates`. If the test asserts
on a config dataclass, use the spec-loader pattern from
`workflow/tests/scripts/conftest.py` — scripts are not importable as modules (no
`__init__.py`, and `src/figures.py` shadows `scripts/figures/`).

- [ ] **Step 2: Run the suite**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests -q
```
Expected: the existing 68 still pass, plus the moved test. **A skip is a failure here** —
Phase 2 shipped 26 silently-skipped tests because fixtures pointed one directory above the
real data. If a test skips, find the data (check `18_figure_data/arc/` and the numbered
stage dirs) and fix the fixture; do not accept the skip.

- [ ] **Step 3: Confirm nothing is left behind, then commit**

```bash
ls workflow/scripts/depletion_scoring/test_*.py 2>/dev/null && echo "FAIL: tests remain" || echo "ok"
git add -A workflow/tests workflow/scripts/depletion_scoring
git commit -m "test(depletion): move insertion-level test to workflow/tests/depletion/

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3.5: Phase 3 verification

- [ ] **Step 1: Byte-compare every TSV, baseline vs current**

```bash
fail=0
for b in $(find /tmp/refactor_baseline/phase3 -name "*.tsv" | sort); do
  n=${b#/tmp/refactor_baseline/phase3/}; a="/tmp/refactor_after/phase3/$n"
  if   [ ! -f "$a" ];    then echo "MISSING:   $n"; fail=1
  elif cmp -s "$b" "$a"; then echo "IDENTICAL: $n ($(stat -c%s "$a") bytes)"
  else                        echo "DIFFERS:   $n"; fail=1
  fi
done
exit $fail
```
Every file must be `IDENTICAL`. A `DIFFERS` on a depletion TSV means the refactor changed
a computed result — stop and fix, never park it.

- [ ] **Step 2: Structural checks**

```bash
grep -rn "^def main\|^def parse_args\|setup_logger(" workflow/src/depletion/ && echo "FAIL" || echo "ok: library-only"
grep -rn "^def sigmoid_function" workflow/ --include="*.py"   # expect exactly ONE hit, in curve_model.py
wc -l workflow/scripts/depletion_scoring/*.py | sort -rn
```

- [ ] **Step 3: Suite plus report**

Run `pytest workflow/tests -q`, then write the before/after line-count table and the
per-file TSV comparison into the report. State explicitly that the Snakemake DAG was not
verified if `snakemake` is unavailable in your env — do not imply a check you did not run.

---

## Phase 4: Read Processing Domain Extraction

Nine scripts, 3,046 lines. This is the hardest phase for two reasons: **every** script is
config-coupled (unlike Phase 2, where one function was), and the upstream scripts consume
BAM/parquet inputs that may no longer exist for the reference project.

### Phase 4 Conventions

Same library rules as Phases 2–3 (verbatim moves; IMPORTS → CONSTANTS → CORE LOGIC; no
`main`/`parse_args`/`setup_logger`; bare sibling imports; config dataclasses stay in
scripts). Three additions specific to this phase:

**1. The config-decoupling rule.** Every `src` function that currently reads `config.x`
must take the values it uses as explicit parameters instead, in the order they appear in
the config dataclass. A `src` module must never import a config dataclass from a script.
Apply this mechanically:

```python
# before, in the script
def build_filter_mask(chunk: pd.DataFrame, config: InputOutputConfig) -> pd.Series:
    ... config.thresholds ... config.paired ...

# after, in src/read_processing/filtering.py
def build_filter_mask(chunk: pd.DataFrame, thresholds: FilterThresholds, paired: bool) -> pd.Series:
    ... thresholds ... paired ...
```

The call site in `main()` passes `config.thresholds, config.paired`. Do not invent new
defaults, do not reorder the remaining parameters, and do not bundle scalars into a new
dataclass — that would just re-create the coupling under another name.

**2. Orchestration functions dissolve into `main()`.** Several functions are pure
orchestration plus file I/O and have no domain logic to extract:
`main_processing_function`, `merge_timepoints`, `concatenate`, `filter_read_pairs`,
`count_insertion_sites`, `process_chunks`, `process_bam_file`, `save_*`. Move their
*inner computation* into `src`, and inline the read/write/orchestration into `main()`,
preserving the exact call order, the exact output paths, and every non-default argument to
`to_csv`/`to_parquet`/`write_table`. A changed separator, index flag, or compression level
silently corrupts output that downstream stages then consume.

When a function is *entirely* I/O plus a single computation call, leave nothing behind in
`src` for it — put the computation in `src` and the I/O in `main()`.

**3. Chunked-streaming scripts keep their loop in the script.** `filter_aligned_reads.py`,
`extract_insertion_sites.py` and `parse_bam_to_tsv.py` stream in chunks/batches. The
per-chunk transform is domain logic and moves to `src`; the chunk loop, the reader, and the
writer stay in `main()`. This keeps `src` functions pure and testable on a single chunk.

**Environment — use `bioinformatics` for everything in this phase.** It is the only env
carrying all four dependencies this phase needs (measured, not assumed):

```bash
PY=/data/a/yangyusheng/miniforge3/envs/bioinformatics/bin/python
$PY -c "import pysam, pyarrow, Bio, pandas; print('all deps present')"
```

Do not reach for `envs/pysam` — despite the name it has pysam but **no pandas and no
pyarrow**, so every script in this phase fails there. There is no `biopython` env; `Bio`
lives in `bioinformatics`, `data_analysis` and `cnsplots_figures`.

For reference, what each script needs: `parse_bam_to_tsv.py` → pysam + pyarrow;
`filter_aligned_reads.py` and `extract_insertion_sites.py` → pyarrow;
`concatenate_timepoint_data.py` → Bio; the rest → pandas only.

**Verification, and an honest statement of its limits.** Capture a baseline worktree first
(`git worktree add /tmp/phase4_baseline <pre-phase-4-sha> --detach`), then for each script
run it in both trees and byte-compare outputs.

Available real data under `projects/HD_DIT_HAP/results/`: `6_filtered/` (36 files),
`7_insertions/` (36), `8_merged/` (18), `9_concatenated/` (9), `10_annotated/` (3),
`11_merged/` (3), `12_concatenated/` (2), `13_filtered/` (4). Note `1_fastp/` is EMPTY and
`3_mapped/` holds a single file.

So the three upstream scripts (`parse_bam_to_tsv`, `filter_aligned_reads`,
`extract_insertion_sites`) may have no runnable input. For those:
- If an input exists, run and byte-compare it.
- If not, **say so explicitly in the report, naming the missing input**, and fall back to
  unit-testing the extracted per-chunk functions on a synthetic chunk you construct.
- Never write "verified" for a script you did not execute. A phase report that overstates
  coverage is worse than one that admits a gap — Phase 1 was reported as verified and had
  shipped 10 broken scripts.

---

### Task 4.1: Extract the four pandas-only scripts

Batched: same shape, same env, no heavy dependencies.

**Files:**
- Create: `workflow/src/read_processing/__init__.py` (`"""Read processing library modules."""`)
- Create: `workflow/src/read_processing/annotation.py` — from `annotate_genomic_features.py`
- Create: `workflow/src/read_processing/merge.py` — from `merge_strand_insertions.py` + `merge_similar_timepoints.py`
- Create: `workflow/src/read_processing/concat.py` — from `concat_counts_and_annotations.py`
- Create: `workflow/src/read_processing/hard_filtering.py` — from `reads_hard_filtering.py`
- Modify: those five scripts

**Interfaces:**
- `annotation.AnalysisResult`, `annotation.load_insertion_data(input_path) -> pd.DataFrame`,
  `load_genome_regions(region_path)`, `calculate_codon_distances_vectorized(annotated_df)`,
  `calculate_affected_residue_vectorized(annotated_df)`,
  `assign_insertion_direction_vectorized(annotated_df) -> pd.Series`,
  `drop_boundary_duplicates(sub_df)`, `annotate_insertions(...)`
- `merge.MergeResult`, `merge.merge_insertion_data(pbl_df, pbr_df) -> pd.DataFrame`
- `concat` — see Step 2; `concatenate(config)` is orchestration
- `hard_filtering.AnalysisResult`, `hard_filtering.load_insertion_data(input_file)`,
  `validate_timepoint_exists(df, timepoint) -> None`,
  `apply_hard_filtering(df, timepoint: str, cutoff: float, ...) -> tuple[pd.DataFrame, AnalysisResult]`
  — **decoupled from `InputOutputConfig` per convention 1; read the dataclass for the exact
  field names and order, and pass exactly the fields the body uses**

`merge_strand_insertions.py` and `merge_similar_timepoints.py` share one module because both
are merge operations on insertion tables; their functions do not collide.

- [ ] **Step 1: Decouple every config-reading function** per convention 1.
- [ ] **Step 2: Dissolve orchestration** — `main_processing_function`, `save_annotations`,
      `save_merged_data`, `merge_timepoints`, `concatenate` — into each script's `main()`
      per convention 2, preserving output paths and all `to_csv` arguments exactly.
- [ ] **Step 3: Move the remaining functions verbatim.**
- [ ] **Step 4: Confirm modules import**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python
for m in annotation merge concat hard_filtering; do
  $PY -c "import sys; sys.path.insert(0,'workflow/src'); import read_processing.$m; print('$m ok')"
done
```

- [ ] **Step 5: Execute all five scripts and byte-compare.** Flags from
      `workflow/rules/read_processing.smk`. Inputs: `9_concatenated/`, `10_annotated/`,
      `11_merged/`, `12_concatenated/`, `13_filtered/`. Report per-file diffs.
- [ ] **Step 6: Commit.**

---

### Task 4.2: Extract the sequence-extraction and chunked-streaming scripts

**Files:**
- Create: `workflow/src/read_processing/sequence_extraction.py` — from `concatenate_timepoint_data.py`
- Create: `workflow/src/read_processing/insertions.py` — from `extract_insertion_sites.py`
- Create: `workflow/src/read_processing/filtering.py` — from `filter_aligned_reads.py`
- Create: `workflow/src/read_processing/bam.py` — from `parse_bam_to_tsv.py`
- Modify: those four scripts

**Interfaces:**
- `sequence_extraction.AnalysisResult`, `load_reference_data(genome_path) -> dict`,
  `extract_target_sequence(chrom: str, coordinate: int, ref_dict: dict) -> str`,
  `process_concatenation_data(...)` — decoupled from config
- `insertions.InsertionCounts` (the `type` alias — it moves with the functions that use it),
  `insertions.ExtractionStats`, `calculate_insertion_coordinates_vectorized(valid_df)`,
  `create_validation_mask(df) -> pd.Series`,
  `count_insertions_vectorized(valid_df) -> InsertionCounts`,
  `extract_insertion_sites(chunk, chunk_num) -> tuple[InsertionCounts, int, int]`,
  `create_output_dataframe(insertion_counts) -> tuple[pd.DataFrame, int, int]`
- `filtering.FilterThresholds`, `filtering.AnalysisResult`,
  `load_config_from_yaml(config_file) -> dict[str, Any]`,
  `strip_read_prefix(column: str) -> str`, `coerce_column_dtypes(chunk) -> pd.DataFrame`,
  `build_filter_mask(...)`, `process_chunk(...)` — both decoupled from config
- `bam.ReadInfo`, `bam.ReadPairInfo`,
  `extract_read_info(read: pysam.AlignedSegment | None, tag_list) -> ReadInfo`,
  `determine_proper_pair_status(...)`, `process_read_pair(...)`,
  `format_output_line(pair_info, tag_list) -> list`, `build_header(tag_list) -> list[str]`,
  `build_schema(header_fields) -> pa.Schema`, `flush_batch(...)`

`process_chunks`, `count_insertion_sites`, `write_empty_output`, `filter_read_pairs` and
`process_bam_file` are the chunk loops and stay in their scripts per convention 3.

- [ ] **Step 1: Move the per-chunk/per-record transforms to `src`, decoupled from config.**
- [ ] **Step 2: Keep the streaming loops, readers and writers in `main()`.**
- [ ] **Step 3: Confirm each module imports**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/bioinformatics/bin/python
for m in sequence_extraction insertions filtering bam; do
  $PY -c "import sys; sys.path.insert(0,'workflow/src'); import read_processing.$m; print('$m ok')"
done
```
Expected: four `ok` lines. `bioinformatics` carries pysam, pyarrow, Bio and pandas together,
so one env covers all four modules.

- [ ] **Step 4: Execute what can be executed; unit-test the rest.**
      `concatenate_timepoint_data.py` has inputs (`8_merged/`, `9_concatenated/`).
      The other three may not — follow the honesty rule in the conventions: name the missing
      input, then unit-test the extracted transforms on a synthetic chunk.
- [ ] **Step 5: Commit.**

---

### Task 4.3: Add tests for the extracted read-processing core

**Files:**
- Create: `workflow/tests/read_processing/__init__.py`
- Create: `workflow/tests/read_processing/test_insertions.py`
- Create: `workflow/tests/read_processing/test_filtering.py`
- Create: `workflow/tests/read_processing/test_annotation.py`

`read_processing` had **zero** tests before this refactor. The extraction is the first point
at which its logic is testable in isolation, and the three upstream scripts may be
unverifiable by output comparison — which makes unit tests the only safety net they get.

- [ ] **Step 1: Test the pure transforms on synthetic frames.** Cover at minimum:
      `create_validation_mask` (a row that passes, a row that fails each condition),
      `calculate_insertion_coordinates_vectorized` (coordinate arithmetic on both strands),
      `count_insertions_vectorized` (two reads at one site aggregate; two sites stay
      separate), `strip_read_prefix`, `coerce_column_dtypes`, `build_filter_mask` (one row
      inside every threshold, one outside each), `assign_insertion_direction_vectorized`,
      `drop_boundary_duplicates`.

      Each test asserts a specific expected value, not merely "does not crash". A test that
      only checks for absence of exceptions is not a test.

- [ ] **Step 2: Run and confirm no skips**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests -q
```
If a test needs pyarrow/pysam and that env lacks it, mark it with
`pytest.importorskip("pyarrow")` at the top — but report how many tests that affects, since
a skipped test protects nothing.

- [ ] **Step 3: Commit.**

---

### Task 4.4: Phase 4 verification

- [ ] **Step 1: Byte-compare every output produced in both trees.** Same loop as Task 3.5,
      with `phase4` paths. Every file `IDENTICAL`; any `DIFFERS` is a STOP.
- [ ] **Step 2: Structural checks**

```bash
grep -rn "^def main\|^def parse_args\|setup_logger(" workflow/src/read_processing/ && echo FAIL || echo "ok: library-only"
grep -rn "config\." workflow/src/read_processing/ && echo "FAIL: config leaked into src" || echo "ok: no config coupling"
wc -l workflow/scripts/read_processing/*.py | sort -rn
```
The second check is the one that matters most this phase: it proves convention 1 held.

- [ ] **Step 3: Report** — the before/after line-count table, the per-file comparison, the
      list of scripts that could NOT be executed with the specific missing input for each,
      and the unit-test count that covers them instead. State plainly whether the Snakemake
      DAG was checked.

---

## Phase 5: Quality Control Domain Extraction

Six scripts, 1,948 lines. All emit TSVs, so outputs are byte-comparable. Config coupling is
moderate (3–12 references per script) — apply Phase 4's convention 1 mechanically.

### Phase 5 Conventions

Same as Phase 4's conventions (verbatim moves, library-only `src` modules, bare sibling
imports, config dataclasses stay in scripts, config-decoupling by explicit scalars,
orchestration into `main()`), with one env simplification:

**Environment.** `/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python`
covers every script here (pandas + numpy only; no pysam, pyarrow or Bio needed).

**Baseline.** `git worktree add /tmp/phase5_baseline <pre-phase-5-sha> --detach` before
extracting. Real inputs: `projects/HD_DIT_HAP/results/{9_concatenated,12_concatenated,13_filtered,14_insertion_level_depletion_analysis,16_gene_level_depletion_analysis}/`
and `projects/HD_DIT_HAP/results/18_figure_data/arc/` (which holds this phase's own prior
outputs: `gene_coverage_stats.tsv`, `pbl_pbr_pairs.tsv`, `strand_pairs.tsv`,
`read_count_distribution.tsv`, `read_count_cutoff_stats.tsv`, `insertion_density_analysis.tsv`).
Do not conclude an input is missing without checking `arc/` — that mistake cost a full
verification round in Phase 2.

`extract_mapping_filtering_statistics.py` reads **log files**, not TSVs; its inputs are
under `projects/HD_DIT_HAP/logs/`. Check what exists there before claiming it is unrunnable.

---

### Task 5.1: Extract the four small QC modules

**Files:**
- Create: `workflow/src/qc/__init__.py` (`"""Quality control library modules."""`)
- Create: `workflow/src/qc/correlation.py` — from `PBL_PBR_correlation_analysis.py`
- Create: `workflow/src/qc/orientation.py` — from `insertion_orientation_analysis.py`
- Create: `workflow/src/qc/coverage.py` — from `gene_coverage_analysis.py`
- Create: `workflow/src/qc/read_counts.py` — from `read_count_distribution_analysis.py`
- Modify: those four scripts

**Interfaces:**
- `correlation.read_tsv_file(file_path) -> pd.DataFrame | None`,
  `correlation.parse_filename(file_path) -> tuple[str, str, str] | None`
- `orientation.extract_strand_pairs(df) -> pd.DataFrame`
- `coverage.ViabilityCol` (StrEnum), `coverage.CoverageStat`,
  `coverage.load_covered_genes(lfc_file, annotation_file) -> set[str]`,
  `coverage.load_gene_viability(viability_file) -> pd.DataFrame`,
  `coverage.compute_coverage_stats(viability, covered_genes) -> list[CoverageStat]`,
  `coverage.write_coverage_table(stats, output_file) -> None` — a format serialiser for
  `list[CoverageStat]`, so it moves with its module rather than dissolving into `main()`
- `read_counts.load_and_validate_data(file_path) -> pd.DataFrame`,
  `read_counts.parse_sample_name(file_path) -> str`,
  `read_counts.calculate_cutoff_statistics(df, initial_time_point: str, cutoff: float) -> dict[str, float | int]`,
  `read_counts.compute_binned_distribution(df, sample: str, bins: int) -> pd.DataFrame`,
  `read_counts.log_summary_table(stats_df) -> None`

Note `qc/correlation.py` and `figure_render/correlation.py` share a filename in different
packages — that is fine and intentional (compute vs render). Same for `qc/coverage.py` /
`figure_render/coverage.py`, `qc/orientation.py` / `figure_render/orientation.py`, and
`qc/read_counts.py` / `figure_render/read_counts.py`. Do not rename to avoid the echo; the
package prefix disambiguates and the parallel naming is a feature.

- [ ] **Step 1: Decouple config-reading functions** per Phase 4 convention 1.
- [ ] **Step 2: Move the listed symbols verbatim; dissolve orchestration into `main()`.**
- [ ] **Step 3: Confirm modules import**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python
for m in correlation orientation coverage read_counts; do
  $PY -c "import sys; sys.path.insert(0,'workflow/src'); import qc.$m; print('$m ok')"
done
```

- [ ] **Step 4: Execute all four and byte-compare** against the baseline worktree. Flags from
      `workflow/rules/quality_control.smk`.
- [ ] **Step 5: Commit.**

---

### Task 5.2: Extract insertion density and mapping statistics

`insertion_density_analysis.py` is this phase's largest script (601 lines) and holds the
most computation: six independent statistic calculators plus a Gini coefficient.

**Files:**
- Create: `workflow/src/qc/density.py` — from `insertion_density_analysis.py`
- Create: `workflow/src/qc/mapping_stats.py` — from `extract_mapping_filtering_statistics.py`
- Modify: those two scripts

**Interfaces:**
- `density.AnalysisResult`, `density.load_insertion_data(...)` (decoupled from config),
  `load_annotation_data(annotations_path) -> pd.DataFrame`,
  `filter_in_gene_insertions(insertion_data, ...) -> pd.DataFrame`,
  `calculate_insertion_statistics(gene_insertions) -> dict[str, int | float]`,
  `calculate_gap_statistics(gene_insertions) -> dict[str, int | float | str]`,
  `calculate_gini_coefficient(values: np.ndarray) -> float`,
  `calculate_read_statistics(gene_insertions) -> dict[str, int | float]`,
  `calculate_strand_statistics(gene_insertions) -> dict[str, int | float]`,
  `analyze_gene_insertions(gene_id: str, gene_insertions) -> dict[str, str | int | float]`,
  `generate_summary_statistics(results_df) -> dict[str, int | float]`
- `mapping_stats.FilteringStatistics`, `mapping_stats.AnalysisResult`,
  `parse_log_file(log_file: Path) -> dict[str, FilteringStatistics]`,
  `extract_summary_data(log_files: list[Path]) -> dict[str, FilteringStatistics]`,
  `create_dataframe(statistics: dict[str, FilteringStatistics]) -> pd.DataFrame`

`save_results` in `extract_mapping_filtering_statistics.py` is orchestration + I/O →
dissolve into `main()`.

- [ ] **Step 1: Move verbatim, decoupling config reads.** `calculate_gini_coefficient` is
      pure maths on an ndarray — it must move unchanged; a subtle edit there would silently
      shift every gene's evenness score.
- [ ] **Step 2: Confirm modules import.**
- [ ] **Step 3: Execute both and byte-compare.** Density inputs:
      `12_concatenated/` or `13_filtered/` plus `10_annotated/`; compare against
      `18_figure_data/arc/insertion_density_analysis.tsv` if the rule writes there.
      Mapping stats reads logs under `projects/HD_DIT_HAP/logs/` — list that directory first
      and report what you find rather than assuming.
- [ ] **Step 4: Commit.**

---

### Task 5.3: Add tests for the extracted QC core, then verify

`quality_control` had **zero** tests. The statistic calculators are pure functions on small
frames — the easiest high-value tests in the whole refactor.

**Files:**
- Create: `workflow/tests/qc/__init__.py`
- Create: `workflow/tests/qc/test_density.py`
- Create: `workflow/tests/qc/test_read_counts.py`
- Create: `workflow/tests/qc/test_coverage.py`

- [ ] **Step 1: Write tests asserting specific values, not absence of exceptions.** Cover at
      minimum:
      - `calculate_gini_coefficient`: a perfectly even array → 0.0 (within 1e-12); a
        maximally skewed array → approaches 1.0; a single-element array (document whatever
        it currently returns rather than asserting a value you think is right).
      - `calculate_gap_statistics`: two insertions with a known coordinate gap → that exact
        gap; a single insertion → whatever the current code yields for the no-gap case.
      - `calculate_strand_statistics`: 3 plus / 1 minus → the exact ratio the code computes.
      - `calculate_cutoff_statistics`: a frame straddling the cutoff → exact retained count.
      - `compute_binned_distribution`: known values and `bins=4` → expected bin counts.
      - `compute_coverage_stats`: 2 covered of 3 viable genes → exact `CoverageStat` fields.

      Derive every expected value by reading the implementation, not by running it and
      pasting the output — a test that merely records current behaviour cannot detect that
      the behaviour was already wrong. Where you cannot derive it confidently, say so in the
      report and assert the looser property you *can* justify.

- [ ] **Step 2: Run the suite; no skips.**
- [ ] **Step 3: Phase 5 verification** — byte-compare loop (as Task 3.5, `phase5` paths),
      then:

```bash
grep -rn "^def main\|^def parse_args\|setup_logger(" workflow/src/qc/ && echo FAIL || echo "ok: library-only"
grep -rn "config\." workflow/src/qc/ && echo "FAIL: config leaked" || echo "ok: decoupled"
wc -l workflow/scripts/quality_control/*.py | sort -rn
```

- [ ] **Step 4: Commit** with the before/after table and per-file comparison in the report.

---

## Phase 6: Reference Data Domain Extraction

One script, 645 lines — the single largest in the tree, and the most config-coupled (18
`config.` references). It builds the PomBase genome-region BED files every other stage
depends on.

### Task 6.1: Extract genome region processing

**Files:**
- Create: `workflow/src/reference_data/__init__.py` (`"""Reference data library modules."""`)
- Create: `workflow/src/reference_data/genome_region.py`
- Modify: `workflow/scripts/reference_data/extract_genome_region.py`

**Interfaces:**
- Consumes: `gene_metadata.resolve_gene_ids` (already wired in Phase 1 — this script's local
  `update_sysID` was replaced there; confirm it is still importing the shared one and did not
  regain a local copy), `io_tables.read_table`, `logging_setup.setup_logger`
- Produces (moved verbatim):
  - `get_gff_transcript_id(...)`
  - `parse_gff_data(gff_file_path: Path) -> pd.DataFrame`
  - `calculate_accumulated_cds_bases(transcript_features_df) -> pd.DataFrame`
  - `gff_features_to_bed(...)`
  - `select_primary_transcripts(all_coding_features_bed_df) -> pd.DataFrame`
  - `build_intergenic_bed(primary_transcripts_bed_df, fai_file_path: Path) -> pd.DataFrame`
  - `annotate_intergenic_region_flanks(...)`
  - `find_overlapping_regions(feature) -> pybedtools.Interval`
  - `build_genome_intervals(...)`

`run_pipeline(config)` is orchestration: it reads 18 config fields and writes 6 output files.
Dissolve it into `main()` per Phase 4 convention 2 — the computation steps move to the module,
the sequencing and every `to_csv`/BED write stay in `main()`, preserving output paths and
arguments exactly. This is the highest-risk dissolution in the refactor because six
downstream stages consume these BEDs; if `run_pipeline` is too tangled to split confidently,
report that rather than guessing.

**Environment.** Needs **pybedtools**: `/data/a/yangyusheng/miniforge3/envs/pybedtools/bin/python`
(has pybedtools + pysam + pandas; verified). `bioinformatics` does NOT have pybedtools.

- [ ] **Step 1: Confirm Phase 1's `resolve_gene_ids` wiring survived**

```bash
grep -n "resolve_gene_ids\|def update_sysID" workflow/scripts/reference_data/extract_genome_region.py
```
Expected: an import/call of `resolve_gene_ids`, and NO local `def update_sysID`. If a local
copy is back, report it — that would mean a later commit reverted Phase 1.

- [ ] **Step 2: Move the listed functions verbatim into the module.**
- [ ] **Step 3: Dissolve `run_pipeline` into `main()`**, decoupling the 18 config reads into
      explicit parameters on the moved functions.
- [ ] **Step 4: Confirm module imports and CLI constructs**

```bash
PY=/data/a/yangyusheng/miniforge3/envs/pybedtools/bin/python
$PY -c "import sys; sys.path.insert(0,'workflow/src'); import reference_data.genome_region; print('ok')"
$PY workflow/scripts/reference_data/extract_genome_region.py --help > /dev/null && echo "help ok"
```

- [ ] **Step 5: Execute and byte-compare the BED outputs**

This script's outputs live in `resources/pombase_data/<release>/genome_region/`. The
`rule all` target is `coding_gene_primary_transcripts.bed`. Check what already exists:

```bash
ls resources/pombase_data/*/genome_region/ 2>/dev/null
```
If the outputs and the GFF/fai inputs are present, run the script in the baseline worktree
and the current tree and `diff` every BED. If the inputs are absent, say so explicitly and
name them — do not report this script as verified on the strength of `--help` alone. That is
precisely the mistake Phase 1 made across 10 scripts.

- [ ] **Step 6: Commit.**

---

### Task 6.2: Final whole-refactor verification

- [ ] **Step 1: Structural invariants across the whole tree**

```bash
echo "--- every src module is library-only ---"
grep -rn "^def main\|^def parse_args\|setup_logger(" workflow/src/ --include="*.py" \
  | grep -v "logging_setup.py:.*def setup_logger" && echo FAIL || echo ok

echo "--- no config coupling anywhere in src ---"
grep -rn "config\." workflow/src/ --include="*.py" && echo FAIL || echo ok

echo "--- setup_logger defined exactly once ---"
test "$(grep -rc "^def setup_logger" workflow/src/logging_setup.py)" = "1" && echo ok || echo FAIL
grep -rn "^def setup_logger\|^def setup_logging" workflow/scripts/ --include="*.py" && echo "FAIL: script-local copy" || echo ok

echo "--- sigmoid_function defined exactly once ---"
grep -rn "^def sigmoid_function" workflow/ --include="*.py"

echo "--- no test files left under scripts/ ---"
find workflow/scripts -name "test_*.py" | grep . && echo FAIL || echo ok

echo "--- every script is thin ---"
find workflow/scripts -name "*.py" | xargs wc -l | sort -rn | head -12
```

- [ ] **Step 2: Full suite**

```bash
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/pytest workflow/tests -q
```
Report the count. Any skip must be justified in the report with the reason and the data
checked (including `18_figure_data/arc/`).

- [ ] **Step 3: Write the final summary table** — per phase: scripts touched, lines before,
      lines after, modules created, outputs byte-compared, outputs that could not be
      compared and why.

- [ ] **Step 4: State the DAG position plainly.** If `snakemake` is unavailable in every
      reachable env, say the DAG was never validated by this refactor and name the mitigation
      (no `.smk` file modified; CLI flags byte-identical pre/post). Do not imply a check that
      did not happen.

- [ ] **Step 5: Deferred-cleanup list.** Collect every minor finding parked across all phases
      into one section of the report so it can be triaged as a follow-up:
      duplicate `SCRIPT_DIR` in `insertion_level_depletion_analysis_no_replicates.py`;
      possibly-unused `INDEX_COLUMNS` in `compute_insertion_weights.py`; unused local `fig`
      in `figure_render/correlation.py`; the three dead matplotlib functions in
      `depletion_scoring/` (`create_fitted_plot`, `generate_fitting_plots`,
      `generate_MA_plots`); the two scripts no rule invokes (`plot_dispersions.py`,
      `plot_ma_plot_replicates.py`); and fixtures hardcoding `18_figure_data/arc/` rather
      than preferring fresh pipeline output.

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

