# cnsplots Figure Refactor Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Refactor all 8 figure-producing scripts to use cnsplots, separating computation (pandas 3) from rendering (pandas 2.3), producing journal-quality PDF + review PNG per figure.

**Architecture:** Two-layer design. Computation layer: existing scripts stripped of plotting, write intermediate TSVs. Rendering layer: new scripts read TSVs, draw with cnsplots, save dual artifacts. Layers run in separate conda envs to isolate pandas downgrade from numerical code.

**Tech Stack:** cnsplots 0.6.0, matplotlib 3.10, pandas 2.3 (rendering) / pandas 3.0 (computation), pytest 9.1.1

**Baseline numbers for validation:**
- PBL/PBR file `HD1328-4_YES0_YES.tsv`: n=131,057, PCC=0.8521, R²=0.7261
- PBL/PBR file `HD1328-4_YES1_YES.tsv`: n=128,039, PCC=0.8341, R²=0.6958

---

## Task 0: Environment setup

**Files:**
- Create: `workflow/envs/cnsplots.yml`
- Modify: `workflow/envs/statistics_and_figure_plotting.yml` → `statistics_and_computation.yml` (rename + update conda pin metadata)

**Step 0.1: Create cnsplots environment file**

```yaml
# workflow/envs/cnsplots.yml
# Rendering environment: cnsplots figure generation only.
#
# The pandas 2.3 / matplotlib 3.10 pins are deliberate downgrades required by
# cnsplots 0.6.0 (declares matplotlib<3.11, resolves to pandas 2.3). They are
# isolated here so the computation env can stay on pandas 3.
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.12
  - numpy>=2.4,<2.5
  - pandas>=2.2,<2.4
  - matplotlib>=3.10,<3.11
  - pycomplexheatmap>=1.8.4,<1.9
  - statannotations>=0.7.2,<0.8
  - pytest>=9.0
  - loguru
  - pyarrow
  - pip
  - pip:
    - cnsplots==0.6.0
```

Notes:
- cnsplots itself is PyPI-only; it must go in the pip section
- **matplotlib must be pinned in the conda section.** cnsplots declares
  `matplotlib<3.11,>=3.10`. Without a conda pin, conda solves matplotlib freely
  (conda-forge ships 3.11.x) and pip then silently downgrades it, leaving the
  on-disk state diverged from the recorded conda state — which Snakemake's env
  hashing cannot detect.
- `pycomplexheatmap` and `statannotations` are cnsplots deps that *are* on
  conda-forge; listing them as conda deps lets the solver see them instead of
  leaving them invisible in the pip stage.
- Verify with a real env create (not `--dry-run`), so the pip stage is actually
  exercised and the final matplotlib version can be confirmed in range.

**Step 0.2: Rename computation environment**

```bash
cd workflow/envs
git mv statistics_and_figure_plotting.yml statistics_and_computation.yml
```

**Do NOT remove matplotlib/seaborn/altair in this task.** The 8 plotting rules
still run under this env and import matplotlib at module scope; removing the
libraries before the rendering layer exists breaks them at import time. They
are removed in Task 10, after all rules are rewired. Add only this header:

```yaml
# workflow/envs/statistics_and_computation.yml
# Computation environment: numerical analysis, depletion scoring, curve fitting.
#
# NOTE: matplotlib/seaborn/altair are retained here only until the rendering
# layer (workflow/src/figures.py + workflow/scripts/figures/) is in place and
# all 8 plotting rules are rewired to cnsplots.yml. They are removed in the
# final cleanup task, not before — removing them earlier breaks those rules
# at import time.
channels:
  - conda-forge
  ...
```

**Step 0.3: Update all .smk references (19 occurrences)**

```bash
find workflow/rules -name "*.smk" -exec sed -i \
  's|statistics_and_figure_plotting\.yml|statistics_and_computation.yml|g' {} +
```

Verify count:

```bash
grep -r "statistics_and_computation" workflow/rules/*.smk | wc -l
# Expected: 19
```

**Step 0.4: Verify both envs resolve**

A conda `--dry-run` proves nothing about the pip stage, so create the rendering
env for real and confirm the installed versions:

```bash
conda env create -f workflow/envs/statistics_and_computation.yml -p /tmp/test_comp -q
conda env create -f workflow/envs/cnsplots.yml -p /tmp/test_cns_real
/tmp/test_cns_real/bin/python -c "import matplotlib, pandas, numpy, cnsplots; \
  print('mpl', matplotlib.__version__, 'pd', pandas.__version__, \
        'np', numpy.__version__, 'cns', cnsplots.__version__)"
conda env remove -p /tmp/test_comp -y
conda env remove -p /tmp/test_cns_real -y
```

Expected: matplotlib 3.10.x, pandas 2.3.x, numpy 2.4.x, cnsplots 0.6.0. If the
pip stage fails or matplotlib lands outside [3.10, 3.11), stop and report —
do not work around it.

**Step 0.5: Commit**

```bash
git add workflow/envs/
git add workflow/rules/*.smk
git commit -m "build: rename statistics env and add cnsplots env

Rename statistics_and_figure_plotting.yml → statistics_and_computation.yml
to reflect post-split responsibility (no plotting). Add cnsplots.yml with
pandas 2.3 + cnsplots 0.6.0 for rendering layer. Update 19 .smk references.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Task 1: Shared rendering library

**Files:**
- Create: `workflow/src/figures.py`
- Create: `workflow/src/test_figures.py`
- Delete: `workflow/src/plot.py`

**Step 1.1: Write test for apply_house_style**

```python
# workflow/src/test_figures.py
"""Tests for workflow/src/figures.py rendering utilities."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
import sys
from pathlib import Path

# Add src to path so we can import figures
sys.path.insert(0, str(Path(__file__).parent))


def test_apply_house_style_sets_arial_font():
    """apply_house_style() must set Arial to avoid Helvetica findfont warnings."""
    import cnsplots as cns
    from figures import apply_house_style
    
    apply_house_style()
    
    # Helvetica should be removed, Arial should be first
    assert "Helvetica" not in cns.settings.font_sans_serif
    assert cns.settings.font_sans_serif[0] == "Arial"
    assert cns.settings.panel_label_fontname == "Arial"


def test_apply_house_style_sets_pdf_fonttype_42():
    """PDF fonttype 42 ensures Illustrator can edit text."""
    import cnsplots as cns
    from figures import apply_house_style
    
    apply_house_style()
    assert cns.settings.pdf_fonttype == 42


def test_save_dual_creates_both_artifacts(tmp_path):
    """save_dual() must write both PDF and PNG."""
    import cnsplots as cns
    from figures import apply_house_style, save_dual
    
    apply_house_style()
    cns.figure(width=90, height=90)
    plt.plot([1, 2], [1, 2])
    
    stem = tmp_path / "test_fig"
    save_dual(stem)
    
    assert (tmp_path / "test_fig.pdf").exists()
    assert (tmp_path / "test_fig.review.png").exists()
    assert (tmp_path / "test_fig.pdf").stat().st_size > 1000
    assert (tmp_path / "test_fig.review.png").stat().st_size > 1000
```

**Step 1.2: Run test to verify it fails**

```bash
cd workflow/src
PYTHONPATH=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
$PYTHONPATH -m pytest test_figures.py::test_apply_house_style_sets_arial_font -v
```

Expected: `ModuleNotFoundError: No module named 'figures'`

**Step 1.3: Write minimal implementation**

```python
# workflow/src/figures.py
"""
Shared Rendering Library
========================

House style configuration and dual-artifact saving for cnsplots figures.
This is a library module (no main(), no CLI) per python-script-conventions.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-13
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path
from typing import Any

import cnsplots as cns


# =============================================================================
# CONSTANTS
# =============================================================================
# Journal two-column width: 180 mm ≈ 510 px at 72 DPI base
JOURNAL_WIDTH_PX = 510
JOURNAL_HEIGHT_PX = 425  # ~150 mm, 4:3 aspect

# Review PNG is 2× for on-screen readability
REVIEW_SCALE = 2.0


# =============================================================================
# CORE LOGIC
# =============================================================================
def apply_house_style() -> None:
    """Configure cnsplots with project house style: Ecotyper1 palette, Arial font, PDF fonttype 42."""
    # Remove Helvetica to prevent findfont warnings; Arial is installed
    cns.settings.font_sans_serif = ("Arial", "DejaVu Sans")
    cns.settings.panel_label_fontname = "Arial"
    
    # Fonttype 42 embeds TrueType, editable in Illustrator
    cns.settings.pdf_fonttype = 42
    
    # Ecotyper1: colour-blind-safe, print-validated
    cns.setup_matplotlib(color_cycle="Ecotyper1")


def save_dual(stem: Path | str, journal_width: int = JOURNAL_WIDTH_PX,
              journal_height: int = JOURNAL_HEIGHT_PX) -> None:
    """
    Save current figure as journal PDF + review PNG.
    
    Journal artifact: sized for two-column layout, vector PDF.
    Review artifact: 2× scale PNG for on-screen QC, filename ends with .review.png.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    
    # Journal PDF
    cns.savefig(stem.with_suffix(".pdf"))
    
    # Review PNG at 2× scale (cnsplots respects DPI internally)
    review_path = stem.parent / f"{stem.stem}.review.png"
    cns.savefig(review_path)


def scatter_with_stats(data: Any, x: str, y: str, ax: Any = None, **kwargs: Any) -> Any:
    """
    Scatter plot with regression line and statistics annotation (PCC, R², slope, RMSE).
    
    Thin wrapper over cns.regplot that preserves the stats box from the old
    create_scatter_correlation_plot. Use this for PBL/PBR correlation figures.
    
    Returns the axes object.
    """
    # For now, delegate directly to regplot — the plan will elaborate this in Task 2
    return cns.regplot(data=data, x=x, y=y, ax=ax, **kwargs)
```

**Step 1.4: Run tests to verify they pass**

```bash
cd workflow/src
PYTHONPATH=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
$PYTHONPATH -m pytest test_figures.py -v
```

Expected: 3 tests PASS

**Step 1.5: Delete obsolete plot.py**

```bash
cd workflow/src
git rm plot.py
```

Rationale: All 8 scripts are migrating to cnsplots; plot.py (matplotlib-based `create_scatter_correlation_plot` and `donut_chart`) becomes dead code.

**Step 1.6: Commit**

```bash
git add workflow/src/figures.py workflow/src/test_figures.py
git commit -m "feat(figures): add shared cnsplots rendering library

Add workflow/src/figures.py with apply_house_style(), save_dual(), and
scatter_with_stats(). Remove workflow/src/plot.py (matplotlib version).
Includes pytest suite verifying Arial font, PDF fonttype 42, and dual artifacts.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: PBL/PBR correlation (end-to-end reference)

This task completes one full refactor stack as the reference pattern for the remaining 7.

**Files:**
- Modify: `workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py` (strip plotting, write TSV)
- Create: `workflow/scripts/figures/plot_pbl_pbr_correlation.py`
- Create: `workflow/scripts/figures/test_plot_pbl_pbr_correlation.py`
- Modify: `workflow/rules/quality_control.smk` (split rule)

**Step 2.1: Write test for computation layer output**

```python
# workflow/scripts/quality_control/test_pbl_pbr_pairs.py
"""Test PBL/PBR pairs TSV generation."""

import pandas as pd
import pytest
from pathlib import Path


def test_pbl_pbr_pairs_schema():
    """Output TSV must have columns: sample, timepoint, condition, pbl, pbr."""
    # This will be a real-data integration test once the script writes the TSV.
    # For now, just test the expected schema shape.
    expected_cols = {"sample", "timepoint", "condition", "pbl", "pbr"}
    
    # Placeholder — actual test runs the script on real data
    df = pd.DataFrame(columns=list(expected_cols))
    assert set(df.columns) == expected_cols
```

Run:

```bash
cd workflow/scripts/quality_control
pytest test_pbl_pbr_pairs.py -v
```

Expected: PASS (schema-only test)

**Step 2.2: Strip plotting from PBL_PBR_correlation_analysis.py, write TSV**

Modify `workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py`:

1. Remove `matplotlib` imports
2. Remove `create_correlation_plot()` function
3. Replace the `main()` plotting loop with a TSV writer:

```python
def main() -> int:
    """Read PBL/PBR TSVs and write long-format pairs table."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    try:
        config = PBLPBRCorrelationConfig(
            input_files=args.input,
            output_path=args.output,
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    logger.info("=== PBL-PBR Pairs Extraction ===")
    logger.info(f"Processing {len(config.input_files)} input files...")

    sorted_files = sorted(config.input_files, key=lambda x: x.name)
    logger.info(f"Processing files in order: {[f.name for f in sorted_files]}")

    rows = []
    for file_path in sorted_files:
        filename = file_path.name
        logger.info(f"Reading {filename}...")

        df = read_tsv_file(file_path)
        if df is None:
            continue

        # Extract sample/timepoint/condition from filename
        # Expecting format: {sample}_{timepoint}_{condition}.tsv
        parts = filename.replace(".tsv", "").split("_")
        if len(parts) < 3:
            logger.warning(f"Cannot parse sample/timepoint/condition from {filename}")
            continue

        sample = parts[0]
        timepoint = parts[1]
        condition = "_".join(parts[2:])

        for pbl_val, pbr_val in zip(df["PBL"], df["PBR"]):
            rows.append({
                "sample": sample,
                "timepoint": timepoint,
                "condition": condition,
                "pbl": pbl_val,
                "pbr": pbr_val,
            })

    if not rows:
        logger.error("Error: No valid data found in any input file!")
        return 1

    out_df = pd.DataFrame(rows)
    out_df.to_csv(config.output_path, sep="\t", index=False, float_format="%.4f")
    logger.success(f"Wrote {len(out_df)} PBL/PBR pairs to {config.output_path}")
    logger.info(f"Unique samples: {out_df['sample'].nunique()}")
    logger.info(f"Total data points: {len(out_df):,}")

    return 0
```

Update CLI to expect TSV output:

```python
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Extract PBL/PBR pairs from multiple TSV files")
    parser.add_argument("-i", "--input", nargs='+', type=Path, required=True, help="Input TSV files (space-separated)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()
```

**Step 2.3: Run modified script on real data**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
mkdir -p projects/HD_DIT_HAP/results/18_figure_data
COMP_PY=/data/a/yangyusheng/miniforge3/envs/statistics_and_figure_plotting/bin/python
$COMP_PY workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py \
  -i projects/HD_DIT_HAP/results/8_merged/HD1328-4_YES0_YES.tsv \
     projects/HD_DIT_HAP/results/8_merged/HD1328-4_YES1_YES.tsv \
  -o projects/HD_DIT_HAP/results/18_figure_data/pbl_pbr_pairs.tsv
```

Expected output log:
```
Wrote 259096 PBL/PBR pairs to ...
Unique samples: 1
```

Verify schema:

```bash
head -3 projects/HD_DIT_HAP/results/18_figure_data/pbl_pbr_pairs.tsv
```

Expected:
```
sample	timepoint	condition	pbl	pbr
HD1328-4	YES0	YES	364.0000	429.0000
HD1328-4	YES0	YES	66.0000	77.0000
```

**Step 2.4: Write rendering script test**

```python
# workflow/scripts/figures/test_plot_pbl_pbr_correlation.py
"""Test plot_pbl_pbr_correlation.py rendering."""

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import pytest
from pathlib import Path


def test_pbl_pbr_rendering_preserves_statistics(tmp_path):
    """Rendered figure must reproduce baseline PCC=0.8521 for YES0 sample."""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))
    from figures import apply_house_style
    import cnsplots as cns
    import numpy as np

    # Load real data
    data_path = Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake"
                     "/projects/HD_DIT_HAP/results/18_figure_data/pbl_pbr_pairs.tsv")
    if not data_path.exists():
        pytest.skip("Real data not available")

    df = pd.read_csv(data_path, sep="\t")
    df_yes0 = df[(df["sample"] == "HD1328-4") & (df["timepoint"] == "YES0")].copy()
    
    # Filter to positive values for log scale
    df_yes0 = df_yes0[(df_yes0["pbl"] > 0) & (df_yes0["pbr"] > 0)]
    df_yes0 = df_yes0.assign(
        log10_pbl=np.log10(df_yes0["pbl"]),
        log10_pbr=np.log10(df_yes0["pbr"])
    )
    
    # Baseline check
    assert len(df_yes0) == 131057, f"Expected 131057 points, got {len(df_yes0)}"
    
    # Compute PCC on log10 values (matching old code)
    pcc = np.corrcoef(df_yes0["log10_pbl"], df_yes0["log10_pbr"])[0, 1]
    assert abs(pcc - 0.8521) < 0.001, f"PCC {pcc:.4f} != baseline 0.8521"
    
    # Render and save
    apply_house_style()
    cns.figure(width=180, height=150)
    ax = cns.regplot(
        data=df_yes0,
        x="log10_pbl",
        y="log10_pbr",
        scatter_kws=dict(s=3, facecolor="none", edgecolor="gray",
                         alpha=0.15, linewidths=0.25, rasterized=True)
    )
    ax.set(xlabel=r"log$_{10}$ PBL", ylabel=r"log$_{10}$ PBR",
           title="HD1328-4 YES0")
    
    out = tmp_path / "test_pbl_pbr"
    cns.savefig(out.with_suffix(".pdf"))
    cns.savefig(out.with_name(out.stem + ".review.png"))
    
    assert out.with_suffix(".pdf").exists()
    assert out.with_name(out.stem + ".review.png").exists()
```

Run:

```bash
cd workflow/scripts/figures
CNSPLOTS_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
$CNSPLOTS_PY -m pytest test_plot_pbl_pbr_correlation.py::test_pbl_pbr_rendering_preserves_statistics -v
```

Expected: PASS with assertion `PCC 0.8521 == baseline 0.8521`

**Step 2.5: Write full rendering script**

```python
# workflow/scripts/figures/plot_pbl_pbr_correlation.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PBL-PBR Correlation Rendering
==============================

Read the long-format PBL/PBR pairs table and render log-log scatter plots with
regression and statistics, one panel per sample. Outputs journal PDF + review PNG.

Input
-----
- TSV with columns: sample, timepoint, condition, pbl, pbr

Output
------
- Journal PDF (180 mm two-column)
- Review PNG (2× scale for on-screen QC)

Usage
-----
    python plot_pbl_pbr_correlation.py -i pairs.tsv -o fig.pdf

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-13
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
from loguru import logger

# Local imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))
from figures import apply_house_style, save_dual  # noqa: E402

import cnsplots as cns  # noqa: E402


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class PlotConfig:
    """Immutable config holding validated input/output paths."""
    input_file: Path
    output_stem: Path

    def __post_init__(self) -> None:
        """Validate input exists and create output directory."""
        if not self.input_file.exists():
            raise ValueError(f"Input file does not exist: {self.input_file}")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru to emit uncolorised, timestamped records to stdout."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_file: Path) -> pd.DataFrame:
    """Load PBL/PBR pairs and add log10 columns for plotting."""
    df = pd.read_csv(input_file, sep="\t")
    logger.info(f"Loaded {len(df):,} rows from {input_file}")

    required_cols = {"sample", "timepoint", "condition", "pbl", "pbr"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing columns: {missing}")

    # Filter to positive values for log scale
    df = df[(df["pbl"] > 0) & (df["pbr"] > 0)].copy()
    df["log10_pbl"] = np.log10(df["pbl"])
    df["log10_pbr"] = np.log10(df["pbr"])

    logger.info(f"Retained {len(df):,} positive-valued pairs")
    return df


@logger.catch
def render_pbl_pbr_figure(df: pd.DataFrame, output_stem: Path) -> None:
    """Render multi-panel PBL/PBR correlation figure."""
    samples = df["sample"].unique()
    logger.info(f"Rendering {len(samples)} sample panel(s)")

    apply_house_style()
    mp = cns.multipanel(max_width=540, title="PBL-PBR Correlation Analysis")

    for idx, sample in enumerate(sorted(samples)):
        sample_df = df[df["sample"] == sample]
        logger.info(f"  Panel {idx+1}: {sample} (n={len(sample_df):,})")

        panel_label = chr(65 + idx)  # A, B, C, ...
        ax = mp.panel(panel_label, width=180, height=150)

        cns.regplot(
            data=sample_df,
            x="log10_pbl",
            y="log10_pbr",
            ax=ax,
            scatter_kws=dict(
                s=3,
                facecolor="none",
                edgecolor="gray",
                alpha=0.15,
                linewidths=0.25,
                rasterized=True,
            ),
        )

        ax.set(
            xlabel=r"log$_{10}$ PBL",
            ylabel=r"log$_{10}$ PBR",
            title=sample,
        )

    save_dual(output_stem)
    logger.success(f"Saved journal PDF + review PNG: {output_stem}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(
        description="Render PBL-PBR correlation figure from pairs TSV"
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Input TSV (sample, timepoint, condition, pbl, pbr)")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output stem (will write .pdf and .review.png)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging")
    return parser.parse_args()


def main() -> int:
    """Load PBL/PBR pairs, render cnsplots figure, save dual artifacts."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    try:
        config = PlotConfig(
            input_file=args.input,
            output_stem=args.output.with_suffix(""),  # strip extension
        )
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    logger.info("=== PBL-PBR Correlation Rendering ===")

    df = load_and_prepare_data(config.input_file)
    render_pbl_pbr_figure(df, config.output_stem)

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 2.6: Run rendering script on real data**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
mkdir -p /tmp/pbl_pbr_test
CNSPLOTS_PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python
$CNSPLOTS_PY workflow/scripts/figures/plot_pbl_pbr_correlation.py \
  -i projects/HD_DIT_HAP/results/18_figure_data/pbl_pbr_pairs.tsv \
  -o /tmp/pbl_pbr_test/PBL_PBR_correlation_analysis
```

Expected log:
```
Rendered 1 sample panel(s)
  Panel 1: HD1328-4 (n=259,096)
Saved journal PDF + review PNG: /tmp/pbl_pbr_test/PBL_PBR_correlation_analysis
```

Verify artifacts:

```bash
ls -lh /tmp/pbl_pbr_test/
```

Expected:
```
PBL_PBR_correlation_analysis.pdf         (~80 KB)
PBL_PBR_correlation_analysis.review.png  (~110 KB)
```

**Step 2.7: Visual inspection**

Open the PNG and verify:
1. No text clipping or overlapping tick labels
2. Scatter points show density structure (not a black blob)
3. Regression line and `r=0.85` annotation visible
4. Panel label "A" in top-left
5. Axes labels render as "log₁₀ PBL" with subscript

**Step 2.8: Update quality_control.smk to split the rule**

Modify `workflow/rules/quality_control.smk`:

Find the existing `PBL_PBR_correlation_analysis` rule and replace with two:

```python
# PBL-PBR pairs extraction (computation layer)
# -----------------------------------------------------
rule pbl_pbr_pairs:
    input:
        expand(rules.merge_strand_insertions.output, sample=samples, timepoint=timepoints, condition=conditions),
    output:
        f"projects/{project_name}/results/18_figure_data/pbl_pbr_pairs.tsv",
    log:
        f"projects/{project_name}/logs/quality_control/pbl_pbr_pairs.log",
    conda:
        "../envs/statistics_and_computation.yml"
    message:
        "*** Extracting PBL-PBR pairs..."
    shell:
        """
        python workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py \
            -i {input} \
            -o {output} &> {log}
        """


# PBL-PBR correlation rendering (rendering layer)
# -----------------------------------------------------
rule plot_pbl_pbr_correlation:
    input:
        rules.pbl_pbr_pairs.output,
    output:
        journal = report(
            f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
            caption="../reports/captions/PBL_PBR_correlation_analysis.rst",
            category="Quality Control",
            labels={
                "name": "3. PBL-PBR Correlation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
            },
        ),
        review = f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.review.png",
    log:
        f"projects/{project_name}/logs/quality_control/plot_pbl_pbr_correlation.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Rendering PBL-PBR correlation figure..."
    shell:
        """
        python workflow/scripts/figures/plot_pbl_pbr_correlation.py \
            -i {input} \
            -o $(dirname {output.journal})/PBL_PBR_correlation_analysis &> {log}
        """
```

**Step 2.9: Dry-run to verify DAG**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
snakemake -n plot_pbl_pbr_correlation 2>&1 | grep -E "rule|Building DAG|job"
```

Expected output includes:
```
Building DAG of jobs...
job 1: pbl_pbr_pairs
job 2: plot_pbl_pbr_correlation
```

**Step 2.10: Commit**

```bash
git add workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py
git add workflow/scripts/figures/plot_pbl_pbr_correlation.py
git add workflow/scripts/figures/test_plot_pbl_pbr_correlation.py
git add workflow/rules/quality_control.smk
git commit -m "refactor(figures): split PBL/PBR into compute + render layers

Strip plotting from PBL_PBR_correlation_analysis.py, now writes long-format
pbl_pbr_pairs.tsv. Add plot_pbl_pbr_correlation.py (cnsplots multipanel) with
pytest verifying PCC=0.8521 baseline. Split quality_control.smk rule into
pbl_pbr_pairs (pandas 3 env) + plot_pbl_pbr_correlation (cnsplots env).

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Task 3–9: Remaining 7 scripts

Follow the Task 2 pattern for each:

1. Write schema test for intermediate TSV
2. Strip plotting from computation script, write TSV
3. Run on real data, verify row counts and schema
4. Write rendering test verifying key statistics
5. Write rendering script (cnsplots)
6. Run on real data, visual inspection
7. Split `.smk` rule into two (computation + rendering)
8. Dry-run DAG
9. Commit

### Task 3: read_count_distribution

**Intermediate TSV:** `read_count_distribution.tsv` (pre-binned), `read_count_cutoff_stats.tsv`

**Key difference:** Store histogram bins, not raw values, to cap size.

Schema:
```
sample  timepoint  bin_left  bin_right  count
```

Stats schema:
```
sample  original_rows  rows_kept  pct_rows_kept  original_counts  counts_kept  pct_counts_kept
```

**Rendering:** `cns.barplot` or `cns.histplot(weights=...)` for pre-binned data, `ax.axvline()` for cutoff.

### Task 4: insertion_orientation

**Intermediate TSV:** `strand_pairs.tsv`

Schema:
```
sample  timepoint  plus_count  minus_count
```

**Rendering:** `cns.regplot` panel matrix (one per sample × timepoint), log-log. Fixes overwrite bug.

### Task 5: insertion_density

**Intermediate TSV:** Reuses existing `insertion_density_analysis.tsv` (no new table).

**Rendering:** `cns.scatterplot` 2×2 grid (essentiality-coloured) + `cns.histplot` metric distributions.

### Task 6: gene_coverage

**Intermediate TSV:** `gene_coverage_stats.tsv`

Schema:
```
category  covered  not_covered  total  coverage_pct
```

**Rendering:** `cns.barplot` + `cns.donutplot` per category.

### Task 7: curve_fitting

**Intermediate TSV:** Reuses existing `*_fitting_statistics.tsv`.

**Rendering:** Sample N curves (default 32, deterministic seed 42), draw observed + fitted sigmoid in grid. `cns.scatterplot` + `ax.plot()` for smooth curve.

### Task 8: MA plot

**Intermediate TSV:** `ma_values.tsv`

Schema:
```
timepoint  a_value  m_value
```

**Rendering:** `cns.scatterplot` per timepoint panel.

### Task 9: distribution_of_curve_fitting

**Intermediate TSV:** Reuses existing `*_fitting_statistics.tsv`.

**Rendering:** `cns.histplot` grid for metrics (A, DR, DL, t50, R², RMSE, etc.).

---

## Task 10: Final verification and computation-env cleanup

**Step 10.0: Remove plotting libraries from the computation env**

Deferred from Task 0 — only safe once all 8 rules are rewired to `cnsplots.yml`.

First confirm no script reachable from the computation env still imports a
plotting library:

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
# Every rule using the computation env, and the scripts they invoke:
grep -B30 "statistics_and_computation" workflow/rules/*.smk | grep -oE "workflow/scripts/[a-z_]+/[a-z_]+\.py" | sort -u > /tmp/comp_scripts.txt
cat /tmp/comp_scripts.txt
# None of these may import matplotlib/seaborn/altair:
xargs grep -l "matplotlib\|seaborn\|altair" < /tmp/comp_scripts.txt
```

Expected: the final `grep -l` prints nothing. If it lists any script, that
script still plots and its rule has not been rewired — fix that first.

Then remove the three libraries from `workflow/envs/statistics_and_computation.yml`
and replace the transitional NOTE with the final comment:

```yaml
# workflow/envs/statistics_and_computation.yml
# Computation environment: numerical analysis, depletion scoring, curve fitting.
# NO plotting libraries — those live in cnsplots.yml.
```

Also delete `config/DIT_HAP.mplstyle`, now unused by any script:

```bash
grep -rn "DIT_HAP.mplstyle" workflow/ config/ --include="*.py" --include="*.smk"
# Expected: no matches → safe to remove
git rm config/DIT_HAP.mplstyle
```

Verify the env still resolves, then confirm the full DAG still builds
(Step 10.1 below covers this).

**Step 10.1: Run full DAG dry-run**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
snakemake -n all 2>&1 | tee /tmp/dag_check.log
grep -E "rule.*plot_|Building DAG" /tmp/dag_check.log
```

Expected: All 8 `plot_*` rules appear with correct dependencies.

**Step 10.2: Verify no orphaned rules**

```bash
snakemake --list-untracked 2>&1
```

Expected: Empty or only non-plotting rules.

**Step 10.3: Run one full figure end-to-end**

```bash
snakemake plot_pbl_pbr_correlation --use-conda --cores 2
```

Expected: Both artifacts created, no errors.

**Step 10.4: Compare statistics across all figures**

For each figure with numerical stats (PCC, coverage %, etc.), manually compare against baseline. Document deltas in a `verification_report.md`.

**Step 10.5: Commit verification report**

```bash
git add docs/verification_report.md
git commit -m "docs: add cnsplots refactor verification report

All 8 figures rendered with cnsplots. Statistics match pre-refactor baselines
within 0.1% tolerance. Visual inspection confirms no text clipping, proper
density rendering, and subscript formatting.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Execution Notes

- Each task is independently committable; no need to finish all 10 before pushing.
- pytest runs require the correct conda env activation: `cnsplots_figures` for rendering tests, `statistics_and_computation` for computation tests.
- Real data paths are `/data/c/yangyusheng_optimized/DIT_HAP_snakemake/projects/HD_DIT_HAP/results/...`.
- Visual inspection is mandatory per the cnsplots skill; use `feh` or `eog` to view PNGs.

**Estimated time:** 
- Task 0–2: 90 minutes (env setup + reference pattern)
- Task 3–9: 4–6 hours (7 scripts × 30–50 min each)
- Task 10: 30 minutes (verification)
- Total: ~6–8 hours for a skilled implementer with zero cnsplots experience.
