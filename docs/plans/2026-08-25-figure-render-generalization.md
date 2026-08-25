# Figure Render Generalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize `workflow/src/figure_render/` from ten figure-specific modules into six general plot-grammar modules, moving all figure-specific column names and labels into the `workflow/scripts/figures/` entrypoints.

**Architecture:** Each `src/figure_render/` module owns one *kind* of plot and receives every figure-specific value as a keyword argument. Two shared helpers (`_layout.py`, `_schema.py`) absorb logic currently duplicated across five and nine modules respectively. Loaders move to the entrypoint scripts, because a required-column list is figure-specific knowledge. Correctness is protected by an exact-pixel baseline captured before any code changes.

**Tech Stack:** Python 3.12, pandas 2.3, cnsplots 0.6.0, matplotlib 3.10.9 (render env), numpy 2.4.6, PIL 12.3.0, pytest 9.1.1, Snakemake.

**Spec:** `docs/plans/2026-08-25-figure-render-generalization-design.md` (commit `99bad30`)

## Global Constraints

- **Render env python:** `/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python`. All test and render commands in this plan use it explicitly. Do not use bare `python`.
- **Repo root:** `/data/c/yangyusheng_optimized/DIT_HAP_snakemake`. All relative paths are from here.
- **`src/` modules are library modules:** IMPORTS → CONSTANTS → CORE LOGIC only. No `main()`, no `parse_args()`, no `setup_logger()`. Per the `python-script-conventions` skill.
- **No module in `src/figure_render/` may name a biological quantity** (`pbl`, `pbr`, `plus_count`, `minus_count`, `FYPOviability`, `baseMean`) **or a timepoint literal** (`YES0`). This is the acceptance test for the whole plan.
- **`x`, `y`, `xlabel`, `ylabel` are keyword-only with no default** on every renderer, so a stale positional call raises `TypeError` instead of silently drawing wrong labels.
- **cnsplots 0.6.0 gotchas that silently corrupt output:**
  - `regplot` resets axis labels it manages → apply labels *after* drawing.
  - `regplot` drops its own `s=` when `scatter_kws` is supplied → marker size goes *inside* `scatter_kws`.
  - `regplot` computes `r`/`P` from the columns given → callers pass explicit `log10_*` columns (project convention is PCC on log10; raw-space r ≈ -0.03 vs log-space r = 0.85 for PBL/PBR).
  - `cns.scatterplot` always passes `edgecolor=None` internally → never supply `edgecolor` to it.
  - Sizes in `cns.figure()` / `mp.panel()` are **pixels** (`inches = px / 72`).
  - `cns.multipanel` is a class, not a function.
- **Zero `findfont` warnings** is a pass criterion. `apply_house_style()` in `workflow/src/figures.py` sets Arial; Helvetica is absent on this machine.
- **House style palette** (Cell, from `apply_house_style()`): `['#c84c3a', '#2f7e8f', '#e1a22e', '#4e5a8a', '#5f9862', '#d07a6a', '#8b6fa8', '#7b8c9e', '#b85f7a', '#6b6b6b']`.
- **Existing loader convention:** loaders are decorated `@logger.catch`, so a `ValueError` is logged and the function returns `None` rather than propagating. Three existing tests assert `is None`. Preserve this.
- **No Snakefile / rule / env changes.** Rules keep invoking `python workflow/scripts/figures/plot_*.py` with identical CLI flags.

---

## File Structure

**Created in `workflow/src/figure_render/`:**

| File | Responsibility |
|---|---|
| `_layout.py` | `panel_labels()`, `grid_panel_size()`, `PANEL_DECORATION_PX`. No plotting. |
| `_schema.py` | `require_columns()`. No plotting, no I/O. |
| `scatter.py` | `render_grouped_regression_figure()` (grouped log-log regplot grid), `render_scatter_grid_figure()` (fixed multi-panel scatter with hue) |
| `histogram.py` | `render_histogram_grid_figure()` — raw-values mode and pre-binned-weights mode |
| `series.py` | `render_series_scatter_figure()` — N series against one shared x on log-log axes |
| `ma.py` | `render_ma_figure()` — wide (per-timepoint columns) and long (per-row timepoint) inputs, optional per-point colours, vertical/horizontal orientation |
| `curves.py` | `render_fitted_curves_figure()` — observed points + fitted model overlay grid |
| `composition.py` | `render_composition_figure()` — bar overview + per-category donuts |

**Deleted from `workflow/src/figure_render/`** (all ten): `correlation.py`, `orientation.py`, `density.py`, `distribution.py`, `read_counts.py`, `ma_plot.py`, `ma_plot_replicates.py`, `dispersions.py`, `curve_fitting.py`, `coverage.py`.

**Modified — the ten entrypoints** in `workflow/scripts/figures/` each gain a CONSTANTS block (column names, axis labels, required-column list) and a loader function moved down from `src`:

`plot_pbl_pbr_correlation.py`, `plot_insertion_orientation.py`, `plot_insertion_density.py`, `plot_distribution_of_curve_fitting.py`, `plot_read_count_distribution.py`, `plot_ma_plot.py`, `plot_ma_plot_replicates.py`, `plot_dispersions.py`, `plot_curve_fitting.py`, `plot_gene_coverage.py`

**Test files:** ten under `workflow/tests/figure_render/` are rewritten to import from the new grammar modules; two new files `test_layout.py` and `test_schema.py`; one new `test_pixel_baseline.py`.

**Why loaders move to scripts:** every loader is `read_csv` + required-column validation + a few derived columns. The required-column list names biological quantities, so under the spec's boundary it belongs in `scripts/`. `src` keeps only the generic `require_columns()` primitive.

**Left untouched:** `workflow/tests/scripts/test_plot_configs.py` and `workflow/tests/scripts/conftest.py`. That test parametrizes over five scripts' `PlotConfig` dataclasses, asserting each rejects a non-existent input path. No task modifies a `PlotConfig` — only imports, a new CONSTANTS block, a loader function, and the `main()` render call change — so it must keep passing untouched. Task 15's full-suite run covers it. Its `load_script` fixture is the same `importlib.util` pattern the rewritten figure tests adopt, because `workflow/scripts/figures/` has no `__init__.py` and `workflow/src/figures.py` shadows the directory name on `sys.path`.

---

## Phase 1 — Shared helpers and pixel baseline

### Task 1: Capture the pixel baseline

This must happen **before any code change** or the baseline is contaminated. It produces no source changes — only a stored set of reference PNGs and a harness script.

**Files:**
- Create: `workflow/tests/figure_render/pixel_baseline.py` (harness, importable by later tasks)
- Create: `tmp/pixel_baseline/*.png` (git-ignored artifacts)

**Interfaces:**
- Consumes: nothing (first task)
- Produces: `render_all_baseline_figures(out_dir: Path) -> dict[str, Path]` mapping figure name → PNG path; `compare_png(a: Path, b: Path) -> tuple[bool, int]` returning `(identical, max_abs_diff)`

- [ ] **Step 1: Confirm PNG rendering is deterministic**

The whole strategy depends on this. Run:

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src')
from matplotlib import use; use('Agg')
from pathlib import Path
import numpy as np
from PIL import Image
from figure_render.coverage import load_coverage_data, render_coverage_figure
df = load_coverage_data(Path('projects/HD_DIT_HAP/results/18_figure_data/arc/gene_coverage_stats.tsv'))
frames = []
for i in (1, 2):
    render_coverage_figure(df, Path(f'/tmp/_det{i}'))
    frames.append(np.asarray(Image.open(f'/tmp/_det{i}.review.png').convert('RGB')))
print('identical:', np.array_equal(frames[0], frames[1]))
"
```

Expected: `identical: True`. If this prints `False`, STOP — pixel comparison is not a valid strategy and the plan needs revisiting.

- [ ] **Step 2: Write the baseline harness**

Create `workflow/tests/figure_render/pixel_baseline.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Pixel baseline harness for figure render refactoring.

Renders every figure from ``HD_DIT_HAP`` data and compares PNGs byte-exactly.
PNG output was verified deterministic (max abs diff 0 across repeat renders),
so comparison needs no tolerance. PDFs embed timestamps and are not compared.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
from PIL import Image

# =============================================================================
# CONSTANTS
# =============================================================================
PROJECT_ROOT = Path("/data/c/yangyusheng_optimized/DIT_HAP_snakemake")
ARC_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/18_figure_data/arc"
DEPLETION_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/14_insertion_level_depletion_analysis"
FITTING_DIR = PROJECT_ROOT / "projects/HD_DIT_HAP/results/15_insertion_level_curve_fitting"

# Figures whose pixels are expected to change by design (spec section
# "Expected pixel changes"): correlation adopts orientation's grid layout,
# curve_fitting gets the house-style palette.
EXPECTED_TO_CHANGE = frozenset({"correlation", "curve_fitting"})


# =============================================================================
# CORE LOGIC
# =============================================================================
def compare_png(baseline: Path, candidate: Path) -> tuple[bool, int]:
    """Compare two PNGs pixel-exactly, returning (identical, max_abs_diff)."""
    a = np.asarray(Image.open(baseline).convert("RGB")).astype(int)
    b = np.asarray(Image.open(candidate).convert("RGB")).astype(int)

    if a.shape != b.shape:
        return False, -1

    diff = int(np.abs(a - b).max())
    return diff == 0, diff
```

- [ ] **Step 3: Add the render driver to the harness**

Append to `workflow/tests/figure_render/pixel_baseline.py`:

```python
def render_all_baseline_figures(out_dir: Path) -> dict[str, Path]:
    """Render all ten figures into out_dir, returning {figure_name: png_path}.

    Imports are function-local so this module stays importable while modules are
    mid-refactor: a figure whose module has been replaced raises ImportError,
    which is recorded rather than aborting the whole sweep.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rendered: dict[str, Path] = {}

    def _run(name: str, fn) -> None:
        stem = out_dir / name
        try:
            fn(stem)
        except ImportError as exc:
            print(f"SKIP {name}: {exc}")
            return
        png = stem.parent / f"{stem.name}.review.png"
        if png.exists():
            rendered[name] = png
        else:
            print(f"MISSING {name}: no PNG produced")

    _run("correlation", _render_correlation)
    _run("orientation", _render_orientation)
    _run("read_counts", _render_read_counts)
    _run("density", _render_density)
    _run("coverage", _render_coverage)
    _run("dispersions", _render_dispersions)
    _run("ma_plot_replicates", _render_ma_replicates)
    _run("ma_plot", _render_ma_plot)
    _run("distribution", _render_distribution)
    _run("curve_fitting", _render_curve_fitting)

    return rendered
```

- [ ] **Step 4: Add the ten per-figure render shims**

Append to `workflow/tests/figure_render/pixel_baseline.py`. These call whatever API currently exists; Task 15 updates them to the new API.

```python
def _render_correlation(stem: Path) -> None:
    from figure_render.correlation import load_and_prepare_data, render_correlation_figure
    df = load_and_prepare_data(ARC_DIR / "pbl_pbr_pairs.tsv")
    render_correlation_figure(df, stem)


def _render_orientation(stem: Path) -> None:
    from figure_render.orientation import load_and_prepare_data, render_orientation_figure
    df = load_and_prepare_data(ARC_DIR / "strand_pairs.tsv")
    render_orientation_figure(df, stem)


def _render_read_counts(stem: Path) -> None:
    from figure_render.read_counts import (
        load_cutoff_stats, load_distribution_data, render_distribution_figure,
    )
    df = load_distribution_data(ARC_DIR / "read_count_distribution.tsv")
    stats = load_cutoff_stats(ARC_DIR / "read_count_cutoff_stats.tsv")
    render_distribution_figure(df, stats, stem, "YES0", 8.0)


def _render_density(stem: Path) -> None:
    from figure_render.density import load_density_data, render_density_figure
    df = load_density_data(ARC_DIR / "insertion_density_analysis.tsv")
    render_density_figure(df, stem, "YES0", "YES4")


def _render_coverage(stem: Path) -> None:
    from figure_render.coverage import load_coverage_data, render_coverage_figure
    df = load_coverage_data(ARC_DIR / "gene_coverage_stats.tsv")
    render_coverage_figure(df, stem)


def _render_dispersions(stem: Path) -> None:
    from figure_render.dispersions import load_dispersion_data, render_dispersion_figure
    df = load_dispersion_data(ARC_DIR / "dispersion_data.tsv")
    render_dispersion_figure(df, stem)


def _render_ma_replicates(stem: Path) -> None:
    from figure_render.ma_plot_replicates import load_ma_data, render_ma_figure
    df = load_ma_data(ARC_DIR / "ma_values.tsv")
    render_ma_figure(df, stem)


def _render_ma_plot(stem: Path) -> None:
    from figure_render.ma_plot import Orientation, load_ma_data, render_ma_figure
    basemean, lfc = load_ma_data(DEPLETION_DIR / "baseMean.tsv", DEPLETION_DIR / "LFC.tsv")
    render_ma_figure(basemean, lfc, stem, Orientation.VERTICAL)


def _render_distribution(stem: Path) -> None:
    from figure_render.distribution import load_fitting_stats, render_distribution_figure
    df, metric_cols = load_fitting_stats(FITTING_DIR / "insertion_level_fitting_statistics.tsv")
    render_distribution_figure(df, metric_cols, stem, 30)


def _render_curve_fitting(stem: Path) -> None:
    from figure_render.curve_fitting import load_and_sample_data, render_curve_fitting_figure
    stats, lfc, time_points = load_and_sample_data(
        FITTING_DIR / "insertion_level_fitting_statistics.tsv",
        DEPLETION_DIR / "LFC.tsv",
        32,
        42,
    )
    render_curve_fitting_figure(stats, lfc, time_points, stem)
```

- [ ] **Step 5: Capture the baseline**

Run:

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
from pixel_baseline import render_all_baseline_figures
got = render_all_baseline_figures(Path('tmp/pixel_baseline'))
print('captured', len(got), 'figures:', sorted(got))
" 2>&1 | grep -v "INFO\|SUCCESS" | tail -20
```

Expected: `captured 10 figures:` with all ten names. Record which (if any) fail — `ma_plot` reads `baseMean.tsv`/`LFC.tsv` and `curve_fitting` is known to work for `HD_DIT_HAP`. If a figure fails to render *before* any changes, note it and treat "still fails identically" as its pass criterion.

- [ ] **Step 6: Verify no findfont warnings in the baseline run**

Run the Step 5 command again, filtering for the warning:

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
from pixel_baseline import render_all_baseline_figures
render_all_baseline_figures(Path('tmp/_findfont_check'))
" 2>&1 | grep -c findfont
```

Expected: `0`.

- [ ] **Step 7: Ignore the baseline artifacts and commit the harness**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
grep -qxF 'tmp/' .gitignore || echo 'tmp/' >> .gitignore
git add .gitignore workflow/tests/figure_render/pixel_baseline.py
git commit -m "test(figures): add pixel baseline harness for render refactor

PNG output verified deterministic (max abs diff 0 across repeat renders),
so comparison is byte-exact with no tolerance."
```

---

### Task 2: `_layout.py` — panel labels and grid sizing

**Files:**
- Create: `workflow/src/figure_render/_layout.py`
- Test: `workflow/tests/figure_render/test_layout.py`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `PANEL_DECORATION_PX: int = 40`
  - `panel_labels(n: int) -> list[str]` — `["A", ..., "Z", "A1", "A2", ...]`
  - `grid_panel_size(total_width_px, total_height_px, n_cols, n_rows, decoration_px=PANEL_DECORATION_PX, min_width=45, min_height=55) -> tuple[int, int]`

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_layout.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for shared figure layout helpers.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pytest

from figure_render._layout import PANEL_DECORATION_PX, grid_panel_size, panel_labels


# =============================================================================
# TESTS
# =============================================================================
def test_panel_labels_are_letters_below_26() -> None:
    """Assert the first 26 labels are plain uppercase letters."""
    labels = panel_labels(26)

    assert labels[0] == "A"
    assert labels[25] == "Z"
    assert all(label.isalpha() and label.isupper() for label in labels)


def test_panel_labels_do_not_emit_punctuation_past_z() -> None:
    """Assert overflow past Z stays alphanumeric, guarding the chr(65+i) bug.

    The old bare ``chr(65 + i)`` produced '[', '\\\\' and ']' for i >= 26.
    curve_fitting.py renders 32 panels by default, so it was already emitting
    those characters as panel labels.
    """
    labels = panel_labels(32)

    assert len(labels) == 32
    for label in labels:
        assert label.isalnum(), f"Non-alphanumeric panel label: {label!r}"
    assert "[" not in labels
    assert "\\" not in labels
    assert "]" not in labels


def test_panel_labels_are_unique() -> None:
    """Assert no label repeats, so panels stay individually citable."""
    labels = panel_labels(60)

    assert len(set(labels)) == 60


def test_panel_labels_zero_returns_empty() -> None:
    """Assert a zero-panel figure yields no labels rather than raising."""
    assert panel_labels(0) == []


def test_grid_panel_size_subtracts_decoration() -> None:
    """Assert per-panel size is the share of the figure minus decoration."""
    width, height = grid_panel_size(510, 425, n_cols=5, n_rows=3, decoration_px=40)

    assert width == 510 // 5 - 40
    assert height == 425 // 3 - 40


def test_grid_panel_size_enforces_floor() -> None:
    """Assert crowded grids clamp to the minimum rather than going negative."""
    width, height = grid_panel_size(510, 425, n_cols=50, n_rows=50, decoration_px=40)

    assert width == 45
    assert height == 55


def test_grid_panel_size_rejects_zero_cols() -> None:
    """Assert a zero column count is refused instead of dividing by zero."""
    with pytest.raises(ValueError, match="n_cols"):
        grid_panel_size(510, 425, n_cols=0, n_rows=1)


def test_grid_panel_size_reproduces_orientation_values() -> None:
    """Assert the helper matches orientation.py's arithmetic for its real grid.

    orientation.py rendered a 5-timepoint by 3-sample grid at
    JOURNAL_WIDTH_PX=510 / JOURNAL_HEIGHT_PX=425 with PANEL_DECORATION_PX=40.
    Reproducing it exactly is what keeps that figure pixel-identical.
    """
    width, height = grid_panel_size(510, 425, n_cols=5, n_rows=3, decoration_px=PANEL_DECORATION_PX)

    assert width == max(45, int(510 / 5) - 40)
    assert height == max(55, int(425 / 3) - 40)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_layout.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render._layout'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/_layout.py`:

```python
"""Shared panel layout helpers for figure rendering.

Panel-letter generation and grid sizing were duplicated across five figure
modules in two incompatible variants; this module is the single implementation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================

# =============================================================================
# CONSTANTS
# =============================================================================
# Pixels each panel spends on its label, tick labels and axis title, on top of
# the axes box itself. Subtracted from the available width/height so a whole row
# of panels fits inside the journal width. Measured at 40 px for log10 scatter
# panels; frequency histograms carry wider decorations and pass 55 explicitly.
# Overestimating shrinks the axes enough for an extra panel to fit in a row,
# which silently reflows the grid.
PANEL_DECORATION_PX = 40

# Uppercase A-Z, then A1, A2, ... A bare chr(65 + i) emits '[', '\', ']' past Z.
_ALPHABET_SIZE = 26


# =============================================================================
# CORE LOGIC
# =============================================================================
def panel_labels(n: int) -> list[str]:
    """Return n unique panel labels: A-Z, then A1, A2, ... for n > 26."""
    return [
        chr(65 + i) if i < _ALPHABET_SIZE else f"A{i - _ALPHABET_SIZE + 1}"
        for i in range(n)
    ]


def grid_panel_size(
    total_width_px: int,
    total_height_px: int,
    n_cols: int,
    n_rows: int,
    decoration_px: int = PANEL_DECORATION_PX,
    min_width: int = 45,
    min_height: int = 55,
) -> tuple[int, int]:
    """Return the (width, height) in pixels for one panel of an n_cols x n_rows grid."""
    if n_cols < 1:
        raise ValueError(f"n_cols must be at least 1, got {n_cols}")
    if n_rows < 1:
        raise ValueError(f"n_rows must be at least 1, got {n_rows}")

    width = max(min_width, int(total_width_px / n_cols) - decoration_px)
    height = max(min_height, int(total_height_px / n_rows) - decoration_px)

    return width, height
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_layout.py -v
```

Expected: 8 passed.

- [ ] **Step 5: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/_layout.py workflow/tests/figure_render/test_layout.py
git commit -m "feat(figure_render): add shared panel layout helpers

panel_labels() replaces five duplicated implementations, two of which used a
bare chr(65+i) that emitted '[', '\\\\', ']' past 26 panels. curve_fitting
renders 32 panels by default and was already affected."
```

---

### Task 3: `_schema.py` — required-column validation

**Files:**
- Create: `workflow/src/figure_render/_schema.py`
- Test: `workflow/tests/figure_render/test_schema.py`

**Interfaces:**
- Consumes: nothing
- Produces: `require_columns(df: pd.DataFrame, required: Sequence[str], context: str = "input") -> None` — raises `ValueError` listing every missing column

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_schema.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for shared figure input schema validation.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
import pandas as pd
import pytest

from figure_render._schema import require_columns


# =============================================================================
# TESTS
# =============================================================================
def test_passes_when_all_columns_present() -> None:
    """Assert a frame holding every required column validates silently."""
    df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})

    require_columns(df, ["a", "b"])


def test_passes_on_empty_frame_with_correct_columns() -> None:
    """Assert an empty frame with the right schema is accepted.

    Empty-input handling is exercised by every figure's empty-data test, which
    writes a header-only TSV; those must not fail validation.
    """
    df = pd.DataFrame(columns=["a", "b"])

    require_columns(df, ["a", "b"])


def test_raises_listing_every_missing_column() -> None:
    """Assert the error names all missing columns, not just the first."""
    df = pd.DataFrame({"a": [1]})

    with pytest.raises(ValueError) as excinfo:
        require_columns(df, ["a", "b", "c"])

    message = str(excinfo.value)
    assert "b" in message
    assert "c" in message


def test_error_message_includes_context() -> None:
    """Assert the caller-supplied context appears, so the failing file is identifiable."""
    df = pd.DataFrame({"a": [1]})

    with pytest.raises(ValueError, match="strand pairs"):
        require_columns(df, ["plus"], context="strand pairs TSV")


def test_preserves_required_order_in_message() -> None:
    """Assert missing columns are reported in the order the caller declared them."""
    df = pd.DataFrame({"z": [1]})

    with pytest.raises(ValueError) as excinfo:
        require_columns(df, ["a", "b"])

    message = str(excinfo.value)
    assert message.index("'a'") < message.index("'b'")
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_schema.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render._schema'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/_schema.py`:

```python
"""Shared input schema validation for figure rendering.

Required-column checks were duplicated verbatim across nine figure loaders.
The column *lists* are figure-specific and live in the entrypoint scripts; only
the check itself is shared.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence

import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
def require_columns(
    df: pd.DataFrame,
    required: Sequence[str],
    context: str = "input",
) -> None:
    """Raise ValueError naming every column in `required` that `df` is missing."""
    missing = [column for column in required if column not in df.columns]

    if missing:
        rendered = ", ".join(repr(column) for column in missing)
        raise ValueError(f"Missing required columns in {context}: {rendered}")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_schema.py -v
```

Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/_schema.py workflow/tests/figure_render/test_schema.py
git commit -m "feat(figure_render): add shared required-column validation

Replaces nine verbatim copies. Column lists stay in the entrypoint scripts;
only the check is shared."
```

---

## Phase 2 — `scatter.py`

### Task 4: `render_grouped_regression_figure()` — merge correlation and orientation

This is the originally reported problem. `orientation.py`'s layout is the
baseline (spec decision), so the generic renderer must reproduce its arithmetic
exactly or the orientation figure loses pixel-identity. `correlation.py`'s
`plt.tight_layout()` is dropped — the current test run already warns
`UserWarning: This figure includes Axes that are not compatible with
tight_layout, so results might be incorrect` at `correlation.py:134`.

**Files:**
- Create: `workflow/src/figure_render/scatter.py`
- Test: `workflow/tests/figure_render/test_scatter.py`

**Interfaces:**
- Consumes: `panel_labels`, `grid_panel_size`, `PANEL_DECORATION_PX` from `figure_render._layout` (Task 2)
- Produces:

```python
REGPLOT_SCATTER_KWS: dict[str, object]  # s=3, facecolor="none", edgecolor="gray",
                                        # alpha=0.15, linewidths=0.25, rasterized=True

def render_grouped_regression_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    row_key: str = "sample",
    col_key: str = "timepoint",
    scatter_kws: Mapping[str, object] | None = None,
    panel_decoration_px: int = PANEL_DECORATION_PX,
) -> None
```

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_scatter.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic grouped regression scatter renderer.

Asserts generality (arbitrary column names, arbitrary group keys) rather than
one figure's biology. Figure-specific baselines live in the entrypoint tests.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.scatter import render_grouped_regression_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def arbitrary_frame() -> pd.DataFrame:
    """A frame whose columns share no names with any pipeline figure."""
    rng = np.random.default_rng(0)
    n = 200
    return pd.DataFrame(
        {
            "batch": ["b1"] * n + ["b2"] * n,
            "stage": (["s1"] * (n // 2) + ["s2"] * (n // 2)) * 2,
            "left_signal": np.abs(rng.normal(3, 1, 2 * n)),
            "right_signal": np.abs(rng.normal(3, 1, 2 * n)),
        }
    )


# =============================================================================
# TESTS
# =============================================================================
def test_renders_with_arbitrary_column_and_group_names(
    arbitrary_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the renderer works on columns unrelated to PBL/PBR or strand counts."""
    stem = tmp_path / "arbitrary"

    render_grouped_regression_figure(
        arbitrary_frame,
        stem,
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    assert (tmp_path / "arbitrary.pdf").exists()
    assert (tmp_path / "arbitrary.review.png").exists()
    assert (tmp_path / "arbitrary.pdf").stat().st_size > 0


def test_axis_labels_survive_regplot(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert labels are applied after drawing, since regplot resets the ones it manages."""
    import matplotlib.pyplot as plt

    render_grouped_regression_figure(
        arbitrary_frame,
        stem := tmp_path / "labels",
        x="left_signal",
        y="right_signal",
        xlabel="MY-X-LABEL",
        ylabel="MY-Y-LABEL",
        row_key="batch",
        col_key="stage",
    )

    axes = plt.gcf().get_axes()
    assert axes, "No axes on the rendered figure"
    assert any(ax.get_xlabel() == "MY-X-LABEL" for ax in axes), (
        "xlabel was lost — regplot reset it and it was not reapplied"
    )
    assert any("MY-Y-LABEL" in ax.get_ylabel() for ax in axes), (
        "ylabel was lost — regplot reset it and it was not reapplied"
    )
    assert stem.with_suffix(".pdf").exists() or (tmp_path / "labels.pdf").exists()


def test_row_key_value_appears_in_first_column_ylabel(
    arbitrary_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the row key's value is carried in the leading panel's ylabel.

    orientation.py put the sample name there rather than in the title: a full
    "{sample} {timepoint}" title is wider than the axes and pushes the measured
    panel width past the point where n_cols panels fit in a row, silently
    reflowing the grid.
    """
    import matplotlib.pyplot as plt

    render_grouped_regression_figure(
        arbitrary_frame,
        tmp_path / "rowlabel",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    ylabels = [ax.get_ylabel() for ax in plt.gcf().get_axes()]
    assert any("b1" in label for label in ylabels), f"Row key value missing from ylabels: {ylabels}"


def test_marker_size_is_inside_scatter_kws() -> None:
    """Assert the default kwargs carry `s`, since regplot drops its own `s=`.

    cnsplots regplot builds scatter_kws={"s": s, ...} internally, so supplying
    scatter_kws replaces it wholesale and loses the size.
    """
    from figure_render.scatter import REGPLOT_SCATTER_KWS

    assert "s" in REGPLOT_SCATTER_KWS, "marker size must live inside scatter_kws"
    assert REGPLOT_SCATTER_KWS["rasterized"] is True, "large point clouds must rasterize"


def test_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty input returns without raising and without writing artifacts."""
    empty = pd.DataFrame(columns=["batch", "stage", "left_signal", "right_signal"])

    render_grouped_regression_figure(
        empty,
        tmp_path / "empty",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )


def test_missing_x_column_raises(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert a bad column name fails loudly rather than drawing an empty panel."""
    with pytest.raises(ValueError, match="absent_column"):
        render_grouped_regression_figure(
            arbitrary_frame,
            tmp_path / "bad",
            x="absent_column",
            y="right_signal",
            xlabel="left",
            ylabel="right",
            row_key="batch",
            col_key="stage",
        )


def test_labels_are_keyword_only_without_defaults(
    arbitrary_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert a stale positional call fails instead of drawing wrong labels."""
    with pytest.raises(TypeError):
        render_grouped_regression_figure(arbitrary_frame, tmp_path / "stale")  # type: ignore[call-arg]


def test_panel_count_equals_group_count(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert one panel per row-key/col-key combination."""
    import matplotlib.pyplot as plt

    render_grouped_regression_figure(
        arbitrary_frame,
        tmp_path / "count",
        x="left_signal",
        y="right_signal",
        xlabel="left",
        ylabel="right",
        row_key="batch",
        col_key="stage",
    )

    expected = arbitrary_frame.groupby(["batch", "stage"]).ngroups
    assert len(plt.gcf().get_axes()) == expected
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_scatter.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render.scatter'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/scatter.py`. The layout arithmetic is copied
from `orientation.py:110-133` unchanged, so that figure stays pixel-identical.

```python
"""Generic scatter and regression figure rendering.

Supersedes ``correlation.py`` (PBL vs PBR) and ``orientation.py`` (plus vs minus
strand), which rendered the same grouped log-log regression grid for different
columns. Column names and axis labels are supplied by the caller.

Callers pass explicitly log-transformed columns when they want log-space
statistics: cnsplots ``regplot`` computes and annotates r and P from whatever
columns it is handed, so passing raw values reports a different statistic while
looking equally plausible.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Mapping
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import PANEL_DECORATION_PX, grid_panel_size, panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Verified to render ~50k points per panel legibly: density stays visible instead
# of collapsing into a solid block. Marker size must live inside these kwargs --
# regplot's own `s` argument is silently dropped when scatter_kws is supplied.
REGPLOT_SCATTER_KWS: dict[str, object] = {
    "s": 3,
    "facecolor": "none",
    "edgecolor": "gray",
    "alpha": 0.15,
    "linewidths": 0.25,
    "rasterized": True,
}


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_grouped_regression_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    row_key: str = "sample",
    col_key: str = "timepoint",
    scatter_kws: Mapping[str, object] | None = None,
    panel_decoration_px: int = PANEL_DECORATION_PX,
) -> None:
    """Render a row_key x col_key grid of regression panels for the x/y columns."""
    logger.info("Rendering grouped regression figure...")

    if df.empty:
        logger.warning("No data to plot!")
        return

    require_columns(df, [row_key, col_key, x, y], context="scatter input")

    apply_house_style()

    grouped = df.groupby([row_key, col_key], sort=True)
    n_panels = grouped.ngroups

    if n_panels == 0:
        logger.warning("No groups to plot!")
        return

    logger.info(f"Creating figure with {n_panels} panels...")

    # One row per row_key value and one column per col_key value, so a row stays
    # on one line instead of wrapping mid-group.
    n_cols = df[col_key].nunique()
    n_rows = df[row_key].nunique()
    panel_width, panel_height = grid_panel_size(
        JOURNAL_WIDTH_PX, JOURNAL_HEIGHT_PX, n_cols, n_rows, decoration_px=panel_decoration_px
    )

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    labels = panel_labels(n_panels)
    kws = dict(REGPLOT_SCATTER_KWS if scatter_kws is None else scatter_kws)

    for panel_index, ((row_value, col_value), group_df) in enumerate(grouped):
        label = labels[panel_index]
        column = panel_index % max(n_cols, 1)

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            pad_left=2,
            pad_top=2,
            margin_right=6,
            margin_bottom=18,
        )

        # The row_key value rides in the first column's ylabel rather than the
        # title: a full "{row} {col}" title is wider than the axes and pushes the
        # measured panel width past the point where n_cols panels fit in a row,
        # which silently reflows the grid.
        panel_ylabel = f"{row_value}\n{ylabel}" if column == 0 else ylabel

        def _apply_labels() -> None:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(panel_ylabel)
            ax.set_title(str(col_value))

        _apply_labels()

        if group_df.empty:
            logger.warning(f"  Panel {label}: {row_value} {col_value} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            continue

        logger.info(f"  Panel {label}: {row_value} {col_value} (n={len(group_df)})")

        cns.regplot(data=group_df, x=x, y=y, scatter_kws=kws, ax=ax)

        # regplot resets the axis labels it manages, so reapply after drawing.
        _apply_labels()

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_scatter.py -v
```

Expected: 8 passed. If `test_missing_x_column_raises` fails because
`@logger.catch` swallowed the `ValueError`, confirm the decorator carries
`reraise=True` — renderers must propagate, unlike loaders.

- [ ] **Step 5: Verify the orientation figure is still pixel-identical**

The renderer must reproduce `orientation.py` exactly. Render through the new
path and compare against the Task 1 baseline:

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
import numpy as np, pandas as pd
from pixel_baseline import ARC_DIR, compare_png
from figure_render.scatter import render_grouped_regression_figure

df = pd.read_csv(ARC_DIR / 'strand_pairs.tsv', sep='\t')
df = df[(df['plus_count'] > 0) & (df['minus_count'] > 0)].copy()
df['log10_plus_count'] = np.log10(df['plus_count'])
df['log10_minus_count'] = np.log10(df['minus_count'])

render_grouped_regression_figure(
    df, Path('tmp/check/orientation'),
    x='log10_plus_count', y='log10_minus_count',
    xlabel='log\$_{10}\$ (+) strand', ylabel='log\$_{10}\$ (-) strand',
    row_key='sample', col_key='timepoint',
)
same, diff = compare_png(Path('tmp/pixel_baseline/orientation.review.png'),
                         Path('tmp/check/orientation.review.png'))
print('IDENTICAL' if same else f'DIFFERS (max abs diff {diff})')
" 2>&1 | tail -3
```

Expected: `IDENTICAL`. If it differs, the layout arithmetic diverged from
`orientation.py:110-133` — fix before proceeding; do not accept the diff.

- [ ] **Step 6: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/scatter.py workflow/tests/figure_render/test_scatter.py
git commit -m "feat(figure_render): add generic grouped regression renderer

Merges correlation.py and orientation.py, which rendered the same grouped
log-log regression grid for different columns. Adopts orientation's explicit
grid arithmetic; correlation's plt.tight_layout() is dropped (it already
warned that the axes were incompatible with it).

Verified pixel-identical to the orientation baseline."
```

---

### Task 5: `render_scatter_grid_figure()` — the density four-panel comparison

**Files:**
- Modify: `workflow/src/figure_render/scatter.py` (append)
- Test: `workflow/tests/figure_render/test_scatter_grid.py`

**Interfaces:**
- Consumes: `panel_labels` from `figure_render._layout`; `require_columns` from `figure_render._schema`
- Produces:

```python
SCATTERPLOT_KWS: dict[str, object]  # s=6, alpha=0.4, rasterized=True — no edgecolor

@dataclass(kw_only=True, slots=True, frozen=True)
class ScatterPanel:
    x: str
    y: str
    xlabel: str
    ylabel: str
    title: str
    reference: Literal["identity", "zero", "unit_identity", "none"] = "none"
    log_scale: bool = False

def render_scatter_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    panels: Sequence[ScatterPanel],
    hue: str | None = None,
    hue_order: Sequence[str] | None = None,
    scatter_kws: Mapping[str, object] | None = None,
) -> None
```

`ScatterPanel` replaces `density.py`'s four hand-written blocks. `reference`
encodes the reference line each panel drew: `"identity"` = `y = x` scaled to the
data max (panel A), `"zero"` = `axhline(0)` (panel B), `"unit_identity"` =
`y = x` over `[0, 1]` (panel D), `"none"` = no line (panel C).

**Note on `scatter_kws`:** these defaults are deliberately *not* shared with
`REGPLOT_SCATTER_KWS`. `cns.scatterplot` always passes `edgecolor=None` to
seaborn internally, so supplying `edgecolor` here raises on the duplicate
keyword.

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_scatter_grid.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic multi-panel scatter grid renderer.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.scatter import ScatterPanel, render_scatter_grid_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def arbitrary_frame() -> pd.DataFrame:
    """A frame with generic metric names and an optional grouping column."""
    rng = np.random.default_rng(1)
    n = 300
    return pd.DataFrame(
        {
            "metric_a": rng.uniform(0, 10, n),
            "metric_b": rng.uniform(0, 10, n),
            "delta": rng.normal(0, 1, n),
            "group": rng.choice(["alpha", "beta"], n),
        }
    )


@pytest.fixture
def two_panels() -> list[ScatterPanel]:
    """Two panels exercising an identity reference line and a zero line."""
    return [
        ScatterPanel(
            x="metric_a", y="metric_b",
            xlabel="A", ylabel="B", title="A vs B",
            reference="identity",
        ),
        ScatterPanel(
            x="metric_a", y="delta",
            xlabel="A", ylabel="delta", title="Delta vs A",
            reference="zero",
        ),
    ]


# =============================================================================
# TESTS
# =============================================================================
def test_renders_one_panel_per_spec(
    arbitrary_frame: pd.DataFrame, two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert the panel count matches the number of specs supplied."""
    import matplotlib.pyplot as plt

    render_scatter_grid_figure(arbitrary_frame, tmp_path / "grid", panels=two_panels)

    assert len(plt.gcf().get_axes()) == 2
    assert (tmp_path / "grid.pdf").exists()
    assert (tmp_path / "grid.review.png").exists()


def test_scatter_kws_omit_edgecolor() -> None:
    """Assert the defaults omit edgecolor, which cns.scatterplot supplies itself.

    cns.scatterplot always passes edgecolor=None to seaborn internally, so
    supplying it here raises on the duplicate keyword.
    """
    from figure_render.scatter import SCATTERPLOT_KWS

    assert "edgecolor" not in SCATTERPLOT_KWS
    assert SCATTERPLOT_KWS["rasterized"] is True


def test_hue_is_optional(
    arbitrary_frame: pd.DataFrame, two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert rendering works both with and without a hue column."""
    render_scatter_grid_figure(arbitrary_frame, tmp_path / "nohue", panels=two_panels)
    assert (tmp_path / "nohue.pdf").exists()

    render_scatter_grid_figure(
        arbitrary_frame, tmp_path / "hue", panels=two_panels,
        hue="group", hue_order=["alpha", "beta"],
    )
    assert (tmp_path / "hue.pdf").exists()


def test_hue_order_filters_absent_levels(
    arbitrary_frame: pd.DataFrame, two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert a hue_order naming absent levels does not raise.

    density.py filtered hue_order against the observed values; the generic
    renderer must keep doing so or a project missing a viability level crashes.
    """
    render_scatter_grid_figure(
        arbitrary_frame, tmp_path / "extra_levels", panels=two_panels,
        hue="group", hue_order=["alpha", "beta", "never_present"],
    )

    assert (tmp_path / "extra_levels.pdf").exists()


def test_log_scale_panel_sets_symlog(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert a log_scale panel applies symlog to both axes, as density panel C did."""
    import matplotlib.pyplot as plt

    panels = [
        ScatterPanel(
            x="metric_a", y="metric_b",
            xlabel="A", ylabel="B", title="log",
            log_scale=True,
        )
    ]
    render_scatter_grid_figure(arbitrary_frame, tmp_path / "log", panels=panels)

    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "symlog"
    assert ax.get_yscale() == "symlog"


def test_empty_frame_renders_placeholder(
    two_panels: list[ScatterPanel], tmp_path: Path
) -> None:
    """Assert an empty frame still produces artifacts with placeholder panels.

    density.py rendered 'No valid data' panels rather than returning early, so
    the Snakemake output file always exists.
    """
    empty = pd.DataFrame(columns=["metric_a", "metric_b", "delta", "group"])

    render_scatter_grid_figure(empty, tmp_path / "empty", panels=two_panels)

    assert (tmp_path / "empty.pdf").exists(), "Snakemake requires the output to exist"


def test_missing_panel_column_raises(arbitrary_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert a panel naming an absent column fails loudly."""
    panels = [
        ScatterPanel(x="absent", y="metric_b", xlabel="A", ylabel="B", title="bad")
    ]

    with pytest.raises(ValueError, match="absent"):
        render_scatter_grid_figure(arbitrary_frame, tmp_path / "bad", panels=panels)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_scatter_grid.py -v
```

Expected: `ImportError: cannot import name 'ScatterPanel'`.

- [ ] **Step 3: Add the imports and constants**

Modify the IMPORTS and CONSTANTS sections of
`workflow/src/figure_render/scatter.py`. Replace the import block's first lines:

```python
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
```

Then append to the CONSTANTS section, after `REGPLOT_SCATTER_KWS`:

```python
# Verified for ~5k points with categorical colour-coding. Deliberately NOT shared
# with REGPLOT_SCATTER_KWS: cns.scatterplot always passes edgecolor=None to
# seaborn internally, so supplying 'edgecolor' here raises on the duplicate.
SCATTERPLOT_KWS: dict[str, object] = {
    "s": 6,
    "alpha": 0.4,
    "rasterized": True,
}


@dataclass(kw_only=True, slots=True, frozen=True)
class ScatterPanel:
    """One panel of a scatter grid: which columns, how to label, which reference line."""

    x: str
    y: str
    xlabel: str
    ylabel: str
    title: str
    reference: Literal["identity", "zero", "unit_identity", "none"] = "none"
    log_scale: bool = False
```

- [ ] **Step 4: Add the reference-line helper**

Append to the CORE LOGIC section of `workflow/src/figure_render/scatter.py`:

```python
def _draw_reference_line(ax, df: pd.DataFrame, panel: ScatterPanel) -> None:
    """Draw the panel's reference line, sized to the data where the mode requires it."""
    match panel.reference:
        case "identity":
            max_val = max(df[panel.x].max(), df[panel.y].max(), 1)
            ax.plot([0, max_val], [0, max_val], color="red", linestyle="--", alpha=0.6, linewidth=1)
        case "unit_identity":
            ax.plot([0, 1], [0, 1], color="red", linestyle="--", alpha=0.6, linewidth=1)
        case "zero":
            ax.axhline(0, color="red", linestyle="--", alpha=0.6, linewidth=1)
        case "none":
            pass
```

- [ ] **Step 5: Add the grid renderer**

Append to the CORE LOGIC section of `workflow/src/figure_render/scatter.py`:

```python
@logger.catch(reraise=True)
def render_scatter_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    panels: Sequence[ScatterPanel],
    hue: str | None = None,
    hue_order: Sequence[str] | None = None,
    scatter_kws: Mapping[str, object] | None = None,
) -> None:
    """Render one scatter panel per ScatterPanel spec, with optional hue colouring."""
    logger.info(f"Rendering scatter grid with {len(panels)} panels...")

    if not panels:
        logger.warning("No panels requested!")
        return

    # Validated even when df is empty, so a typo'd column name is caught rather
    # than silently rendering a placeholder panel.
    for panel in panels:
        require_columns(df, [panel.x, panel.y], context=f"scatter panel '{panel.title}'")

    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)

    kws = dict(SCATTERPLOT_KWS if scatter_kws is None else scatter_kws)
    labels = panel_labels(len(panels))

    for label, panel in zip(labels, panels, strict=True):
        ax = multipanel.panel(label=label)

        if df.empty:
            logger.warning(f"  Panel {label}: {panel.title} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlabel(panel.xlabel)
            ax.set_ylabel(panel.ylabel)
            ax.set_title(panel.title)
            continue

        logger.info(f"  Panel {label}: {panel.title} (n={len(df)})")

        if hue is not None and hue in df.columns:
            # Levels absent from the data must be dropped, or a project lacking
            # one of them raises inside seaborn.
            observed = set(df[hue].unique())
            order = [level for level in (hue_order or sorted(observed)) if level in observed]
            cns.scatterplot(data=df, x=panel.x, y=panel.y, hue=hue, hue_order=order, ax=ax, **kws)
            ax.legend(fontsize=6, loc="best", frameon=False)
        else:
            cns.scatterplot(data=df, x=panel.x, y=panel.y, ax=ax, **kws)

        ax.set_xlabel(panel.xlabel)
        ax.set_ylabel(panel.ylabel)
        ax.set_title(panel.title)

        _draw_reference_line(ax, df, panel)

        if panel.log_scale:
            ax.set_xscale("symlog")
            ax.set_yscale("symlog")

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 6: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_scatter_grid.py -v
```

Expected: 7 passed.

- [ ] **Step 7: Verify the density figure is still pixel-identical**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
import pandas as pd
from pixel_baseline import ARC_DIR, compare_png
from figure_render.scatter import ScatterPanel, render_scatter_grid_figure

df = pd.read_csv(ARC_DIR / 'insertion_density_analysis.tsv', sep='\t')
i, f = 'YES0', 'YES4'
panels = [
    ScatterPanel(x='insertion_density_per_kb_initial', y='insertion_density_per_kb_final',
                 xlabel=f'Insertion density per kb ({i})', ylabel=f'Insertion density per kb ({f})',
                 title='Initial vs. Final Insertion Density', reference='identity'),
    ScatterPanel(x='insertion_density_per_kb_initial', y='insertion_density_log2fc',
                 xlabel=f'Insertion density per kb ({i})', ylabel=f'log2FC density ({f} / {i})',
                 title='Density Depletion vs. Initial Coverage', reference='zero'),
    ScatterPanel(x='total_reads_initial', y='total_reads_final',
                 xlabel=f'Total reads ({i})', ylabel=f'Total reads ({f})',
                 title='Initial vs. Final Read Depth', log_scale=True),
    ScatterPanel(x='gini_coefficient_of_depth_initial', y='gini_coefficient_of_depth_final',
                 xlabel=f'Gini coefficient of depth ({i})', ylabel=f'Gini coefficient of depth ({f})',
                 title='Initial vs. Final Depth Inequality', reference='unit_identity'),
]
render_scatter_grid_figure(df, Path('tmp/check/density'), panels=panels,
                           hue='FYPOviability',
                           hue_order=['viable','inviable','condition-dependent','unknown'])
same, diff = compare_png(Path('tmp/pixel_baseline/density.review.png'),
                         Path('tmp/check/density.review.png'))
print('IDENTICAL' if same else f'DIFFERS (max abs diff {diff})')
" 2>&1 | tail -3
```

Expected: `IDENTICAL`.

**If it differs:** the likely cause is panel ordering or the reference-line
sizing. `density.py` drew the reference line *after* the scatter and labels for
panels A/B/D, and applied `symlog` after labels for panel C — the
implementation above preserves that ordering. Do not accept a diff here; density
is not on the expected-to-change list.

- [ ] **Step 8: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/scatter.py workflow/tests/figure_render/test_scatter_grid.py
git commit -m "feat(figure_render): add generic scatter grid renderer

Replaces density.py's four hand-written panel blocks with a ScatterPanel spec
sequence. Verified pixel-identical to the density baseline."
```

---

## Phase 3 — `curves.py`

### Task 6: `render_fitted_curves_figure()` — with the timepoint whitelist fix

Highest-risk task: the only figure with a *reproduced* crash. Current behaviour,
measured on this repo's projects:

| Project | LFC columns | Whitelist matches | Result today |
|---|---|---|---|
| `HD_DIT_HAP` | 5 | 5 | works |
| `Spore2YES6_1328` | 7 | 5 | `ValueError: x and y must be the same size` |
| `LD_haploid` | 8 (incl. `0h`) | 5 | same, plus `KeyError: 'DR'` |
| `Spikein` | 6 (`Spikein0..5`) | 0 | guaranteed failure |

Two pixel-affecting changes land here, both approved in the spec: the
hardcoded `#1f77b4` / `#ff7f0e` are replaced by house-style palette colours, and
`set_ylim` becomes a parameter **defaulting to `(-1.5, 8.5)`** so the default
preserves current framing.

**Files:**
- Create: `workflow/src/figure_render/curves.py`
- Test: `workflow/tests/figure_render/test_curves.py`

**Interfaces:**
- Consumes: `panel_labels` from `figure_render._layout`; `require_columns` from `figure_render._schema`
- Produces:

```python
OBSERVED_COLOR: str   # "#c84c3a" — house palette index 0
FITTED_COLOR: str     # "#2f7e8f" — house palette index 1
DEFAULT_YLIM: tuple[float, float]  # (-1.5, 8.5)

def render_fitted_curves_figure(
    observed: pd.DataFrame,
    output_stem: Path,
    *,
    x_values: Sequence[float],
    value_columns: Sequence[str],
    model: Callable[..., np.ndarray],
    model_params: Sequence[str],
    annotations: Sequence[str] = (),
    xlabel: str = "Time",
    ylabel: str = "LFC",
    ylim: tuple[float, float] | None = DEFAULT_YLIM,
    n_cols: int = 4,
) -> None
```

`observed` is one frame holding both the per-row model parameters and the
observed value columns — the caller joins the stats and LFC tables before
calling. `model` is the curve function (the pipeline passes
`depletion.curve_model.sigmoid_function`), keeping `src/figure_render/` free of
any specific curve model. `model_params` names the columns supplying that
function's arguments in order; `annotations` names extra columns to print in the
per-panel text box.

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_curves.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic fitted-curve grid renderer.

The timepoint-count tests are regression tests for a reproduced crash: the old
renderer selected value columns with a literal ['YES0'..'YES4'] whitelist and
raised "x and y must be the same size" for any project with a different number
of timepoints, or a different naming scheme.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.curves import DEFAULT_YLIM, render_fitted_curves_figure


# =============================================================================
# HELPERS
# =============================================================================
def linear_model(x: np.ndarray, slope: float, intercept: float) -> np.ndarray:
    """A trivial stand-in for the pipeline's sigmoid, so tests need no domain import."""
    return slope * np.asarray(x) + intercept


def build_frame(timepoint_names: list[str], n_rows: int = 3) -> pd.DataFrame:
    """Build an observed+params frame for the given timepoint column names."""
    rng = np.random.default_rng(2)
    data = {name: rng.normal(2, 1, n_rows) for name in timepoint_names}
    data["slope"] = rng.uniform(0.1, 0.5, n_rows)
    data["intercept"] = rng.uniform(-1, 1, n_rows)
    data["quality"] = rng.uniform(0, 1, n_rows)
    return pd.DataFrame(data)


# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize(
    "timepoint_names",
    [
        ["YES0", "YES1", "YES2", "YES3", "YES4"],
        ["YES0", "YES1", "YES2", "YES3", "YES4", "YES5", "YES6"],
        ["0h", "YES0", "YES1", "YES2", "YES3", "YES4", "YES5", "YES6"],
        ["Spikein0", "Spikein1", "Spikein2", "Spikein3", "Spikein4", "Spikein5"],
    ],
    ids=["five_yes", "seven_yes", "eight_with_0h", "spikein_naming"],
)
def test_renders_for_any_timepoint_count_and_naming(
    timepoint_names: list[str], tmp_path: Path
) -> None:
    """Assert rendering works regardless of timepoint count or naming scheme.

    Regression test: the old whitelist matched 5 of 7 columns for
    Spore2YES6_1328, 5 of 8 for LD_haploid, and 0 of 6 for Spikein.
    """
    df = build_frame(timepoint_names)
    x_values = list(range(len(timepoint_names)))

    render_fitted_curves_figure(
        df,
        tmp_path / f"curves_{len(timepoint_names)}",
        x_values=x_values,
        value_columns=timepoint_names,
        model=linear_model,
        model_params=["slope", "intercept"],
    )

    assert (tmp_path / f"curves_{len(timepoint_names)}.pdf").exists()


def test_mismatched_x_and_value_columns_raises_readable_error(tmp_path: Path) -> None:
    """Assert a length mismatch is caught at the boundary, not inside matplotlib.

    The old failure surfaced as ValueError("x and y must be the same size") from
    deep inside Axes.scatter, with no indication of which lengths disagreed.
    """
    df = build_frame(["t1", "t2", "t3"])

    with pytest.raises(ValueError) as excinfo:
        render_fitted_curves_figure(
            df,
            tmp_path / "mismatch",
            x_values=[0.0, 1.0],
            value_columns=["t1", "t2", "t3"],
            model=linear_model,
            model_params=["slope", "intercept"],
        )

    message = str(excinfo.value)
    assert "2" in message and "3" in message, f"Error must report both lengths: {message}"


def test_uses_house_palette_not_matplotlib_defaults() -> None:
    """Assert the point and curve colours come from the house palette.

    The old renderer hardcoded '#1f77b4' and '#ff7f0e', matplotlib's defaults,
    bypassing the Cell palette that apply_house_style() installs.
    """
    from figure_render.curves import FITTED_COLOR, OBSERVED_COLOR

    assert OBSERVED_COLOR.lower() != "#1f77b4"
    assert FITTED_COLOR.lower() != "#ff7f0e"
    assert OBSERVED_COLOR == "#c84c3a"
    assert FITTED_COLOR == "#2f7e8f"


def test_default_ylim_preserves_legacy_framing() -> None:
    """Assert the default y-limits stay (-1.5, 8.5) so panels remain comparable."""
    assert DEFAULT_YLIM == (-1.5, 8.5)


def test_ylim_is_overridable(tmp_path: Path) -> None:
    """Assert an explicit ylim is applied, and None means autoscale."""
    import matplotlib.pyplot as plt

    df = build_frame(["t1", "t2", "t3"])
    common = dict(
        x_values=[0.0, 1.0, 2.0],
        value_columns=["t1", "t2", "t3"],
        model=linear_model,
        model_params=["slope", "intercept"],
    )

    render_fitted_curves_figure(df, tmp_path / "fixed", ylim=(0.0, 5.0), **common)
    assert plt.gcf().get_axes()[0].get_ylim() == (0.0, 5.0)

    render_fitted_curves_figure(df, tmp_path / "auto", ylim=None, **common)
    assert plt.gcf().get_axes()[0].get_ylim() != (0.0, 5.0)


def test_panel_labels_stay_alphanumeric_past_26(tmp_path: Path) -> None:
    """Assert a 32-panel figure emits no punctuation labels.

    The pipeline samples 32 curves by default, so the old bare chr(65+i) was
    already producing '[', '\\\\' and ']' as panel labels in production.
    """
    df = build_frame(["t1", "t2", "t3"], n_rows=32)

    render_fitted_curves_figure(
        df,
        tmp_path / "many",
        x_values=[0.0, 1.0, 2.0],
        value_columns=["t1", "t2", "t3"],
        model=linear_model,
        model_params=["slope", "intercept"],
    )

    assert (tmp_path / "many.pdf").exists()


def test_annotations_are_rendered(tmp_path: Path) -> None:
    """Assert annotation columns reach the panel text box."""
    import matplotlib.pyplot as plt

    df = build_frame(["t1", "t2", "t3"], n_rows=1)

    render_fitted_curves_figure(
        df,
        tmp_path / "annot",
        x_values=[0.0, 1.0, 2.0],
        value_columns=["t1", "t2", "t3"],
        model=linear_model,
        model_params=["slope", "intercept"],
        annotations=["quality"],
    )

    texts = [t.get_text() for ax in plt.gcf().get_axes() for t in ax.texts]
    assert any("quality" in text for text in texts), f"Annotation missing from: {texts}"


def test_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty input returns without raising."""
    empty = pd.DataFrame(columns=["t1", "t2", "slope", "intercept"])

    render_fitted_curves_figure(
        empty,
        tmp_path / "empty",
        x_values=[0.0, 1.0],
        value_columns=["t1", "t2"],
        model=linear_model,
        model_params=["slope", "intercept"],
    )


def test_missing_param_column_raises(tmp_path: Path) -> None:
    """Assert an absent model-parameter column is reported by name.

    LD_haploid's stats file predates the current curve model and carries um/lam
    instead of A/DR/DL; it must fail with a readable message.
    """
    df = build_frame(["t1", "t2", "t3"])

    with pytest.raises(ValueError, match="DR"):
        render_fitted_curves_figure(
            df,
            tmp_path / "badparams",
            x_values=[0.0, 1.0, 2.0],
            value_columns=["t1", "t2", "t3"],
            model=linear_model,
            model_params=["A", "DR", "DL"],
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_curves.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render.curves'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/curves.py`:

```python
"""Generic observed-vs-fitted curve grid rendering.

Supersedes ``curve_fitting.py``. The curve model is injected by the caller, so
this module holds no specific model; the pipeline passes
``depletion.curve_model.sigmoid_function``.

Value columns come from the caller rather than a name whitelist. The old
renderer selected them with a literal ['YES0'..'YES4'] test, which silently
picked the wrong subset for any project with a different timepoint count and
raised "x and y must be the same size" from inside matplotlib.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Callable, Sequence
from pathlib import Path

import cnsplots as cns
import numpy as np
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

from ._layout import panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# House palette (Cell) indices 0 and 1. The old renderer hardcoded matplotlib's
# '#1f77b4' / '#ff7f0e', bypassing the palette apply_house_style() installs.
OBSERVED_COLOR = "#c84c3a"
FITTED_COLOR = "#2f7e8f"

# Preserved from the legacy renderer so panels stay comparable across projects.
# Callers may override, or pass None to autoscale.
DEFAULT_YLIM: tuple[float, float] = (-1.5, 8.5)

PANEL_WIDTH_PX = 180
PANEL_HEIGHT_PX = 150

# Points used to draw the fitted curve smoothly between the observed x values.
FITTED_CURVE_RESOLUTION = 100


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_fitted_curves_figure(
    observed: pd.DataFrame,
    output_stem: Path,
    *,
    x_values: Sequence[float],
    value_columns: Sequence[str],
    model: Callable[..., np.ndarray],
    model_params: Sequence[str],
    annotations: Sequence[str] = (),
    xlabel: str = "Time",
    ylabel: str = "LFC",
    ylim: tuple[float, float] | None = DEFAULT_YLIM,
    n_cols: int = 4,
) -> None:
    """Render one panel per row: observed points plus the model curve from its parameters."""
    logger.info("Rendering fitted curves figure...")

    # Checked before the empty guard so a mismatch is reported even for an empty
    # frame, and before any column access so the message names the real problem.
    if len(x_values) != len(value_columns):
        raise ValueError(
            f"x_values and value_columns must have equal length: "
            f"got {len(x_values)} x values and {len(value_columns)} value columns "
            f"({list(value_columns)})"
        )

    require_columns(
        observed,
        [*value_columns, *model_params, *annotations],
        context="fitted curves input",
    )

    if observed.empty:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    n_panels = len(observed)
    n_rows = (n_panels + n_cols - 1) // n_cols
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels} panels)...")

    multipanel = cns.multipanel(max_width=PANEL_WIDTH_PX * n_cols)
    labels = panel_labels(n_panels)

    x_array = np.asarray(x_values, dtype=float)
    x_smooth = np.linspace(x_array.min(), x_array.max(), FITTED_CURVE_RESOLUTION)

    for label, (row_index, row) in zip(labels, observed.iterrows(), strict=True):
        ax = multipanel.panel(label=label, width=PANEL_WIDTH_PX, height=PANEL_HEIGHT_PX)

        y_observed = row[list(value_columns)].to_numpy(dtype=float)
        params = [float(row[name]) for name in model_params]

        ax.scatter(
            x_array,
            y_observed,
            s=30,
            color=OBSERVED_COLOR,
            alpha=0.8,
            edgecolor="white",
            linewidths=0.5,
        )
        ax.plot(x_smooth, model(x_smooth, *params), color=FITTED_COLOR, linewidth=2, label="Fitted")

        if annotations:
            text = "\n".join(f"{name}={row[name]:.3f}" for name in annotations)
            ax.text(0.05, 0.95, text, transform=ax.transAxes, verticalalignment="top", fontsize=8)

        title = row_index if isinstance(row_index, str) else " ".join(map(str, np.atleast_1d(row_index)))
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)

        if ylim is not None:
            ax.set_ylim(*ylim)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_curves.py -v
```

Expected: 12 passed (4 parametrized + 8).

- [ ] **Step 5: Verify the reproduced crash is fixed on real data**

This is the task's real acceptance test — `Spore2YES6_1328` cannot render today.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src')
from matplotlib import use; use('Agg')
from pathlib import Path
import pandas as pd
from depletion.curve_model import sigmoid_function
from figure_render.curves import render_fitted_curves_figure

for project in ['HD_DIT_HAP', 'Spore2YES6_1328']:
    stats = pd.read_csv(
        f'projects/{project}/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv',
        sep='\t', index_col=[0,1,2,3])
    lfc = pd.read_csv(
        f'projects/{project}/results/14_insertion_level_depletion_analysis/LFC.tsv',
        sep='\t', index_col=[0,1,2,3])
    ok = stats[stats['Status'] == 'Success'].sample(n=8, random_state=42)
    time_points = [float(t) for t in ok['time_points'].iloc[0].split(',')]
    joined = ok.join(lfc[list(lfc.columns)], rsuffix='_lfc')
    render_fitted_curves_figure(
        joined, Path(f'tmp/check/curves_{project}'),
        x_values=time_points,
        value_columns=list(lfc.columns),
        model=sigmoid_function,
        model_params=['A','DR','DL'],
        annotations=['R2','RMSE'],
    )
    print(f'{project}: OK ({len(lfc.columns)} timepoints)')
" 2>&1 | grep -E "OK|Error|error" | tail -5
```

Expected: both lines print `OK`, with `HD_DIT_HAP: OK (5 timepoints)` and
`Spore2YES6_1328: OK (7 timepoints)`. The second is impossible before this task.

- [ ] **Step 6: Confirm `LD_haploid` fails with a readable message**

Its stats file predates the current curve model (`um`/`lam` instead of
`A`/`DR`/`DL`). Per the spec this is a documented limitation, not a fix — but the
error must name the missing columns instead of raising `KeyError: 'DR'`.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src')
from matplotlib import use; use('Agg')
from pathlib import Path
import pandas as pd
from depletion.curve_model import sigmoid_function
from figure_render.curves import render_fitted_curves_figure

stats = pd.read_csv('projects/LD_haploid/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv',
                    sep='\t', index_col=[0,1,2,3])
lfc = pd.read_csv('projects/LD_haploid/results/14_insertion_level_depletion_analysis/LFC.tsv',
                  sep='\t', index_col=[0,1,2,3])
ok = stats[stats['Status'] == 'Success'].sample(n=4, random_state=42)
tp = [float(t) for t in ok['time_points'].iloc[0].split(',')]
try:
    render_fitted_curves_figure(ok.join(lfc, rsuffix='_lfc'), Path('tmp/check/ld'),
        x_values=tp, value_columns=list(lfc.columns),
        model=sigmoid_function, model_params=['A','DR','DL'])
    print('UNEXPECTED: rendered without error')
except ValueError as exc:
    print('READABLE ERROR:', exc)
except KeyError as exc:
    print('BAD: still a bare KeyError:', exc)
" 2>&1 | grep -E "READABLE|BAD|UNEXPECTED"
```

Expected: `READABLE ERROR: Missing required columns in fitted curves input: 'DR', 'DL'`.

- [ ] **Step 7: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/curves.py workflow/tests/figure_render/test_curves.py
git commit -m "feat(figure_render): add generic fitted-curve renderer

Fixes a reproduced crash: the old renderer selected value columns with a
literal ['YES0'..'YES4'] whitelist, so the figure failed for four of six
projects in this repo (5 of 7 columns for Spore2YES6_1328, 5 of 8 for
LD_haploid, 0 of 6 for Spikein). Value columns now come from the caller, with
a length check that names both lengths instead of surfacing 'x and y must be
the same size' from inside matplotlib.

Also replaces hardcoded matplotlib default colours with the house palette and
makes the fixed y-limits overridable (default unchanged at (-1.5, 8.5))."
```

---

## Phase 4 — Remaining grammar modules

### Task 7: `histogram.py` — raw-values and pre-binned modes

Merges `distribution.py` (raw values, one panel per metric) and `read_counts.py`
(pre-binned counts replayed via weights, plus a cutoff line and a figure-level
footer). These share the grid-of-histograms shape but not their input form, so
the renderer carries both modes.

**Files:**
- Create: `workflow/src/figure_render/histogram.py`
- Test: `workflow/tests/figure_render/test_histogram.py`

**Interfaces:**
- Consumes: `panel_labels`, `grid_panel_size` from `figure_render._layout`; `require_columns` from `figure_render._schema`
- Produces:

```python
FOOTER_LINE_PX: int          # 9
HISTOGRAM_DECORATION_PX: int # 55 — wider than scatter panels; frequency tick labels

def render_histogram_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_columns: Sequence[str],
    bins: int,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    show_summary_stats: bool = True,
    n_cols: int = 4,
) -> None

def render_prebinned_histogram_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    row_key: str,
    col_key: str,
    left_column: str,
    right_column: str,
    count_column: str,
    xlabel: str,
    ylabel: str,
    marker_value: float | None = None,
    marker_label: str = "",
    marker_on_col_value: str | None = None,
    footer_lines: Sequence[str] = (),
    footer_header: str = "",
) -> None
```

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_histogram.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic histogram grid renderers.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.histogram import (
    render_histogram_grid_figure,
    render_prebinned_histogram_figure,
)


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def metric_frame() -> pd.DataFrame:
    """A frame of arbitrarily named numeric metrics."""
    rng = np.random.default_rng(3)
    return pd.DataFrame({f"metric_{i}": rng.normal(0, 1, 500) for i in range(6)})


@pytest.fixture
def prebinned_frame() -> pd.DataFrame:
    """A pre-binned frame with two groups and five bins each."""
    rows = []
    for sample in ("s1", "s2"):
        for stage in ("early", "late"):
            for edge in range(5):
                rows.append(
                    {
                        "sample": sample,
                        "stage": stage,
                        "bin_left": float(edge),
                        "bin_right": float(edge + 1),
                        "count": 10 * (edge + 1),
                    }
                )
    return pd.DataFrame(rows)


# =============================================================================
# TESTS — raw values mode
# =============================================================================
def test_one_panel_per_value_column(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert the panel count equals the number of value columns requested."""
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "metrics",
        value_columns=["metric_0", "metric_1", "metric_2"], bins=20,
    )

    assert len(plt.gcf().get_axes()) == 3
    assert (tmp_path / "metrics.pdf").exists()


def test_value_columns_are_caller_supplied(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert no metric-name whitelist is applied.

    distribution.py hardcoded a 15-name list ['A','DR','DL','t10',...]; a caller
    naming arbitrary columns must work.
    """
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "arb", value_columns=["metric_5"], bins=10,
    )

    assert len(plt.gcf().get_axes()) == 1


def test_summary_stats_can_be_disabled(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert the per-panel n/mean/std box is optional."""
    import matplotlib.pyplot as plt

    render_histogram_grid_figure(
        metric_frame, tmp_path / "nostats",
        value_columns=["metric_0"], bins=10, show_summary_stats=False,
    )
    assert not [t.get_text() for ax in plt.gcf().get_axes() for t in ax.texts]

    render_histogram_grid_figure(
        metric_frame, tmp_path / "stats",
        value_columns=["metric_0"], bins=10, show_summary_stats=True,
    )
    texts = [t.get_text() for ax in plt.gcf().get_axes() for t in ax.texts]
    assert any("Mean" in text for text in texts)


def test_all_nan_column_renders_placeholder(tmp_path: Path) -> None:
    """Assert a column with no finite values renders a placeholder, not a crash."""
    df = pd.DataFrame({"empty_metric": [np.nan, np.nan], "ok_metric": [1.0, 2.0]})

    render_histogram_grid_figure(
        df, tmp_path / "nan", value_columns=["empty_metric", "ok_metric"], bins=5,
    )

    assert (tmp_path / "nan.pdf").exists()


def test_missing_value_column_raises(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert an absent column is reported by name."""
    with pytest.raises(ValueError, match="absent"):
        render_histogram_grid_figure(
            metric_frame, tmp_path / "bad", value_columns=["absent"], bins=5,
        )


def test_empty_value_columns_does_not_crash(metric_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert requesting no columns returns without raising."""
    render_histogram_grid_figure(metric_frame, tmp_path / "none", value_columns=[], bins=5)


# =============================================================================
# TESTS — pre-binned mode
# =============================================================================
def test_prebinned_bar_heights_equal_stored_counts(
    prebinned_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert stored counts are replayed exactly, with no re-binning.

    This is the core guarantee of the pre-binned mode: the computation layer
    already binned the data, and the renderer must not bin again.
    """
    import matplotlib.pyplot as plt

    render_prebinned_histogram_figure(
        prebinned_frame, tmp_path / "prebinned",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
    )

    first_group = prebinned_frame[
        (prebinned_frame["sample"] == "s1") & (prebinned_frame["stage"] == "early")
    ]
    heights = np.array([p.get_height() for p in plt.gcf().get_axes()[0].patches])
    expected = first_group["count"].to_numpy().astype(float)

    assert np.array_equal(np.sort(heights), np.sort(expected)), (
        f"Bar heights {heights} differ from stored counts {expected}"
    )


def test_marker_only_on_named_column_value(
    prebinned_frame: pd.DataFrame, tmp_path: Path
) -> None:
    """Assert the cutoff marker appears only on panels matching marker_on_col_value.

    read_counts.py drew the cutoff line only on the initial-timepoint panel,
    because that is the timepoint the cutoff is applied to.
    """
    import matplotlib.pyplot as plt

    render_prebinned_histogram_figure(
        prebinned_frame, tmp_path / "marker",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
        marker_value=2.5, marker_label="Cutoff", marker_on_col_value="early",
    )

    axes = plt.gcf().get_axes()
    with_marker = sum(1 for ax in axes if ax.get_legend() is not None)

    assert with_marker == 2, f"Expected 2 'early' panels to carry the marker, got {with_marker}"


def test_footer_is_rendered(prebinned_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert figure-level footer lines reach the figure.

    Retention numbers are per sample, not per panel; read_counts.py moved them
    out of a cramped in-axes box that overlapped the histograms.
    """
    import matplotlib.pyplot as plt

    render_prebinned_histogram_figure(
        prebinned_frame, tmp_path / "footer",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
        footer_lines=["s1: 5/10 rows kept"], footer_header="Retention:",
    )

    figure_texts = [t.get_text() for t in plt.gcf().texts]
    assert any("5/10 rows kept" in text for text in figure_texts)


def test_prebinned_all_nan_group_renders_placeholder(tmp_path: Path) -> None:
    """Assert a group whose bin edges are all NaN renders a placeholder panel."""
    df = pd.DataFrame(
        {
            "sample": ["s1", "s1"],
            "stage": ["early", "late"],
            "bin_left": [0.0, np.nan],
            "bin_right": [1.0, np.nan],
            "count": [10, 0],
        }
    )

    render_prebinned_histogram_figure(
        df, tmp_path / "marker_rows",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
    )

    assert (tmp_path / "marker_rows.pdf").exists()


def test_prebinned_empty_frame_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty pre-binned frame returns without raising."""
    empty = pd.DataFrame(columns=["sample", "stage", "bin_left", "bin_right", "count"])

    render_prebinned_histogram_figure(
        empty, tmp_path / "empty",
        row_key="sample", col_key="stage",
        left_column="bin_left", right_column="bin_right", count_column="count",
        xlabel="value", ylabel="Frequency",
    )
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_histogram.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render.histogram'`.

- [ ] **Step 3: Write the module header and raw-values renderer**

Create `workflow/src/figure_render/histogram.py`:

```python
"""Generic histogram grid rendering.

Supersedes ``distribution.py`` (one panel per metric, raw values) and
``read_counts.py`` (pre-binned counts replayed via weights, with a cutoff marker
and a figure-level footer). Both are grids of histograms, but their inputs
differ in form, so each has its own entry point here.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from pathlib import Path

import cnsplots as cns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import grid_panel_size, panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Height reserved per line of the footer at the bottom of the figure.
FOOTER_LINE_PX = 9

# Frequency histograms carry wider decorations than log10 scatter panels, whose
# tick labels are single-digit.
HISTOGRAM_DECORATION_PX = 55

PANEL_WIDTH_PX = 180
PANEL_HEIGHT_PX = 150


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_histogram_grid_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    value_columns: Sequence[str],
    bins: int,
    xlabel: str = "Value",
    ylabel: str = "Frequency",
    show_summary_stats: bool = True,
    n_cols: int = 4,
) -> None:
    """Render one histogram panel per named value column."""
    logger.info("Rendering histogram grid figure...")

    require_columns(df, value_columns, context="histogram input")

    if df.empty or not value_columns:
        logger.warning("No data to plot!")
        return

    apply_house_style()

    n_panels = len(value_columns)
    n_rows = (n_panels + n_cols - 1) // n_cols
    logger.info(f"Creating figure with {n_rows}x{n_cols} grid ({n_panels} panels)...")

    multipanel = cns.multipanel(max_width=PANEL_WIDTH_PX * n_cols)
    labels = panel_labels(n_panels)

    for label, column in zip(labels, value_columns, strict=True):
        data = df[column].dropna()
        logger.info(f"  Panel {label}: {column} (n={len(data)})")

        ax = multipanel.panel(label=label, width=PANEL_WIDTH_PX, height=PANEL_HEIGHT_PX)

        if data.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(column)
            continue

        # ax.hist rather than cns.histplot: histplot requires a DataFrame, and
        # this path already holds a plain Series.
        ax.hist(data, bins=bins, alpha=0.8, edgecolor="white", linewidth=0.5)

        if show_summary_stats:
            stats_text = f"n = {len(data):,}\nMean = {data.mean():.3f}\nStd = {data.std():.3f}"
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment="top", fontsize=8)

        ax.set_title(column, fontsize=10)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 4: Add the pre-binned renderer**

Append to the CORE LOGIC section of `workflow/src/figure_render/histogram.py`:

```python
@logger.catch(reraise=True)
def render_prebinned_histogram_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    row_key: str,
    col_key: str,
    left_column: str,
    right_column: str,
    count_column: str,
    xlabel: str,
    ylabel: str,
    marker_value: float | None = None,
    marker_label: str = "",
    marker_on_col_value: str | None = None,
    footer_lines: Sequence[str] = (),
    footer_header: str = "",
) -> None:
    """Replay pre-computed bins as one histogram panel per row_key/col_key group."""
    logger.info("Rendering pre-binned histogram figure...")

    require_columns(
        df,
        [row_key, col_key, left_column, right_column, count_column],
        context="pre-binned histogram input",
    )

    if df.empty:
        logger.warning("No bins to plot!")
        return

    grouped = df.groupby([row_key, col_key], sort=True)
    if grouped.ngroups == 0:
        logger.warning("No groups to plot!")
        return

    apply_house_style()

    n_cols = df[col_key].nunique()
    n_rows = df[row_key].nunique()
    panel_width, panel_height = grid_panel_size(
        JOURNAL_WIDTH_PX, JOURNAL_HEIGHT_PX, n_cols, n_rows,
        decoration_px=HISTOGRAM_DECORATION_PX,
    )

    # multipanel resizes the figure to fit its panels, so the footer strip must be
    # reserved inside the layout via the last row's bottom margin. Shrinking the
    # requested figure height instead lands the footer on the last row's labels.
    footer_reserve_px = FOOTER_LINE_PX * (n_rows + 2) if footer_lines else 0
    last_row_index = n_rows - 1

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    labels = panel_labels(grouped.ngroups)

    for panel_index, ((row_value, col_value), group_df) in enumerate(grouped):
        label = labels[panel_index]
        is_last_row = (panel_index // max(n_cols, 1)) == last_row_index

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            pad_left=2,
            pad_top=2,
            margin_right=4,
            margin_bottom=20 + (footer_reserve_px if is_last_row else 0),
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{row_value} {col_value}")

        valid = group_df.dropna(subset=[left_column, right_column])
        total = float(valid[count_column].sum()) if not valid.empty else 0.0

        if valid.empty or total <= 0:
            logger.warning(f"  Panel {label}: {row_value} {col_value} has no valid data")
            ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
            continue

        logger.info(f"  Panel {label}: {row_value} {col_value} ({len(valid)} bins, n={int(total):,})")

        # Bin centres weighted by their counts, with bin count and range from the
        # stored edges, so nothing is re-binned. An explicit edge array cannot be
        # used: seaborn compares `bins == "auto"` before binning, which raises on
        # an array when weights are supplied.
        centers = (valid[left_column] + valid[right_column]) / 2.0
        plot_df = valid.assign(_bin_center=centers)
        binrange = (float(valid[left_column].min()), float(valid[right_column].max()))

        cns.histplot(
            data=plot_df,
            x="_bin_center",
            weights=count_column,
            bins=len(valid),
            binrange=binrange,
            ax=ax,
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{row_value} {col_value}")

        if marker_value is None:
            continue
        if marker_on_col_value is not None and col_value != marker_on_col_value:
            continue

        ax.axvline(marker_value, color="firebrick", linestyle="--", linewidth=0.8, label=marker_label)
        ax.legend(frameon=False, loc="upper right", fontsize=4)

    if footer_lines:
        footer = "\n".join([footer_header, *footer_lines]) if footer_header else "\n".join(footer_lines)
        plt.gcf().text(0.02, 0.004, footer, ha="left", va="bottom", fontsize=5, linespacing=1.4)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 5: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_histogram.py -v
```

Expected: 11 passed.

- [ ] **Step 6: Verify both figures are still pixel-identical**

Note the `_bin_center` column name: `read_counts.py`'s loader added `bin_center`
to the frame, and `x="bin_center"` was passed to `histplot`. The generic version
derives it internally under a leading-underscore name. The values are identical,
so pixels must be too.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
import numpy as np, pandas as pd
from pixel_baseline import ARC_DIR, FITTING_DIR, compare_png
from figure_render.histogram import render_histogram_grid_figure, render_prebinned_histogram_figure

# --- read_counts
dist = pd.read_csv(ARC_DIR / 'read_count_distribution.tsv', sep='\t')
stats = pd.read_csv(ARC_DIR / 'read_count_cutoff_stats.tsv', sep='\t').set_index('sample')
lines = [
    f\"{s}: {int(r['rows_kept']):,}/{int(r['original_rows']):,} rows kept ({r['pct_rows_kept']:.1f}%), \"
    f\"{int(r['counts_kept']):,}/{int(r['original_counts']):,} counts kept ({r['pct_counts_kept']:.1f}%)\"
    for s, r in stats.iterrows() if s in set(dist['sample'])
]
render_prebinned_histogram_figure(
    dist, Path('tmp/check/read_counts'),
    row_key='sample', col_key='timepoint',
    left_column='bin_left', right_column='bin_right', count_column='count',
    xlabel='log\$_{10}\$(read count)', ylabel='Frequency',
    marker_value=float(np.log10(8.0)), marker_label='Cutoff = 8', marker_on_col_value='YES0',
    footer_lines=lines, footer_header=\"Cutoff applied to 'YES0' (>= 8):\",
)
same, diff = compare_png(Path('tmp/pixel_baseline/read_counts.review.png'),
                         Path('tmp/check/read_counts.review.png'))
print('read_counts:', 'IDENTICAL' if same else f'DIFFERS ({diff})')

# --- distribution
METRICS = ['A','DR','DL','t10','t50','t90','t_window','t_inflection','y_inflection','auc','R2','RMSE','normalized_RMSE','AIC','BIC']
fit = pd.read_csv(FITTING_DIR / 'insertion_level_fitting_statistics.tsv', sep='\t', index_col=[0,1,2,3])
ok = fit[fit['Status'] == 'Success']
cols = [c for c in METRICS if c in ok.columns]
render_histogram_grid_figure(ok[cols], Path('tmp/check/distribution'), value_columns=cols, bins=30)
same, diff = compare_png(Path('tmp/pixel_baseline/distribution.review.png'),
                         Path('tmp/check/distribution.review.png'))
print('distribution:', 'IDENTICAL' if same else f'DIFFERS ({diff})')
" 2>&1 | grep -E "read_counts:|distribution:"
```

Expected: both print `IDENTICAL`.

- [ ] **Step 7: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/histogram.py workflow/tests/figure_render/test_histogram.py
git commit -m "feat(figure_render): add generic histogram grid renderers

Merges distribution.py (raw values) and read_counts.py (pre-binned weights,
cutoff marker, figure-level footer). The 15-name metric whitelist moves to the
caller. Verified pixel-identical to both baselines."
```

---

### Task 8: `series.py` — N series against a shared x

Supersedes `dispersions.py`, which plotted three dispersion series against
`normed_mean` on log-log axes.

**Files:**
- Create: `workflow/src/figure_render/series.py`
- Test: `workflow/tests/figure_render/test_series.py`

**Interfaces:**
- Consumes: `require_columns` from `figure_render._schema`
- Produces:

```python
@dataclass(kw_only=True, slots=True, frozen=True)
class Series:
    column: str
    label: str
    color: str

def render_series_scatter_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    series: Sequence[Series],
    xlabel: str,
    ylabel: str,
    title: str,
    log_x: bool = True,
    log_y: bool = True,
) -> None
```

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_series.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic multi-series scatter renderer.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.series import Series, render_series_scatter_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def series_frame() -> pd.DataFrame:
    """A frame with one shared x column and three y series."""
    rng = np.random.default_rng(4)
    n = 400
    return pd.DataFrame(
        {
            "shared_x": np.abs(rng.normal(100, 30, n)),
            "series_one": np.abs(rng.normal(0.5, 0.2, n)),
            "series_two": np.abs(rng.normal(0.4, 0.2, n)),
            "series_three": np.abs(rng.normal(0.3, 0.1, n)),
        }
    )


@pytest.fixture
def three_series() -> list[Series]:
    """Three series with explicit labels and colours."""
    return [
        Series(column="series_one", label="First", color="k"),
        Series(column="series_two", label="Second", color="b"),
        Series(column="series_three", label="Third", color="r"),
    ]


# =============================================================================
# TESTS
# =============================================================================
def test_renders_all_series_on_one_axes(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert every series lands on a single shared axes."""
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "series",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="Three series",
    )

    axes = plt.gcf().get_axes()
    assert len(axes) == 1
    assert len(axes[0].collections) == 3
    assert (tmp_path / "series.pdf").exists()


def test_legend_labels_match_series(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert legend text comes from the Series labels."""
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "legend",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )

    legend = plt.gcf().get_axes()[0].get_legend()
    assert legend is not None
    assert [t.get_text() for t in legend.get_texts()] == ["First", "Second", "Third"]


def test_legend_markers_are_scaled_and_opaque(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert legend handles are made visible despite tiny, low-alpha scatter points.

    dispersions.py used markerscale=6 and reset handle alpha to 1.0, because the
    handles otherwise inherit s=1.0 and alpha=0.12 and are invisible.
    """
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "handles",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )

    legend = plt.gcf().get_axes()[0].get_legend()
    assert all(handle.get_alpha() == 1.0 for handle in legend.legend_handles)


def test_log_scales_are_toggleable(
    series_frame: pd.DataFrame, three_series: list[Series], tmp_path: Path
) -> None:
    """Assert log_x / log_y control the axis scales."""
    import matplotlib.pyplot as plt

    render_series_scatter_figure(
        series_frame, tmp_path / "loglog",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"

    render_series_scatter_figure(
        series_frame, tmp_path / "linear",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
        log_x=False, log_y=False,
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "linear" and ax.get_yscale() == "linear"


def test_empty_frame_renders_placeholder(three_series: list[Series], tmp_path: Path) -> None:
    """Assert an empty frame still writes artifacts with a placeholder panel.

    dispersions.py rendered the placeholder and saved, so the Snakemake output
    always exists.
    """
    empty = pd.DataFrame(columns=["shared_x", "series_one", "series_two", "series_three"])

    render_series_scatter_figure(
        empty, tmp_path / "empty",
        x="shared_x", series=three_series,
        xlabel="x", ylabel="y", title="t",
    )

    assert (tmp_path / "empty.pdf").exists(), "Snakemake requires the output to exist"


def test_missing_series_column_raises(series_frame: pd.DataFrame, tmp_path: Path) -> None:
    """Assert an absent series column is reported by name rather than partially plotted."""
    bad = [Series(column="absent_series", label="X", color="k")]

    with pytest.raises(ValueError, match="absent_series"):
        render_series_scatter_figure(
            series_frame, tmp_path / "bad",
            x="shared_x", series=bad,
            xlabel="x", ylabel="y", title="t",
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_series.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render.series'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/series.py`:

```python
"""Generic multi-series scatter rendering.

Supersedes ``dispersions.py``, which plotted three named dispersion series
against a shared x column on log-log axes. Series columns, labels and colours
are supplied by the caller.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================
# Rasterize: ~94k points per series would bloat the vector PDF. Alpha is kept low
# because overlapping clouds otherwise paint over each other -- at 0.5 the second
# series hides the first entirely.
SERIES_SCATTER_KWS: dict[str, object] = {
    "s": 1.0,
    "alpha": 0.12,
    "linewidths": 0,
    "rasterized": True,
}

# Legend handles inherit the scatter's tiny size and low alpha, so they are
# scaled up and forced opaque or the legend reads as blank.
LEGEND_MARKER_SCALE = 6


@dataclass(kw_only=True, slots=True, frozen=True)
class Series:
    """One y series drawn against the shared x column."""

    column: str
    label: str
    color: str


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_series_scatter_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    series: Sequence[Series],
    xlabel: str,
    ylabel: str,
    title: str,
    log_x: bool = True,
    log_y: bool = True,
) -> None:
    """Render every series against the shared x column on one axes."""
    logger.info(f"Rendering series scatter figure ({len(series)} series)...")

    require_columns(
        df,
        [x, *(item.column for item in series)],
        context="series scatter input",
    )

    apply_house_style()

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    ax = multipanel.panel()

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if df.empty:
        logger.warning("No data to plot!")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        save_dual(output_stem)
        return

    x_values = df[x]
    for item in series:
        logger.info(f"  Series {item.label}: {item.column} (n={df[item.column].notna().sum()})")
        ax.scatter(x_values, df[item.column], c=item.color, label=item.label, **SERIES_SCATTER_KWS)

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    legend = ax.legend(loc="best", frameon=False, markerscale=LEGEND_MARKER_SCALE)
    for handle in legend.legend_handles:
        handle.set_alpha(1.0)

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_series.py -v
```

Expected: 6 passed.

- [ ] **Step 5: Verify the dispersions figure is still pixel-identical**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
import pandas as pd
from pixel_baseline import ARC_DIR, compare_png
from figure_render.series import Series, render_series_scatter_figure

df = pd.read_csv(ARC_DIR / 'dispersion_data.tsv', sep='\t', index_col=[0,1,2,3])
render_series_scatter_figure(
    df, Path('tmp/check/dispersions'),
    x='normed_mean',
    series=[
        Series(column='genewise_dispersion', label='Estimated', color='k'),
        Series(column='MAP_dispersion', label='Final', color='b'),
        Series(column='fitted_dispersion', label='Fitted', color='r'),
    ],
    xlabel='mean of normalized counts', ylabel='dispersion',
    title='DESeq2 dispersion estimates',
)
same, diff = compare_png(Path('tmp/pixel_baseline/dispersions.review.png'),
                         Path('tmp/check/dispersions.review.png'))
print('IDENTICAL' if same else f'DIFFERS (max abs diff {diff})')
" 2>&1 | tail -2
```

Expected: `IDENTICAL`.

- [ ] **Step 6: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/series.py workflow/tests/figure_render/test_series.py
git commit -m "feat(figure_render): add generic multi-series scatter renderer

Replaces dispersions.py's hardcoded three-series tuple with a caller-supplied
Series sequence. Verified pixel-identical to the dispersions baseline."
```

---

### Task 9: `ma.py` — wide and long inputs, optional significance colouring

Merges `ma_plot.py` (two wide tables, one column per timepoint, with a
vertical/horizontal orientation switch) and `ma_plot_replicates.py` (one long
table with `padj`, points coloured by significance). The two are selected by
`has_replicates` in `depletion_scoring.smk` and never run together, so neither
may regress in service of the other.

**Files:**
- Create: `workflow/src/figure_render/ma.py`
- Test: `workflow/tests/figure_render/test_ma.py`

**Interfaces:**
- Consumes: `panel_labels` from `figure_render._layout`; `require_columns` from `figure_render._schema`
- Produces:

```python
class Orientation(StrEnum):
    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"

NONSIGNIFICANT_COLOR: str  # "gray"

def render_ma_figure(
    panels: Sequence[tuple[str, pd.Series, pd.Series]],
    output_stem: Path,
    *,
    abundance_label: str,
    effect_label: str,
    title_prefix: str,
    orientation: Orientation = Orientation.VERTICAL,
    stack: bool | None = None,
    point_colors: Sequence[pd.Series | str] | None = None,
    panel_width: int = 220,
    panel_height: int = 220,
    share_axes: bool = True,
    null_effect: float = 0.0,
) -> None
```

`panels` is a sequence of `(title, abundance, effect)` triples — the caller
reshapes either input form into it, so the renderer needs no knowledge of wide
vs long tables. `point_colors`, when given, is one entry per panel: either a
per-point `Series` (replicate branch) or a single colour string.

**`orientation` and `stack` are independent, and must stay so.** Reading the
two legacy modules:

| Figure | Axis assignment | Panel arrangement |
|---|---|---|
| `ma_plot.py` VERTICAL | effect on x, abundance log-y, `axvline` | stacked (`max_width=panel_size`) |
| `ma_plot.py` HORIZONTAL | abundance log-x, effect on y, `axhline` | single row |
| `ma_plot_replicates.py` | abundance log-x, effect on y, `axhline` | **stacked** (`max_width=510`) |

So the replicate branch pairs *horizontal axes* with a *vertical stack* — a
combination a single enum cannot express. `orientation` therefore controls axis
assignment only; `stack` controls arrangement, defaulting to `True` for
`VERTICAL` and `False` for `HORIZONTAL` so `ma_plot.py` needs no extra argument,
while the replicate branch passes `stack=True` explicitly.

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_ma.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic MA plot renderer.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from figure_render.ma import Orientation, render_ma_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def two_panels() -> list[tuple[str, pd.Series, pd.Series]]:
    """Two (title, abundance, effect) triples with positive abundances."""
    rng = np.random.default_rng(5)
    return [
        (
            f"stage{i}",
            pd.Series(np.abs(rng.normal(500, 200, 300)) + 1),
            pd.Series(rng.normal(0, 1.5, 300)),
        )
        for i in (1, 2)
    ]


# =============================================================================
# TESTS
# =============================================================================
@pytest.mark.parametrize("orientation", [Orientation.VERTICAL, Orientation.HORIZONTAL])
def test_renders_in_both_orientations(
    two_panels: list[tuple[str, pd.Series, pd.Series]],
    orientation: Orientation,
    tmp_path: Path,
) -> None:
    """Assert both orientations produce artifacts."""
    stem = tmp_path / f"ma_{orientation.value}"

    render_ma_figure(
        two_panels, stem,
        abundance_label="mean of normalized counts",
        effect_label="log2 fold change",
        title_prefix="MA plot",
        orientation=orientation,
    )

    assert stem.with_name(f"{stem.name}.pdf").exists()


def test_orientation_controls_axis_assignment_and_log_axis(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert abundance goes on the log axis in whichever direction it is drawn.

    Vertical put effect on x and log-scaled y; horizontal put abundance on x and
    log-scaled x. Abundance spans orders of magnitude; the effect is symmetric.
    """
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "vert",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        orientation=Orientation.VERTICAL,
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_yscale() == "log"
    assert ax.get_xlabel() == "effect"
    assert ax.get_ylabel() == "abundance"

    render_ma_figure(
        two_panels, tmp_path / "horiz",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        orientation=Orientation.HORIZONTAL,
    )
    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log"
    assert ax.get_xlabel() == "abundance"
    assert ax.get_ylabel() == "effect"


def test_stack_is_independent_of_orientation(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert horizontal axes can be stacked vertically.

    ma_plot_replicates.py used max_width=510 (a vertical stack) while assigning
    axes horizontally: scatter(abundance, effect) with a log x-axis and axhline.
    A single enum cannot express that combination.
    """
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "horiz_stacked",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        orientation=Orientation.HORIZONTAL, stack=True,
    )

    ax = plt.gcf().get_axes()[0]
    assert ax.get_xscale() == "log", "Axis assignment must stay horizontal"
    assert ax.get_xlabel() == "abundance"
    assert (tmp_path / "horiz_stacked.pdf").exists()


def test_stack_defaults_from_orientation(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert omitting `stack` stacks for VERTICAL and rows for HORIZONTAL."""
    render_ma_figure(
        two_panels, tmp_path / "def_vert",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.VERTICAL,
    )
    assert (tmp_path / "def_vert.pdf").exists()

    render_ma_figure(
        two_panels, tmp_path / "def_horiz",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.HORIZONTAL,
    )
    assert (tmp_path / "def_horiz.pdf").exists()


def test_per_point_colors_are_accepted(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert a per-point colour Series works, as the replicate branch needs."""
    colors = [
        pd.Series(np.where(np.arange(len(abundance)) % 2 == 0, "darkred", "gray"))
        for _, abundance, _ in two_panels
    ]

    render_ma_figure(
        two_panels, tmp_path / "colored",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        point_colors=colors,
    )

    assert (tmp_path / "colored.pdf").exists()


def test_single_color_string_is_accepted(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert a scalar colour per panel works, as the no-replicate branch needs."""
    render_ma_figure(
        two_panels, tmp_path / "single",
        abundance_label="abundance", effect_label="effect", title_prefix="MA",
        point_colors=["gray", "gray"],
    )

    assert (tmp_path / "single.pdf").exists()


def test_null_effect_reference_line_orientation(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert the reference line is vertical when effect is on x, horizontal otherwise."""
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "refv",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.VERTICAL,
    )
    assert len(plt.gcf().get_axes()[0].lines) >= 1

    render_ma_figure(
        two_panels, tmp_path / "refh",
        abundance_label="a", effect_label="e", title_prefix="MA",
        orientation=Orientation.HORIZONTAL,
    )
    assert len(plt.gcf().get_axes()[0].lines) >= 1


def test_titles_use_prefix_and_panel_name(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert each panel's title combines the prefix with the panel name."""
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "titles",
        abundance_label="a", effect_label="e", title_prefix="MA plot",
    )

    titles = [ax.get_title() for ax in plt.gcf().get_axes()]
    assert "MA plot - stage1" in titles
    assert "MA plot - stage2" in titles


def test_panel_count_matches_input(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert one panel per supplied triple."""
    import matplotlib.pyplot as plt

    render_ma_figure(
        two_panels, tmp_path / "count",
        abundance_label="a", effect_label="e", title_prefix="MA",
    )

    assert len(plt.gcf().get_axes()) == 2


def test_empty_panels_does_not_crash(tmp_path: Path) -> None:
    """Assert an empty panel sequence returns without raising."""
    render_ma_figure(
        [], tmp_path / "empty",
        abundance_label="a", effect_label="e", title_prefix="MA",
    )


def test_mismatched_color_count_raises(
    two_panels: list[tuple[str, pd.Series, pd.Series]], tmp_path: Path
) -> None:
    """Assert a point_colors length mismatch is refused rather than silently zipped short."""
    with pytest.raises(ValueError, match="point_colors"):
        render_ma_figure(
            two_panels, tmp_path / "badcolors",
            abundance_label="a", effect_label="e", title_prefix="MA",
            point_colors=["gray"],
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_ma.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render.ma'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/ma.py`:

```python
"""Generic MA plot rendering.

Supersedes ``ma_plot.py`` (wide baseMean/LFC tables, one column per timepoint)
and ``ma_plot_replicates.py`` (one long table with padj-driven point colours).
Callers reshape either form into a sequence of (title, abundance, effect)
triples, so this module needs no knowledge of the input layout.

The two pipeline branches are selected by ``has_replicates`` and never run
together; both must keep working.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from collections.abc import Sequence
from enum import StrEnum
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import apply_house_style, save_dual

from ._layout import panel_labels

# =============================================================================
# CONSTANTS
# =============================================================================
class Orientation(StrEnum):
    """Panel arrangement and axis assignment for the MA plot."""

    VERTICAL = "vertical"
    HORIZONTAL = "horizontal"


NONSIGNIFICANT_COLOR = "gray"

# Rasterize: tens of thousands of points per panel would bloat the vector PDF.
MA_SCATTER_KWS: dict[str, object] = {
    "s": 1.5,
    "alpha": 0.4,
    "linewidths": 0,
    "rasterized": True,
}

# A stack needs more bottom margin than the 10 px default: each panel's title
# would otherwise collide with the x-axis label of the panel above it.
STACK_MARGIN_BOTTOM = 32
ROW_MARGIN_BOTTOM = 10

# Each panel's rendered width exceeds `panel_width`: log-scale tick labels add a
# measured ~30-40 px left reserve plus the 10 px margin_right. Layout wraps once
# the running width sum exceeds max_width, so a single row needs generous
# headroom per panel. Oversizing max_width is safe -- the figure is rendered at
# exactly max_width regardless of how much the panels fill.
ROW_WIDTH_HEADROOM_PX = 80


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def render_ma_figure(
    panels: Sequence[tuple[str, pd.Series, pd.Series]],
    output_stem: Path,
    *,
    abundance_label: str,
    effect_label: str,
    title_prefix: str,
    orientation: Orientation = Orientation.VERTICAL,
    point_colors: Sequence[pd.Series | str] | None = None,
    panel_width: int = 220,
    panel_height: int = 220,
    share_axes: bool = True,
    null_effect: float = 0.0,
) -> None:
    """Render one MA panel per (title, abundance, effect) triple."""
    logger.info(f"Rendering MA plot figure ({orientation})...")

    if not panels:
        logger.warning("No data to plot!")
        return

    if point_colors is not None and len(point_colors) != len(panels):
        raise ValueError(
            f"point_colors must have one entry per panel: "
            f"got {len(point_colors)} colours for {len(panels)} panels"
        )

    apply_house_style()

    n_panels = len(panels)
    logger.info(f"Creating figure with {n_panels} panels...")

    # Arrangement is independent of axis assignment: the replicate branch draws
    # horizontal axes in a vertical stack. Defaulting stack from orientation keeps
    # the no-replicate branch's calls unchanged.
    stacked = (orientation is Orientation.VERTICAL) if stack is None else stack

    if stacked:
        # A one-panel-wide max_width forces one panel per row. The default 10 px
        # bottom margin is too tight for a stack: each panel's title would collide
        # with the x-axis label of the panel above it.
        multipanel = cns.multipanel(max_width=panel_width)
        margin_bottom = STACK_MARGIN_BOTTOM
    else:
        multipanel = cns.multipanel(
            max_width=(panel_width + ROW_WIDTH_HEADROOM_PX) * n_panels
        )
        margin_bottom = ROW_MARGIN_BOTTOM

    labels = panel_labels(n_panels)
    colors = point_colors if point_colors is not None else [NONSIGNIFICANT_COLOR] * n_panels

    first_ax = None
    for label, (title, abundance, effect), color in zip(labels, panels, colors, strict=True):
        logger.info(f"  Panel {label}: {title} (n={len(effect)})")

        ax = multipanel.panel(
            label=label,
            width=panel_width,
            height=panel_height,
            margin_bottom=margin_bottom,
        )

        if share_axes:
            if first_ax is None:
                first_ax = ax
            else:
                ax.sharex(first_ax)
                ax.sharey(first_ax)

        match orientation:
            case Orientation.VERTICAL:
                # Effect on x, abundance (log) on y: reference line is vertical.
                ax.scatter(effect, abundance, c=color, **MA_SCATTER_KWS)
                ax.set_yscale("log")
                ax.axvline(null_effect, color="red", alpha=0.5, linestyle="--", linewidth=1, zorder=3)
                ax.set_xlabel(effect_label)
                ax.set_ylabel(abundance_label)
            case Orientation.HORIZONTAL:
                # Abundance (log) on x, effect on y: reference line is horizontal.
                ax.scatter(abundance, effect, c=color, **MA_SCATTER_KWS)
                ax.set_xscale("log")
                ax.axhline(null_effect, color="red", alpha=0.5, linestyle="--", linewidth=1, zorder=3)
                ax.set_xlabel(abundance_label)
                ax.set_ylabel(effect_label)

        ax.set_title(f"{title_prefix} - {title}")

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_ma.py -v
```

Expected: 12 passed (2 parametrized + 10).

- [ ] **Step 5: Verify both MA figures are still pixel-identical**

The replicate branch renders at `panel_width=510, panel_height=200`,
`share_axes=False`, `orientation=HORIZONTAL` **and `stack=True`**; the
no-replicate branch at `220x220` with sharing and default stacking. Both
baselines must reproduce.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
import pandas as pd
from pixel_baseline import ARC_DIR, DEPLETION_DIR, compare_png
from figure_render.ma import Orientation, render_ma_figure

# --- replicate branch (long table, padj colours)
df = pd.read_csv(ARC_DIR / 'ma_values.tsv', sep='\t')
panels, colors = [], []
for tp in sorted(df['timepoint'].unique()):
    g = df[df['timepoint'] == tp]
    panels.append((tp, g['baseMean'], g['log2FoldChange']))
    colors.append((g['padj'] < 0.05).map({True: 'darkred', False: 'gray'}))
render_ma_figure(panels, Path('tmp/check/ma_plot_replicates'),
    abundance_label='mean of normalized counts', effect_label='log2 fold change',
    title_prefix='MA plot', orientation=Orientation.HORIZONTAL, stack=True,
    point_colors=colors, panel_width=510, panel_height=200, share_axes=False)
same, diff = compare_png(Path('tmp/pixel_baseline/ma_plot_replicates.review.png'),
                         Path('tmp/check/ma_plot_replicates.review.png'))
print('ma_plot_replicates:', 'IDENTICAL' if same else f'DIFFERS ({diff})')

# --- no-replicate branch (wide tables)
bm = pd.read_csv(DEPLETION_DIR / 'baseMean.tsv', sep='\t', index_col=[0,1,2,3])
lfc = pd.read_csv(DEPLETION_DIR / 'LFC.tsv', sep='\t', index_col=[0,1,2,3])
wide = [(tp, bm[tp], lfc[tp]) for tp in lfc.columns]
render_ma_figure(wide, Path('tmp/check/ma_plot'),
    abundance_label='mean of normalized counts', effect_label='log2 fold change',
    title_prefix='MA plot', orientation=Orientation.VERTICAL)
same, diff = compare_png(Path('tmp/pixel_baseline/ma_plot.review.png'),
                         Path('tmp/check/ma_plot.review.png'))
print('ma_plot:', 'IDENTICAL' if same else f'DIFFERS ({diff})')
" 2>&1 | grep -E "ma_plot"
```

Expected: both print `IDENTICAL`.

**If `ma_plot_replicates` differs:** check three things in order — the old code
(a) stacked panels via `max_width=510` despite horizontal axes, so `stack=True`
is required; (b) did *not* share axes, so `share_axes=False` is required; and
(c) grouped with `df[df['timepoint'] == tp]` over `sorted(unique())`.

- [ ] **Step 6: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/ma.py workflow/tests/figure_render/test_ma.py
git commit -m "feat(figure_render): add generic MA plot renderer

Merges ma_plot.py (wide tables, orientation switch) and ma_plot_replicates.py
(long table, padj colours) behind a (title, abundance, effect) triple sequence.
Both pipeline branches verified pixel-identical to their baselines."
```

---

### Task 10: `composition.py` — bar overview plus per-category donuts

Supersedes `coverage.py`.

**Files:**
- Create: `workflow/src/figure_render/composition.py`
- Test: `workflow/tests/figure_render/test_composition.py`

**Interfaces:**
- Consumes: `panel_labels` from `figure_render._layout`; `require_columns` from `figure_render._schema`
- Produces:

```python
def render_composition_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    category_column: str,
    percentage_column: str,
    part_column: str,
    whole_column: str,
    total_column: str,
    part_label: str,
    whole_label: str,
    xlabel: str,
    ylabel: str,
    title: str,
    donut_unit: str = "items",
) -> None
```

`part_column` / `whole_column` are the two complementary counts (covered vs
not-covered); `part_label` / `whole_label` are their display names and fix the
donut segment order.

- [ ] **Step 1: Write the failing tests**

Create `workflow/tests/figure_render/test_composition.py`:

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tests for the generic composition (bar plus donuts) renderer.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import pandas as pd
import pytest

from figure_render.composition import render_composition_figure


# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture
def composition_frame() -> pd.DataFrame:
    """Three categories with complementary part/whole counts."""
    return pd.DataFrame(
        {
            "bucket": ["low", "mid", "high"],
            "present": [80, 50, 20],
            "absent": [20, 50, 80],
            "total": [100, 100, 100],
            "present_pct": [80.0, 50.0, 20.0],
        }
    )


@pytest.fixture
def render_kwargs() -> dict[str, str]:
    """The column mapping shared by most tests."""
    return dict(
        category_column="bucket",
        percentage_column="present_pct",
        part_column="present",
        whole_column="absent",
        total_column="total",
        part_label="Present",
        whole_label="Absent",
        xlabel="Bucket",
        ylabel="Present (%)",
        title="Presence by bucket",
    )


# =============================================================================
# TESTS
# =============================================================================
def test_panel_count_is_one_plus_categories(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert one bar overview panel plus one donut per category."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "comp", **render_kwargs)

    assert len(plt.gcf().get_axes()) == 1 + len(composition_frame)
    assert (tmp_path / "comp.pdf").exists()


def test_bar_panel_is_percentage_scaled(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert the overview panel's y-axis is pinned to 0-100."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "bar", **render_kwargs)

    assert plt.gcf().get_axes()[0].get_ylim() == (0.0, 100.0)


def test_percentage_labels_are_annotated(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert each bar carries its percentage as text."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "annot", **render_kwargs)

    texts = [t.get_text() for t in plt.gcf().get_axes()[0].texts]
    assert "80.0%" in texts


def test_zero_total_category_renders_placeholder(
    render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert a category with no items renders a placeholder instead of an empty donut."""
    df = pd.DataFrame(
        {
            "bucket": ["empty"],
            "present": [0],
            "absent": [0],
            "total": [0],
            "present_pct": [0.0],
        }
    )

    render_composition_figure(df, tmp_path / "zero", **render_kwargs)

    assert (tmp_path / "zero.pdf").exists()


def test_empty_frame_renders_placeholder(
    render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert an empty frame still writes artifacts.

    coverage.py rendered a single placeholder panel and saved, so the Snakemake
    output always exists.
    """
    empty = pd.DataFrame(columns=["bucket", "present", "absent", "total", "present_pct"])

    render_composition_figure(empty, tmp_path / "empty", **render_kwargs)

    assert (tmp_path / "empty.pdf").exists(), "Snakemake requires the output to exist"


def test_category_order_is_preserved(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert donut panels follow the frame's row order, not alphabetical order."""
    import matplotlib.pyplot as plt

    render_composition_figure(composition_frame, tmp_path / "order", **render_kwargs)

    donut_titles = [ax.get_title() for ax in plt.gcf().get_axes()[1:]]
    assert donut_titles[0].startswith("low")
    assert donut_titles[1].startswith("mid")
    assert donut_titles[2].startswith("high")


def test_missing_column_raises(
    composition_frame: pd.DataFrame, render_kwargs: dict[str, str], tmp_path: Path
) -> None:
    """Assert an absent column is reported by name."""
    kwargs = render_kwargs | {"percentage_column": "absent_pct"}

    with pytest.raises(ValueError, match="absent_pct"):
        render_composition_figure(composition_frame, tmp_path / "bad", **kwargs)
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_composition.py -v
```

Expected: collection error — `ModuleNotFoundError: No module named 'figure_render.composition'`.

- [ ] **Step 3: Write the implementation**

Create `workflow/src/figure_render/composition.py`:

```python
"""Generic composition figure rendering: a bar overview plus per-category donuts.

Supersedes ``coverage.py``, which showed gene coverage percentage per viability
category alongside covered/not-covered donuts. Column names and labels are
supplied by the caller.

Author:   Yusheng Yang (guidance) + Claude (implementation)
Date:     2026-08-25
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
from pathlib import Path

import cnsplots as cns
import pandas as pd
from loguru import logger

from figures import JOURNAL_HEIGHT_PX, JOURNAL_WIDTH_PX, apply_house_style, save_dual

from ._layout import panel_labels
from ._schema import require_columns

# =============================================================================
# CONSTANTS
# =============================================================================

# =============================================================================
# CORE LOGIC
# =============================================================================
def _status_counts_frame(part: int, whole: int, part_label: str, whole_label: str) -> pd.DataFrame:
    """Expand two complementary counts into the long-format frame cns.donutplot expects."""
    return pd.DataFrame({"status": [part_label] * part + [whole_label] * whole})


@logger.catch(reraise=True)
def render_composition_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    category_column: str,
    percentage_column: str,
    part_column: str,
    whole_column: str,
    total_column: str,
    part_label: str,
    whole_label: str,
    xlabel: str,
    ylabel: str,
    title: str,
    donut_unit: str = "items",
) -> None:
    """Render a percentage bar chart per category plus one part/whole donut each."""
    logger.info("Rendering composition figure...")

    require_columns(
        df,
        [category_column, percentage_column, part_column, whole_column, total_column],
        context="composition input",
    )

    apply_house_style()

    if df.empty:
        logger.warning("No composition data to plot!")
        cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
        multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
        ax = multipanel.panel(label="A")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax.transAxes)
        save_dual(output_stem)
        return

    n_categories = len(df)
    logger.info(f"Creating figure with {1 + n_categories} panels ({n_categories} categories)...")

    cns.figure(width=JOURNAL_WIDTH_PX, height=JOURNAL_HEIGHT_PX)
    multipanel = cns.multipanel(max_width=JOURNAL_WIDTH_PX)
    labels = panel_labels(1 + n_categories)

    # Panel A: percentage per category
    ax_bar = multipanel.panel(label=labels[0])
    cns.barplot(data=df, x=category_column, y=percentage_column, ax=ax_bar)
    ax_bar.set_ylabel(ylabel)
    ax_bar.set_xlabel(xlabel)
    ax_bar.set_title(title)
    ax_bar.set_ylim(0, 100)
    for index, row in df.reset_index(drop=True).iterrows():
        ax_bar.text(
            index, row[percentage_column], f"{row[percentage_column]:.1f}%",
            ha="center", va="bottom", fontsize=6,
        )

    # Panels B, C, ...: part vs whole per category, in the frame's row order
    for label, (_, row) in zip(labels[1:], df.iterrows(), strict=True):
        ax_donut = multipanel.panel(label=label)

        part, whole = int(row[part_column]), int(row[whole_column])
        if part + whole == 0:
            ax_donut.text(0.5, 0.5, "No valid data", ha="center", va="center", transform=ax_donut.transAxes)
            ax_donut.set_title(str(row[category_column]))
            continue

        status_df = _status_counts_frame(part, whole, part_label, whole_label)
        cns.donutplot(data=status_df, x="status", order=[part_label, whole_label], ax=ax_donut)
        ax_donut.set_title(
            f"{row[category_column]}\n({part:,}/{int(row[total_column]):,} {donut_unit})"
        )

    logger.info(f"Saving figure to {output_stem}...")
    save_dual(output_stem)
    logger.success("Figure rendering complete!")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_composition.py -v
```

Expected: 7 passed.

- [ ] **Step 5: Verify the coverage figure is still pixel-identical**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
import pandas as pd
from pixel_baseline import ARC_DIR, compare_png
from figure_render.composition import render_composition_figure

df = pd.read_csv(ARC_DIR / 'gene_coverage_stats.tsv', sep='\t')
render_composition_figure(
    df, Path('tmp/check/coverage'),
    category_column='category', percentage_column='coverage_pct',
    part_column='covered', whole_column='not_covered', total_column='total',
    part_label='Covered', whole_label='Not covered',
    xlabel='Gene viability', ylabel='Coverage (%)',
    title='Gene coverage by viability', donut_unit='genes',
)
same, diff = compare_png(Path('tmp/pixel_baseline/coverage.review.png'),
                         Path('tmp/check/coverage.review.png'))
print('IDENTICAL' if same else f'DIFFERS (max abs diff {diff})')
" 2>&1 | tail -2
```

Expected: `IDENTICAL`.

- [ ] **Step 6: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/src/figure_render/composition.py workflow/tests/figure_render/test_composition.py
git commit -m "feat(figure_render): add generic composition renderer

Replaces coverage.py's hardcoded Covered/Not-covered status names and gene
wording with caller-supplied labels. Verified pixel-identical to the coverage
baseline."
```

---

## Phase 5 — Wiring

Phases 2–4 added the grammar modules; the ten old modules are still present and
still in use. Phase 5 repoints each entrypoint, then deletes the old modules.

The ten entrypoints are split across three tasks by shape, so a reviewer can
reject one group without blocking the others. Every entrypoint keeps its
existing CLI flags exactly — no rule changes are permitted.

### Task 11: Repoint the two grouped-regression entrypoints

**Files:**
- Modify: `workflow/scripts/figures/plot_pbl_pbr_correlation.py`
- Modify: `workflow/scripts/figures/plot_insertion_orientation.py`
- Test: `workflow/tests/figure_render/test_correlation.py` (rewrite)
- Test: `workflow/tests/figure_render/test_orientation.py` (rewrite)

**Interfaces:**
- Consumes: `render_grouped_regression_figure`, `REGPLOT_SCATTER_KWS` from `figure_render.scatter` (Task 4)
- Produces: in each script, module-level constants plus `load_and_prepare_data(input_path: Path) -> pd.DataFrame`

- [ ] **Step 1: Rewrite the correlation entrypoint's data layer**

In `workflow/scripts/figures/plot_pbl_pbr_correlation.py`, replace the import of
`figure_render.correlation` with:

```python
from figure_render.scatter import render_grouped_regression_figure  # noqa: E402
```

Then add a CONSTANTS section after the imports. Note that `condition` is **not**
required: the old loader demanded it at `correlation.py:45` but never used it.

```python
# =============================================================================
# CONSTANTS
# =============================================================================
# PBL/PBR-specific knowledge lives here, not in the shared renderer.
REQUIRED_COLUMNS = ["sample", "timepoint", "pbl", "pbr"]

VALUE_COLUMNS = ["pbl", "pbr"]
X_COLUMN = "log10_pbl"
Y_COLUMN = "log10_pbr"
X_LABEL = "log$_{10}$ PBL"
Y_LABEL = "log$_{10}$ PBR"
```

- [ ] **Step 2: Add the loader to the correlation entrypoint**

Add to `workflow/scripts/figures/plot_pbl_pbr_correlation.py`, before
`parse_args()`:

```python
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load the pairs TSV, keep strictly-positive pairs, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"pairs TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} rows")

    positive = df[(df["pbl"] > 0) & (df["pbr"] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(positive)} rows")

    if positive.empty:
        logger.warning("No valid data points after filtering!")
        return positive

    # regplot annotates r and P from the columns it is handed, so the log10
    # columns must be explicit: passing raw values would report raw-space
    # correlation (r ~ -0.03) instead of the log-space PCC (r = 0.85) that this
    # figure has always shown.
    positive[X_COLUMN] = np.log10(positive["pbl"])
    positive[Y_COLUMN] = np.log10(positive["pbr"])

    return positive
```

Add these to the script's IMPORTS section:

```python
import numpy as np
import pandas as pd
```

and, after the `sys.path.append` bootstrap:

```python
from figure_render._schema import require_columns  # noqa: E402
```

- [ ] **Step 3: Repoint the correlation `main()`**

In `main()`, replace the `render_correlation_figure(df, config.output_stem)` call:

```python
        df = load_and_prepare_data(config.input_path)

        render_grouped_regression_figure(
            df,
            config.output_stem,
            x=X_COLUMN,
            y=Y_COLUMN,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            row_key="sample",
            col_key="timepoint",
        )
```

- [ ] **Step 4: Apply the same three changes to the orientation entrypoint**

In `workflow/scripts/figures/plot_insertion_orientation.py`, replace the
`figure_render.orientation` import with
`from figure_render.scatter import render_grouped_regression_figure  # noqa: E402`,
add the same `numpy`/`pandas`/`require_columns` imports, then add:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["sample", "timepoint", "plus_count", "minus_count"]

X_COLUMN = "log10_plus_count"
Y_COLUMN = "log10_minus_count"
X_LABEL = "log$_{10}$ (+) strand"
Y_LABEL = "log$_{10}$ (-) strand"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_prepare_data(input_path: Path) -> pd.DataFrame:
    """Load the strand pairs TSV, keep strictly-positive pairs, and add log10 columns."""
    logger.info(f"Loading data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"strand pairs TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} rows")

    # The computation layer already applies min > 0, but re-filtering keeps the
    # renderer safe against a hand-made input and guarantees finite log10 values.
    positive = df[(df["plus_count"] > 0) & (df["minus_count"] > 0)].copy()
    logger.info(f"After filtering to positive values: {len(positive)} rows")

    if positive.empty:
        logger.warning("No valid data points after filtering!")
        return positive

    positive[X_COLUMN] = np.log10(positive["plus_count"])
    positive[Y_COLUMN] = np.log10(positive["minus_count"])

    return positive
```

and in `main()`:

```python
        df = load_and_prepare_data(config.input_path)

        render_grouped_regression_figure(
            df,
            config.output_stem,
            x=X_COLUMN,
            y=Y_COLUMN,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            row_key="sample",
            col_key="timepoint",
        )
```

- [ ] **Step 5: Rewrite `test_orientation.py` to import from the script**

Replace the import block in `workflow/tests/figure_render/test_orientation.py`.
Every baseline assertion in that file is preserved — only the import changes,
because the loader now lives in the script.

Replace:

```python
from figure_render.orientation import load_and_prepare_data, render_orientation_figure
```

with:

```python
import importlib.util
from types import ModuleType

from figure_render.scatter import render_grouped_regression_figure


def _load_orientation_script() -> ModuleType:
    """Load the orientation CLI script by path; workflow/scripts/figures has no __init__.py."""
    path = Path(__file__).resolve().parents[2] / "scripts" / "figures" / "plot_insertion_orientation.py"
    spec = importlib.util.spec_from_file_location("_script_plot_insertion_orientation", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCRIPT = _load_orientation_script()
load_and_prepare_data = _SCRIPT.load_and_prepare_data
```

Then replace the two `render_orientation_figure(...)` call sites:

```python
    render_grouped_regression_figure(
        subset,
        output_stem,
        x="log10_plus_count",
        y="log10_minus_count",
        xlabel="log$_{10}$ (+) strand",
        ylabel="log$_{10}$ (-) strand",
        row_key="sample",
        col_key="timepoint",
    )
```

(in `test_dual_artifacts_created`, and the same call with `df` in
`test_empty_data_handling`).

- [ ] **Step 6: Rewrite `test_correlation.py` the same way**

In `workflow/tests/figure_render/test_correlation.py`, replace:

```python
from figure_render.correlation import load_and_prepare_data, render_correlation_figure
```

with the equivalent script-loading block for `plot_pbl_pbr_correlation`:

```python
import importlib.util
from types import ModuleType

from figure_render.scatter import render_grouped_regression_figure


def _load_correlation_script() -> ModuleType:
    """Load the correlation CLI script by path; workflow/scripts/figures has no __init__.py."""
    path = Path(__file__).resolve().parents[2] / "scripts" / "figures" / "plot_pbl_pbr_correlation.py"
    spec = importlib.util.spec_from_file_location("_script_plot_pbl_pbr_correlation", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SCRIPT = _load_correlation_script()
load_and_prepare_data = _SCRIPT.load_and_prepare_data
```

and replace both `render_correlation_figure(df, output_stem)` calls with:

```python
    render_grouped_regression_figure(
        df,
        output_stem,
        x="log10_pbl",
        y="log10_pbr",
        xlabel="log$_{10}$ PBL",
        ylabel="log$_{10}$ PBR",
        row_key="sample",
        col_key="timepoint",
    )
```

Also delete the `'condition'` entry from the `empty_df` column list in
`test_empty_data_handling` — it is no longer required. Keep the
`test_baseline_statistics` filter on `condition` as-is: that column still exists
in the data, it is simply not mandatory.

- [ ] **Step 7: Run both test files**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_correlation.py \
  workflow/tests/figure_render/test_orientation.py -v
```

Expected: all pass. The row-count and PCC baselines (694169 total rows,
per-group r values, 131057 rows for the PBL/PBR YES0 group with r ≈ 0.8521) must
still hold — those assert the arithmetic, which did not change.

- [ ] **Step 8: Run both figures end-to-end through the CLI**

This exercises the real Snakemake invocation path.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
ARC=projects/HD_DIT_HAP/results/18_figure_data/arc && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python \
  workflow/scripts/figures/plot_pbl_pbr_correlation.py \
  -i $ARC/pbl_pbr_pairs.tsv -o tmp/cli/correlation 2>&1 | tail -3 && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python \
  workflow/scripts/figures/plot_insertion_orientation.py \
  -i $ARC/strand_pairs.tsv -o tmp/cli/orientation 2>&1 | tail -3 && \
ls -la tmp/cli/
```

Expected: both log `Figure rendering complete!` and four files exist
(`correlation.pdf`, `correlation.review.png`, `orientation.pdf`,
`orientation.review.png`), all non-empty.

- [ ] **Step 9: Confirm orientation is pixel-identical and show correlation for review**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/tests/figure_render')
from pathlib import Path
from pixel_baseline import compare_png
same, diff = compare_png(Path('tmp/pixel_baseline/orientation.review.png'),
                         Path('tmp/cli/orientation.review.png'))
print('orientation:', 'IDENTICAL' if same else f'DIFFERS ({diff}) -- MUST NOT DIFFER')
same, diff = compare_png(Path('tmp/pixel_baseline/correlation.review.png'),
                         Path('tmp/cli/correlation.review.png'))
print('correlation:', 'identical' if same else f'differs ({diff}) -- expected, needs visual review')
"
```

Expected: `orientation: IDENTICAL` and `correlation: differs (...)`.

**Stop here and show the user** `tmp/cli/correlation.review.png` alongside
`tmp/pixel_baseline/correlation.review.png`. Correlation adopting orientation's
grid layout is the one approved visual change; it needs sign-off before Task 12.

- [ ] **Step 10: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/scripts/figures/plot_pbl_pbr_correlation.py \
        workflow/scripts/figures/plot_insertion_orientation.py \
        workflow/tests/figure_render/test_correlation.py \
        workflow/tests/figure_render/test_orientation.py
git commit -m "refactor(figures): repoint correlation and orientation to generic renderer

Column names, axis labels and required-column lists move into the entrypoints.
Drops the unused 'condition' requirement from the correlation loader.

orientation verified pixel-identical; correlation changes by design (adopts
orientation's grid layout, approved in the design doc)."
```

---

### Task 12: Repoint the curve-fitting and histogram entrypoints

**Files:**
- Modify: `workflow/scripts/figures/plot_curve_fitting.py`
- Modify: `workflow/scripts/figures/plot_distribution_of_curve_fitting.py`
- Modify: `workflow/scripts/figures/plot_read_count_distribution.py`
- Test: `workflow/tests/figure_render/test_curve_fitting.py` (rewrite)
- Test: `workflow/tests/figure_render/test_distribution.py` (rewrite)
- Test: `workflow/tests/figure_render/test_read_counts.py` (rewrite)

**Interfaces:**
- Consumes: `render_fitted_curves_figure` from `figure_render.curves` (Task 6); `render_histogram_grid_figure`, `render_prebinned_histogram_figure` from `figure_render.histogram` (Task 7)
- Produces: `load_and_sample_data(...) -> tuple[pd.DataFrame, list[float], list[str]]` in `plot_curve_fitting.py`; `load_fitting_stats(...) -> tuple[pd.DataFrame, list[str]]` in `plot_distribution_of_curve_fitting.py`; `load_distribution_data`, `load_cutoff_stats`, `format_retention_caption` in `plot_read_count_distribution.py`

- [ ] **Step 1: Rewrite `plot_curve_fitting.py`'s data layer**

The sampled stats and the LFC values are joined here, and the timepoint columns
come from the LFC table rather than a whitelist. Replace the
`figure_render.curve_fitting` import with:

```python
from depletion.curve_model import sigmoid_function  # noqa: E402
from figure_render.curves import render_fitted_curves_figure  # noqa: E402
```

Add a CONSTANTS section:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

# The sigmoid's parameters, in the order sigmoid_function() accepts them.
MODEL_PARAM_COLUMNS = ["A", "DR", "DL"]

# Printed in each panel's text box.
ANNOTATION_COLUMNS = ["R2", "RMSE"]

SUCCESS_STATUS = "Success"
X_LABEL = "Time"
Y_LABEL = "LFC"
```

- [ ] **Step 2: Add the joining loader to `plot_curve_fitting.py`**

```python
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_and_sample_data(
    fitting_stats_path: Path, lfc_path: Path, n_curves: int, random_seed: int
) -> tuple[pd.DataFrame, list[float], list[str]]:
    """Sample successful fits, join their observed LFC values, and return the x values.

    Timepoint columns are taken from the LFC table, not a name whitelist: the old
    renderer tested `col in ['YES0'..'YES4']`, which selected the wrong subset for
    any project with a different timepoint count or naming scheme.
    """
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    fitting_stats = pd.read_csv(fitting_stats_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)
    logger.info(f"Loaded {len(fitting_stats)} rows")

    successful = fitting_stats[fitting_stats["Status"] == SUCCESS_STATUS].copy()
    logger.info(f"Found {len(successful)} successful fits")

    if successful.empty:
        logger.warning("No successful fits found!")
        return successful, [], []

    sampled = successful.sample(n=min(n_curves, len(successful)), random_state=random_seed)
    logger.info(f"Sampled {len(sampled)} curves for plotting")

    logger.info(f"Loading LFC data from {lfc_path}...")
    lfc_data = pd.read_csv(lfc_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)

    timepoint_columns = list(lfc_data.columns)
    time_points = [float(t) for t in sampled["time_points"].iloc[0].split(",")]
    logger.info(f"Time points: {time_points} for columns {timepoint_columns}")

    # The stats table's own timepoint columns collide with the LFC ones, so the
    # LFC values are joined under a suffix and then renamed back.
    lfc_sampled = lfc_data.loc[sampled.index, timepoint_columns]
    joined = sampled.drop(columns=timepoint_columns, errors="ignore").join(lfc_sampled)

    return joined, time_points, timepoint_columns
```

- [ ] **Step 3: Repoint `plot_curve_fitting.py`'s `main()`**

```python
        joined, time_points, timepoint_columns = load_and_sample_data(
            config.fitting_stats_path, config.lfc_path, config.n_curves, config.random_seed
        )

        render_fitted_curves_figure(
            joined,
            config.output_stem,
            x_values=time_points,
            value_columns=timepoint_columns,
            model=sigmoid_function,
            model_params=MODEL_PARAM_COLUMNS,
            annotations=ANNOTATION_COLUMNS,
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
        )
```

- [ ] **Step 4: Rewrite `plot_distribution_of_curve_fitting.py`**

The 15-name metric whitelist moves here. Replace the
`figure_render.distribution` import with
`from figure_render.histogram import render_histogram_grid_figure  # noqa: E402`,
then add:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

SUCCESS_STATUS = "Success"

# Curve-fitting metrics worth histogramming, in display order. Columns absent
# from a given stats file are skipped.
METRIC_COLUMNS = [
    "A", "DR", "DL", "t10", "t50", "t90", "t_window", "t_inflection",
    "y_inflection", "auc", "R2", "RMSE", "normalized_RMSE", "AIC", "BIC",
]


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_fitting_stats(fitting_stats_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Load successful fits and return them with the metric columns actually present."""
    logger.info(f"Loading fitting statistics from {fitting_stats_path}...")
    df = pd.read_csv(fitting_stats_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)
    logger.info(f"Loaded {len(df)} rows")

    successful = df[df["Status"] == SUCCESS_STATUS].copy()
    logger.info(f"Found {len(successful)} successful fits")

    available = [column for column in METRIC_COLUMNS if column in successful.columns]
    logger.info(f"Found {len(available)} metric columns: {available}")

    return successful[available], available
```

and in `main()`:

```python
        df, metric_cols = load_fitting_stats(config.fitting_stats_path)

        render_histogram_grid_figure(
            df,
            config.output_stem,
            value_columns=metric_cols,
            bins=config.bins,
        )
```

- [ ] **Step 5: Rewrite `plot_read_count_distribution.py`**

Replace the `figure_render.read_counts` import block with
`from figure_render.histogram import render_prebinned_histogram_figure  # noqa: E402`,
add `from figure_render._schema import require_columns  # noqa: E402`, then add:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
DISTRIBUTION_COLUMNS = ["sample", "timepoint", "bin_left", "bin_right", "count"]
STATS_COLUMNS = [
    "sample", "original_rows", "rows_kept", "pct_rows_kept",
    "original_counts", "counts_kept", "pct_counts_kept",
]

X_LABEL = "log$_{10}$(read count)"
Y_LABEL = "Frequency"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_distribution_data(input_path: Path) -> pd.DataFrame:
    """Load the binned distribution TSV and validate its schema."""
    logger.info(f"Loading binned distribution from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, DISTRIBUTION_COLUMNS, context=f"distribution TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} bin rows")
    return df


@logger.catch
def load_cutoff_stats(stats_path: Path) -> pd.DataFrame:
    """Load the cutoff retention statistics TSV and validate its schema."""
    logger.info(f"Loading cutoff statistics from {stats_path}...")
    df = pd.read_csv(stats_path, sep="\t")

    require_columns(df, STATS_COLUMNS, context=f"cutoff stats TSV {stats_path.name}")

    logger.info(f"Loaded statistics for {len(df)} samples")
    return df


def format_retention_caption(sample: str, stats_row: pd.Series) -> str:
    """Format a one-line cutoff retention summary for a sample."""
    return (
        f"{sample}: {int(stats_row['rows_kept']):,}/{int(stats_row['original_rows']):,} rows kept "
        f"({stats_row['pct_rows_kept']:.1f}%), "
        f"{int(stats_row['counts_kept']):,}/{int(stats_row['original_counts']):,} counts kept "
        f"({stats_row['pct_counts_kept']:.1f}%)"
    )


def build_retention_footer(df: pd.DataFrame, stats_df: pd.DataFrame) -> list[str]:
    """Build one retention line per sample present in the distribution data."""
    by_sample = stats_df.set_index("sample")
    present = set(df["sample"].unique())

    return [
        format_retention_caption(sample, row)
        for sample, row in by_sample.iterrows()
        if sample in present
    ]
```

and in `main()`:

```python
        df = load_distribution_data(config.input_path)
        stats_df = load_cutoff_stats(config.stats_path)

        render_prebinned_histogram_figure(
            df,
            config.output_stem,
            row_key="sample",
            col_key="timepoint",
            left_column="bin_left",
            right_column="bin_right",
            count_column="count",
            xlabel=X_LABEL,
            ylabel=Y_LABEL,
            marker_value=float(np.log10(config.cutoff)),
            marker_label=f"Cutoff = {config.cutoff:.2g}",
            marker_on_col_value=config.initial_time_point,
            footer_lines=build_retention_footer(df, stats_df),
            footer_header=f"Cutoff applied to '{config.initial_time_point}' (>= {config.cutoff:.2g}):",
        )
```

Add `import numpy as np` and `import pandas as pd` to the script's IMPORTS section.

- [ ] **Step 6: Rewrite the three test files' import and call sites**

Each keeps every baseline assertion; only imports and render calls change.

In `workflow/tests/figure_render/test_curve_fitting.py`, replace the
`figure_render.curve_fitting` import with a script loader for
`plot_curve_fitting` (same `importlib.util` pattern as Task 11, Step 5), and
update the three-value unpacking — the loader now returns
`(joined, time_points, timepoint_columns)` rather than
`(sampled_stats, lfc_sampled, time_points)`:

```python
    joined, time_points, timepoint_columns = load_and_sample_data(
        real_stats_path, real_lfc_path, 32, 42
    )

    assert len(joined) == 32, f"Expected 32 sampled curves, got {len(joined)}"
    assert len(time_points) == 5, f"Expected 5 time points, got {len(time_points)}"
    assert time_points == [0.0, 2.352, 5.588, 9.104, 12.48], f"Unexpected time points: {time_points}"
    assert len(timepoint_columns) == len(time_points), "One value column per time point"
    assert (joined["Status"] == "Success").all(), "Not all sampled curves are successful"
```

and the render call:

```python
    render_fitted_curves_figure(
        joined,
        output_stem,
        x_values=time_points,
        value_columns=timepoint_columns,
        model=sigmoid_function,
        model_params=["A", "DR", "DL"],
        annotations=["R2", "RMSE"],
    )
```

with these test-file imports:

```python
from depletion.curve_model import sigmoid_function
from figure_render.curves import render_fitted_curves_figure
```

In `test_distribution.py`, load `plot_distribution_of_curve_fitting` for
`load_fitting_stats` and replace the render call with
`render_histogram_grid_figure(df, output_stem, value_columns=metric_cols, bins=30)`.
Keep the 93560-row and 15-column baselines and the mean R2/RMSE assertions.

In `test_read_counts.py`, load `plot_read_count_distribution` for
`load_distribution_data` / `load_cutoff_stats`, and replace both
`render_distribution_figure(df, stats_df, stem, "YES0", 8.0)` calls with the
`render_prebinned_histogram_figure(...)` form from Step 5. Two adjustments:

- `test_renderer_does_not_rebin` computes `bin_center` itself; it must now derive
  the column locally rather than expecting the loader to add it:
  ```python
  group = group.assign(bin_center=(group["bin_left"] + group["bin_right"]) / 2.0)
  ```
- `test_bin_count_is_fifty` and the bin-edge baselines are unchanged.

- [ ] **Step 7: Run the three test files**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_curve_fitting.py \
  workflow/tests/figure_render/test_distribution.py \
  workflow/tests/figure_render/test_read_counts.py -v
```

Expected: all pass.

- [ ] **Step 8: Run all three figures through the CLI**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
ARC=projects/HD_DIT_HAP/results/18_figure_data/arc && \
FIT=projects/HD_DIT_HAP/results/15_insertion_level_curve_fitting && \
DEP=projects/HD_DIT_HAP/results/14_insertion_level_depletion_analysis && \
PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python && \
$PY workflow/scripts/figures/plot_curve_fitting.py \
  -s $FIT/insertion_level_fitting_statistics.tsv -l $DEP/LFC.tsv \
  -o tmp/cli/curve_fitting 2>&1 | tail -2 && \
$PY workflow/scripts/figures/plot_distribution_of_curve_fitting.py \
  -i $FIT/insertion_level_fitting_statistics.tsv -o tmp/cli/distribution 2>&1 | tail -2 && \
$PY workflow/scripts/figures/plot_read_count_distribution.py \
  -i $ARC/read_count_distribution.tsv -s $ARC/read_count_cutoff_stats.tsv \
  -t YES0 -c 8 -o tmp/cli/read_counts 2>&1 | tail -2
```

Expected: three `Figure rendering complete!` lines.

- [ ] **Step 9: Verify pixels — two identical, one expected to change**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/tests/figure_render')
from pathlib import Path
from pixel_baseline import compare_png
for name, must_match in [('distribution', True), ('read_counts', True), ('curve_fitting', False)]:
    same, diff = compare_png(Path(f'tmp/pixel_baseline/{name}.review.png'),
                             Path(f'tmp/cli/{name}.review.png'))
    verdict = 'IDENTICAL' if same else f'differs ({diff})'
    flag = '' if same or not must_match else '  <-- MUST NOT DIFFER'
    print(f'{name}: {verdict}{flag}')
"
```

Expected: `distribution: IDENTICAL`, `read_counts: IDENTICAL`,
`curve_fitting: differs (...)` — the approved palette change.

**Stop and show the user** `tmp/cli/curve_fitting.review.png` next to the
baseline. Two things to confirm: the house-palette colours, and that panel labels
past 26 are now `A1`, `A2`... instead of `[`, `\`, `]`.

- [ ] **Step 10: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/scripts/figures/plot_curve_fitting.py \
        workflow/scripts/figures/plot_distribution_of_curve_fitting.py \
        workflow/scripts/figures/plot_read_count_distribution.py \
        workflow/tests/figure_render/test_curve_fitting.py \
        workflow/tests/figure_render/test_distribution.py \
        workflow/tests/figure_render/test_read_counts.py
git commit -m "refactor(figures): repoint curve fitting and histogram entrypoints

Timepoint columns now come from the LFC table instead of a ['YES0'..'YES4']
whitelist, and the 15-name metric list moves into its entrypoint.

distribution and read_counts verified pixel-identical; curve_fitting changes by
design (house palette, and panel labels past 26 are alphanumeric)."
```

---

### Task 13: Repoint the remaining four entrypoints

**Files:**
- Modify: `workflow/scripts/figures/plot_ma_plot.py`
- Modify: `workflow/scripts/figures/plot_ma_plot_replicates.py`
- Modify: `workflow/scripts/figures/plot_dispersions.py`
- Modify: `workflow/scripts/figures/plot_gene_coverage.py`
- Modify: `workflow/scripts/figures/plot_insertion_density.py`
- Test: rewrite `test_ma_plot.py`, `test_ma_plot_replicates.py`, `test_dispersions.py`, `test_coverage.py`, `test_density.py`

**Interfaces:**
- Consumes: `render_ma_figure`, `Orientation` from `figure_render.ma` (Task 9); `render_series_scatter_figure`, `Series` from `figure_render.series` (Task 8); `render_composition_figure` from `figure_render.composition` (Task 10); `render_scatter_grid_figure`, `ScatterPanel` from `figure_render.scatter` (Task 5)
- Produces: a loader plus a CONSTANTS block in each of the five scripts

- [ ] **Step 1: Rewrite `plot_ma_plot.py`**

Replace the `figure_render.ma_plot` import with
`from figure_render.ma import Orientation, render_ma_figure  # noqa: E402`, then add:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

ABUNDANCE_LABEL = "mean of normalized counts"
EFFECT_LABEL = "log2 fold change"
TITLE_PREFIX = "MA plot"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_ma_data(basemean_path: Path, lfc_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the wide baseMean/LFC tables and validate that their columns match."""
    logger.info(f"Loading baseMean from {basemean_path}...")
    basemean_df = pd.read_csv(basemean_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)

    logger.info(f"Loading LFC from {lfc_path}...")
    lfc_df = pd.read_csv(lfc_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)

    missing_cols = [col for col in lfc_df.columns if col not in basemean_df.columns]
    if missing_cols:
        raise ValueError(f"Timepoints present in LFC but missing from baseMean: {missing_cols}")

    logger.info(f"Loaded {len(lfc_df)} insertions across {len(lfc_df.columns)} timepoints")
    return basemean_df, lfc_df


def build_ma_panels(
    basemean_df: pd.DataFrame, lfc_df: pd.DataFrame
) -> list[tuple[str, pd.Series, pd.Series]]:
    """Reshape the wide tables into (title, abundance, effect) triples, one per timepoint."""
    return [(timepoint, basemean_df[timepoint], lfc_df[timepoint]) for timepoint in lfc_df.columns]
```

and in `main()`, keeping the two-orientation behaviour:

```python
        basemean_df, lfc_df = load_ma_data(config.basemean_path, config.lfc_path)
        panels = build_ma_panels(basemean_df, lfc_df)

        render_ma_figure(
            panels, config.output_stem,
            abundance_label=ABUNDANCE_LABEL, effect_label=EFFECT_LABEL,
            title_prefix=TITLE_PREFIX, orientation=Orientation.VERTICAL,
        )

        horizontal_stem = config.output_stem.with_name(f"{config.output_stem.name}_horizontal")
        render_ma_figure(
            panels, horizontal_stem,
            abundance_label=ABUNDANCE_LABEL, effect_label=EFFECT_LABEL,
            title_prefix=TITLE_PREFIX, orientation=Orientation.HORIZONTAL,
        )
```

- [ ] **Step 2: Rewrite `plot_ma_plot_replicates.py`**

Replace the `figure_render.ma_plot_replicates` import with
`from figure_render.ma import Orientation, render_ma_figure  # noqa: E402` and add
`from figure_render._schema import require_columns  # noqa: E402`:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["timepoint", "baseMean", "log2FoldChange", "padj"]

# pydeseq2's DeseqStats default alpha; the project config sets no override.
PADJ_THRESHOLD = 0.05
SIGNIFICANT_COLOR = "darkred"
NONSIGNIFICANT_COLOR = "gray"

ABUNDANCE_LABEL = "mean of normalized counts"
EFFECT_LABEL = "log2 fold change"
TITLE_PREFIX = "MA plot"

# The legacy replicate figure used a full-width vertical stack of short panels and
# did not share axes between them, while assigning axes horizontally (abundance on
# a log x-axis). Hence orientation=HORIZONTAL together with stack=True.
PANEL_WIDTH = 510
PANEL_HEIGHT = 200


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_ma_data(ma_values_path: Path) -> pd.DataFrame:
    """Load the long-format MA values TSV and validate its schema."""
    logger.info(f"Loading MA values from {ma_values_path}...")
    df = pd.read_csv(ma_values_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"MA values TSV {ma_values_path.name}")

    logger.info(f"Loaded {len(df)} rows")
    return df


def significance_colors(padj: pd.Series) -> pd.Series:
    """Map adjusted p-values to point colours; NaN padj stays non-significant like pydeseq2."""
    # (padj < threshold) is False for NaN, so filtered-out insertions render gray.
    return (padj < PADJ_THRESHOLD).map({True: SIGNIFICANT_COLOR, False: NONSIGNIFICANT_COLOR})


def build_ma_panels(
    df: pd.DataFrame,
) -> tuple[list[tuple[str, pd.Series, pd.Series]], list[pd.Series]]:
    """Split the long table into per-timepoint triples plus their point colours."""
    panels: list[tuple[str, pd.Series, pd.Series]] = []
    colors: list[pd.Series] = []

    for timepoint in sorted(df["timepoint"].unique()):
        group = df[df["timepoint"] == timepoint]
        panels.append((timepoint, group["baseMean"], group["log2FoldChange"]))
        colors.append(significance_colors(group["padj"]))

    return panels, colors
```

and in `main()`:

```python
        df = load_ma_data(config.ma_values_path)
        panels, colors = build_ma_panels(df)

        render_ma_figure(
            panels, config.output_stem,
            abundance_label=ABUNDANCE_LABEL, effect_label=EFFECT_LABEL,
            title_prefix=TITLE_PREFIX,
            orientation=Orientation.HORIZONTAL, stack=True,
            point_colors=colors,
            panel_width=PANEL_WIDTH, panel_height=PANEL_HEIGHT, share_axes=False,
        )
```

- [ ] **Step 3: Rewrite `plot_dispersions.py`**

Replace the `figure_render.dispersions` import with
`from figure_render.series import Series, render_series_scatter_figure  # noqa: E402`
and add `from figure_render._schema import require_columns  # noqa: E402`:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
INSERTION_INDEX_COLUMNS = [0, 1, 2, 3]

X_COLUMN = "normed_mean"

# Series order, labels and colours mirror pydeseq2's plot_dispersions(): it passes
# [genewise, MAP, fitted] with labels ["Estimated", "Final", "Fitted"] and the
# matplotlib colour string "kbr" mapped positionally.
DISPERSION_SERIES = [
    Series(column="genewise_dispersion", label="Estimated", color="k"),
    Series(column="MAP_dispersion", label="Final", color="b"),
    Series(column="fitted_dispersion", label="Fitted", color="r"),
]

REQUIRED_COLUMNS = [X_COLUMN, *(item.column for item in DISPERSION_SERIES)]

X_LABEL = "mean of normalized counts"
Y_LABEL = "dispersion"
TITLE = "DESeq2 dispersion estimates"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_dispersion_data(dispersion_data_path: Path) -> pd.DataFrame:
    """Load the dispersion figure-data TSV and validate its schema."""
    logger.info(f"Loading dispersion data from {dispersion_data_path}...")
    df = pd.read_csv(dispersion_data_path, sep="\t", index_col=INSERTION_INDEX_COLUMNS)

    require_columns(df, REQUIRED_COLUMNS, context=f"dispersion TSV {dispersion_data_path.name}")

    logger.info(f"Loaded {len(df)} rows")
    return df
```

and in `main()`:

```python
        df = load_dispersion_data(config.dispersion_data_path)

        render_series_scatter_figure(
            df, config.output_stem,
            x=X_COLUMN, series=DISPERSION_SERIES,
            xlabel=X_LABEL, ylabel=Y_LABEL, title=TITLE,
        )
```

- [ ] **Step 4: Rewrite `plot_gene_coverage.py`**

Replace the `figure_render.coverage` import with
`from figure_render.composition import render_composition_figure  # noqa: E402`
and add `from figure_render._schema import require_columns  # noqa: E402`:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = ["category", "covered", "not_covered", "total", "coverage_pct"]

COVERED_LABEL = "Covered"
NOT_COVERED_LABEL = "Not covered"

X_LABEL = "Gene viability"
Y_LABEL = "Coverage (%)"
TITLE = "Gene coverage by viability"
DONUT_UNIT = "genes"


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_coverage_data(input_path: Path) -> pd.DataFrame:
    """Load the coverage statistics TSV and validate its schema."""
    logger.info(f"Loading coverage statistics from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"coverage TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} viability categories")
    return df
```

and in `main()`:

```python
        df = load_coverage_data(config.input_path)

        render_composition_figure(
            df, config.output_stem,
            category_column="category", percentage_column="coverage_pct",
            part_column="covered", whole_column="not_covered", total_column="total",
            part_label=COVERED_LABEL, whole_label=NOT_COVERED_LABEL,
            xlabel=X_LABEL, ylabel=Y_LABEL, title=TITLE, donut_unit=DONUT_UNIT,
        )
```

- [ ] **Step 5: Rewrite `plot_insertion_density.py`**

The four panel definitions move here. Replace the `figure_render.density` import
with
`from figure_render.scatter import ScatterPanel, render_scatter_grid_figure  # noqa: E402`
and add `from figure_render._schema import require_columns  # noqa: E402`:

```python
# =============================================================================
# CONSTANTS
# =============================================================================
REQUIRED_COLUMNS = [
    "insertion_density_per_kb_initial",
    "insertion_density_per_kb_final",
    "insertion_density_log2fc",
    "total_reads_initial",
    "total_reads_final",
    "gini_coefficient_of_depth_initial",
    "gini_coefficient_of_depth_final",
]

VIABILITY_COLUMN = "FYPOviability"
VIABILITY_HUE_ORDER = ["viable", "inviable", "condition-dependent", "unknown"]


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_density_data(input_path: Path) -> pd.DataFrame:
    """Load the density statistics TSV and validate its schema."""
    logger.info(f"Loading density statistics from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")

    require_columns(df, REQUIRED_COLUMNS, context=f"density TSV {input_path.name}")

    logger.info(f"Loaded {len(df)} genes")
    return df


def build_density_panels(initial: str, final: str) -> list[ScatterPanel]:
    """Build the four density comparison panels for the given timepoint names."""
    return [
        ScatterPanel(
            x="insertion_density_per_kb_initial", y="insertion_density_per_kb_final",
            xlabel=f"Insertion density per kb ({initial})",
            ylabel=f"Insertion density per kb ({final})",
            title="Initial vs. Final Insertion Density",
            reference="identity",
        ),
        ScatterPanel(
            x="insertion_density_per_kb_initial", y="insertion_density_log2fc",
            xlabel=f"Insertion density per kb ({initial})",
            ylabel=f"log2FC density ({final} / {initial})",
            title="Density Depletion vs. Initial Coverage",
            reference="zero",
        ),
        ScatterPanel(
            x="total_reads_initial", y="total_reads_final",
            xlabel=f"Total reads ({initial})", ylabel=f"Total reads ({final})",
            title="Initial vs. Final Read Depth",
            log_scale=True,
        ),
        ScatterPanel(
            x="gini_coefficient_of_depth_initial", y="gini_coefficient_of_depth_final",
            xlabel=f"Gini coefficient of depth ({initial})",
            ylabel=f"Gini coefficient of depth ({final})",
            title="Initial vs. Final Depth Inequality",
            reference="unit_identity",
        ),
    ]
```

and in `main()`:

```python
        df = load_density_data(config.input_path)

        render_scatter_grid_figure(
            df, config.output_stem,
            panels=build_density_panels(config.initial_timepoint, config.final_timepoint),
            hue=VIABILITY_COLUMN,
            hue_order=VIABILITY_HUE_ORDER,
        )
```

- [ ] **Step 6: Rewrite the five test files' imports and call sites**

Each keeps every baseline assertion. Apply the `importlib.util` script-loader
pattern from Task 11 Step 5, then update render calls:

- `test_ma_plot.py`: load `plot_ma_plot` for `load_ma_data` and
  `build_ma_panels`; import `Orientation, render_ma_figure` from
  `figure_render.ma`. Replace `render_ma_figure(basemean, lfc, stem, orientation)`
  with the panels form:
  ```python
  render_ma_figure(
      _SCRIPT.build_ma_panels(basemean_df, lfc_df), output_stem,
      abundance_label="mean of normalized counts", effect_label="log2 fold change",
      title_prefix="MA plot", orientation=orientation,
  )
  ```
  `test_load_ma_data_rejects_mismatched_timepoints` keeps asserting the loader
  raises on a column mismatch — that check moved into the script unchanged.
- `test_ma_plot_replicates.py`: load `plot_ma_plot_replicates` for
  `load_ma_data`, `significance_colors`, `build_ma_panels`. Replace the render
  call with the panels+colors form from Step 2 (including
  `orientation=Orientation.HORIZONTAL`, `stack=True`, `share_axes=False`,
  `panel_width=510`, `panel_height=200`). Keep
  `test_nan_padj_renders_nonsignificant` as-is: `significance_colors` is
  unchanged, only relocated.
- `test_dispersions.py`: load `plot_dispersions` for `load_dispersion_data` and
  `X_COLUMN`; call `render_series_scatter_figure` with `_SCRIPT.DISPERSION_SERIES`.
- `test_coverage.py`: load `plot_gene_coverage` for `load_coverage_data`; call
  `render_composition_figure` with the Step 4 arguments.
- `test_density.py`: load `plot_insertion_density` for `load_density_data` and
  `build_density_panels`; call `render_scatter_grid_figure` with
  `panels=_SCRIPT.build_density_panels("YES0", "YES4")`.

The three `test_missing_*_returns_none` tests still pass unchanged: the loaders
kept their `@logger.catch` decorator, so a `ValueError` is still logged and
converted to a `None` return.

- [ ] **Step 7: Run the five test files**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests/figure_render/test_ma_plot.py \
  workflow/tests/figure_render/test_ma_plot_replicates.py \
  workflow/tests/figure_render/test_dispersions.py \
  workflow/tests/figure_render/test_coverage.py \
  workflow/tests/figure_render/test_density.py -v
```

Expected: all pass.

- [ ] **Step 8: Run the four figures through the CLI and verify pixels**

All four must be byte-identical — none is on the expected-to-change list.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
ARC=projects/HD_DIT_HAP/results/18_figure_data/arc && \
DEP=projects/HD_DIT_HAP/results/14_insertion_level_depletion_analysis && \
PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python && \
$PY workflow/scripts/figures/plot_dispersions.py \
  -i $ARC/dispersion_data.tsv -o tmp/cli/dispersions 2>&1 | tail -1 && \
$PY workflow/scripts/figures/plot_gene_coverage.py \
  -i $ARC/gene_coverage_stats.tsv -o tmp/cli/coverage 2>&1 | tail -1 && \
$PY workflow/scripts/figures/plot_insertion_density.py \
  -i $ARC/insertion_density_analysis.tsv -t YES0 -f YES4 \
  -o tmp/cli/density 2>&1 | tail -1 && \
$PY workflow/scripts/figures/plot_ma_plot_replicates.py \
  -i $ARC/ma_values.tsv -o tmp/cli/ma_plot_replicates 2>&1 | tail -1 && \
$PY workflow/scripts/figures/plot_ma_plot.py \
  -b $DEP/baseMean.tsv -l $DEP/LFC.tsv -o tmp/cli/ma_plot 2>&1 | tail -1 && \
$PY -c "
import sys; sys.path.insert(0,'workflow/tests/figure_render')
from pathlib import Path
from pixel_baseline import compare_png
for name in ['dispersions','coverage','density','ma_plot_replicates','ma_plot']:
    same, diff = compare_png(Path(f'tmp/pixel_baseline/{name}.review.png'),
                             Path(f'tmp/cli/{name}.review.png'))
    print(f'{name}:', 'IDENTICAL' if same else f'DIFFERS ({diff})  <-- MUST NOT DIFFER')
"
```

Expected: all five print `IDENTICAL`.

- [ ] **Step 9: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add workflow/scripts/figures/plot_ma_plot.py \
        workflow/scripts/figures/plot_ma_plot_replicates.py \
        workflow/scripts/figures/plot_dispersions.py \
        workflow/scripts/figures/plot_gene_coverage.py \
        workflow/scripts/figures/plot_insertion_density.py \
        workflow/tests/figure_render/test_ma_plot.py \
        workflow/tests/figure_render/test_ma_plot_replicates.py \
        workflow/tests/figure_render/test_dispersions.py \
        workflow/tests/figure_render/test_coverage.py \
        workflow/tests/figure_render/test_density.py
git commit -m "refactor(figures): repoint remaining five entrypoints

MA plots (both branches), dispersions, coverage and density now call the
generic renderers, with their column names, labels and panel definitions held
in the entrypoints. All five verified pixel-identical."
```

---

### Task 14: Delete the ten superseded modules

Nothing imports them after Task 13. This task proves that and removes them.

**Files:**
- Delete: all ten old modules in `workflow/src/figure_render/`

**Interfaces:**
- Consumes: everything from Tasks 4–13
- Produces: nothing

- [ ] **Step 1: Prove nothing references the old modules**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
grep -rn "figure_render\.\(correlation\|orientation\|density\|distribution\|read_counts\|ma_plot\|ma_plot_replicates\|dispersions\|curve_fitting\|coverage\)" \
  workflow/ --include="*.py" | grep -v "workflow/src/figure_render/" | grep -v pixel_baseline.py
```

Expected: no output. Any hit must be repointed before deleting.

- [ ] **Step 2: Update the baseline harness to the new API**

`pixel_baseline.py`'s shims still import the old modules. Repoint each to the
entrypoint script's loader plus the generic renderer, mirroring the CLI commands
verified in Tasks 11–13. Replace the whole shim block with:

```python
def _load_script(stem: str):
    """Load a figure CLI script by stem; workflow/scripts/figures has no __init__.py."""
    import importlib.util

    path = PROJECT_ROOT / "workflow" / "scripts" / "figures" / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"_script_{stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _render_correlation(stem: Path) -> None:
    from figure_render.scatter import render_grouped_regression_figure
    script = _load_script("plot_pbl_pbr_correlation")
    df = script.load_and_prepare_data(ARC_DIR / "pbl_pbr_pairs.tsv")
    render_grouped_regression_figure(
        df, stem, x=script.X_COLUMN, y=script.Y_COLUMN,
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL,
        row_key="sample", col_key="timepoint",
    )


def _render_orientation(stem: Path) -> None:
    from figure_render.scatter import render_grouped_regression_figure
    script = _load_script("plot_insertion_orientation")
    df = script.load_and_prepare_data(ARC_DIR / "strand_pairs.tsv")
    render_grouped_regression_figure(
        df, stem, x=script.X_COLUMN, y=script.Y_COLUMN,
        xlabel=script.X_LABEL, ylabel=script.Y_LABEL,
        row_key="sample", col_key="timepoint",
    )
```

Repoint the remaining eight shims the same way, using the exact argument sets
from the Task 11–13 CLI verification commands.

- [ ] **Step 3: Delete the ten modules**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
git rm workflow/src/figure_render/correlation.py \
       workflow/src/figure_render/orientation.py \
       workflow/src/figure_render/density.py \
       workflow/src/figure_render/distribution.py \
       workflow/src/figure_render/read_counts.py \
       workflow/src/figure_render/ma_plot.py \
       workflow/src/figure_render/ma_plot_replicates.py \
       workflow/src/figure_render/dispersions.py \
       workflow/src/figure_render/curve_fitting.py \
       workflow/src/figure_render/coverage.py
```

- [ ] **Step 4: Clear stale bytecode**

The `__pycache__` holds compiled copies that can satisfy an import and mask a
missed reference.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
rm -rf workflow/src/figure_render/__pycache__ && \
ls workflow/src/figure_render/
```

Expected: exactly `__init__.py`, `_layout.py`, `_schema.py`, `composition.py`,
`curves.py`, `histogram.py`, `ma.py`, `scatter.py`, `series.py`.

- [ ] **Step 5: Run the full suite**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
timeout 900 /data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests -q 2>&1 | tail -15
```

Expected: all pass. The pre-refactor baseline was **105 passed, 2 skipped** in
~275 s; the count is now higher because Tasks 2–10 added test files. No failures,
and no more than 2 skips.

- [ ] **Step 6: Commit**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake
git add -A workflow/src/figure_render/ workflow/tests/figure_render/pixel_baseline.py
git commit -m "refactor(figure_render): delete the ten superseded figure modules

Every entrypoint now calls a grammar module. figure_render/ holds six plot
grammars plus two shared helpers, and names no biological quantity."
```

---

### Task 15: Final verification

**Files:**
- None modified — this task only verifies.

**Interfaces:**
- Consumes: everything
- Produces: nothing

- [ ] **Step 1: Assert no biological or timepoint literals remain in `src/figure_render/`**

This is the plan's headline acceptance criterion.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
grep -rniE "pbl|pbr|plus_count|minus_count|FYPOviability|baseMean|YES[0-9]|Spikein|normed_mean|genewise|insertion_density" \
  workflow/src/figure_render/*.py
```

Expected: no output. Prose mentions in docstrings that name a *superseded module*
(e.g. "Supersedes ``correlation.py`` (PBL vs PBR)") are acceptable historical
context — but no *code* may reference these. If the grep hits only docstring
lines, confirm each by eye; any hit in an expression is a failure.

- [ ] **Step 2: Assert the library-module contract holds**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
grep -rnE "^def main|^ *def parse_args|setup_logger|argparse" workflow/src/figure_render/*.py
```

Expected: no output. `src/` modules carry no CLI machinery.

- [ ] **Step 3: Confirm the timepoint whitelist fix on every project with data**

The plan's functional win. `Spore2YES6_1328` could not render before.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
PY=/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python && \
for p in HD_DIT_HAP LD_DIT_HAP Spore2YES6_1328; do
  FIT=projects/$p/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv
  LFC=projects/$p/results/14_insertion_level_depletion_analysis/LFC.tsv
  if [ -f "$FIT" ] && [ -f "$LFC" ]; then
    if $PY workflow/scripts/figures/plot_curve_fitting.py \
         -s "$FIT" -l "$LFC" -o "tmp/generalize/curves_$p" > /dev/null 2>&1; then
      echo "$p: OK"
    else
      echo "$p: FAILED"
    fi
  fi
done
```

Expected: `HD_DIT_HAP: OK`, `LD_DIT_HAP: OK`, `Spore2YES6_1328: OK`.

`LD_haploid` is expected to fail — its stats file predates the current curve
model (`um`/`lam` instead of `A`/`DR`/`DL`), a documented limitation. Confirm the
message names the missing columns:

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python \
  workflow/scripts/figures/plot_curve_fitting.py \
  -s projects/LD_haploid/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv \
  -l projects/LD_haploid/results/14_insertion_level_depletion_analysis/LFC.tsv \
  -o tmp/generalize/curves_LD_haploid 2>&1 | grep -iE "missing required|KeyError" | head -3
```

Expected: a `Missing required columns` message naming `'DR'` and `'DL'`, not a
bare `KeyError`.

- [ ] **Step 4: Confirm synthetic generalization for the seven figures without multi-project data**

Only `HD_DIT_HAP` has `18_figure_data/arc/*.tsv`, so these seven are exercised
with synthetic inputs varying timepoint count, naming and group count. This is
weaker than real data and is not claimed to be equivalent.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src')
from matplotlib import use; use('Agg')
from pathlib import Path
import numpy as np, pandas as pd
from figure_render.scatter import render_grouped_regression_figure

rng = np.random.default_rng(7)
# Six Spikein-named timepoints and two samples: neither matches YES0-4.
names = [f'Spikein{i}' for i in range(6)]
rows = []
for sample in ('A', 'B'):
    for tp in names:
        n = 40
        rows.append(pd.DataFrame({
            'sample': sample, 'timepoint': tp,
            'log10_left': rng.normal(2, 0.5, n), 'log10_right': rng.normal(2, 0.5, n),
        }))
df = pd.concat(rows, ignore_index=True)

render_grouped_regression_figure(
    df, Path('tmp/generalize/synthetic_scatter'),
    x='log10_left', y='log10_right', xlabel='left', ylabel='right',
    row_key='sample', col_key='timepoint',
)
import matplotlib.pyplot as plt
print('panels:', len(plt.gcf().get_axes()), '(expected 12)')
print('artifact:', Path('tmp/generalize/synthetic_scatter.pdf').exists())
" 2>&1 | grep -E "panels:|artifact:"
```

Expected: `panels: 12 (expected 12)` and `artifact: True`. Repeat the pattern for
the other grammars if a reviewer wants broader coverage; the grammar-level unit
tests from Tasks 4–10 already parametrize over timepoint counts and naming.

- [ ] **Step 5: Verify zero findfont warnings across all ten figures**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/src'); sys.path.insert(0,'workflow/tests/figure_render')
from matplotlib import use; use('Agg')
from pathlib import Path
from pixel_baseline import render_all_baseline_figures
got = render_all_baseline_figures(Path('tmp/final'))
print('rendered', len(got), 'figures')
" 2>&1 | tee tmp/final_render.log | grep -E "rendered" ; \
echo "findfont warnings: $(grep -c findfont tmp/final_render.log)"
```

Expected: `rendered 10 figures` and `findfont warnings: 0`.

- [ ] **Step 6: Confirm the eight unchanged figures are still byte-identical**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
/data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -c "
import sys; sys.path.insert(0,'workflow/tests/figure_render')
from pathlib import Path
from pixel_baseline import EXPECTED_TO_CHANGE, compare_png

names = ['correlation','orientation','read_counts','density','coverage',
         'dispersions','ma_plot_replicates','ma_plot','distribution','curve_fitting']
failures = []
for name in names:
    base = Path(f'tmp/pixel_baseline/{name}.review.png')
    new = Path(f'tmp/final/{name}.review.png')
    if not (base.exists() and new.exists()):
        print(f'{name}: MISSING artifact'); failures.append(name); continue
    same, diff = compare_png(base, new)
    if name in EXPECTED_TO_CHANGE:
        print(f'{name}: {\"identical\" if same else f\"differs ({diff})\"} (change approved)')
    elif same:
        print(f'{name}: IDENTICAL')
    else:
        print(f'{name}: DIFFERS ({diff})  <-- REGRESSION'); failures.append(name)
print()
print('REGRESSIONS:', failures if failures else 'none')
"
```

Expected: eight `IDENTICAL`, two `(change approved)`, and `REGRESSIONS: none`.

- [ ] **Step 7: Verify the Snakemake DAG still resolves**

No rule was touched, so the DAG must be unaffected.

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
timeout 600 snakemake -n --use-conda 2>&1 | tail -20
```

Expected: the DAG builds with no error. A `Nothing to be done` or a job listing
both pass; an exception or a `MissingInputException` does not.

- [ ] **Step 8: Run the full suite one final time**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
timeout 900 /data/a/yangyusheng/miniforge3/envs/cnsplots_figures/bin/python -m pytest \
  workflow/tests -q 2>&1 | tail -10
```

Expected: all pass, at most 2 skips (the pre-refactor baseline was 105 passed,
2 skipped).

- [ ] **Step 9: Clean up scratch artifacts and commit any harness fixes**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_snakemake && \
rm -rf tmp/check tmp/cli tmp/final tmp/generalize tmp/_findfont_check tmp/final_render.log && \
git status --short
```

`tmp/pixel_baseline/` may be kept for future refactors — it is git-ignored
either way. If `git status` shows uncommitted harness changes from Task 14
Step 2, commit them:

```bash
git add workflow/tests/figure_render/pixel_baseline.py
git commit -m "test(figures): repoint pixel baseline harness to the generic renderers"
```

---

## Verification Summary

| Check | Command | Expected |
|---|---|---|
| No biological literals in `src/figure_render/` | Task 15 Step 1 | no code hits |
| Library-module contract | Task 15 Step 2 | no CLI machinery |
| Whitelist crash fixed | Task 15 Step 3 | 3 projects OK |
| `LD_haploid` fails readably | Task 15 Step 3 | names `'DR'`, `'DL'` |
| Synthetic generalization | Task 15 Step 4 | 12 panels, artifact written |
| Zero findfont warnings | Task 15 Step 5 | `0` |
| Eight figures byte-identical | Task 15 Step 6 | `REGRESSIONS: none` |
| DAG resolves | Task 15 Step 7 | no error |
| Full suite | Task 15 Step 8 | all pass, ≤2 skips |
| Visual sign-off: correlation | Task 11 Step 9 | user approves |
| Visual sign-off: curve_fitting | Task 12 Step 9 | user approves |





