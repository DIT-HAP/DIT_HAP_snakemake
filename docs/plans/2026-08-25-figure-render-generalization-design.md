# Figure Render Generalization Design

**Date:** 2026-08-25
**Status:** Approved
**Goal:** Turn `workflow/src/figure_render/` from ten figure-specific modules into six general rendering primitives, with all figure-specific knowledge moved to `workflow/scripts/figures/`.

---

## Problem

`workflow/src/` is meant to be a library of general, reusable modules (per the
`python-script-conventions` skill and the 2026-08-19 scripts/src normalization
design). The `figure_render/` subpackage violates this: each module encodes one
specific figure rather than one kind of plot.

**The reported case.** `render_correlation_figure()` in
`workflow/src/figure_render/correlation.py` is a log-log scatter with a
regression line and correlation annotation — a shape that applies to most
correlation figures. But PBL/PBR is welded into it at every level:

- Line 45 requires `['sample', 'timepoint', 'condition', 'pbl', 'pbr']`, of
  which `condition` is demanded but never used by the module.
- Lines 53–62 filter and derive `log10_pbl` / `log10_pbr` by literal name.
- Lines 110–129 hardcode the axis labels `log$_{10}$ PBL` / `log$_{10}$ PBR`
  and the grouping keys `['sample', 'timepoint']`.

Nothing in the module can be reused for any other pair of columns.

**The same figure exists twice.** `orientation.py` renders the identical shape
(sample x timepoint grid of log-log `regplot` panels) for `plus_count` vs
`minus_count`. Its layout implementation is materially better: it computes
explicit grid sizes via `PANEL_DECORATION_PX` and carries the sample name in
the first column's ylabel, where `correlation.py` still relies on
`plt.tight_layout()`. Two copies of one figure are drifting apart.

**Verified defects, not just rigidity.** Generalizing forces these into the
open; they are reproduced against real project data in this repository.

`curve_fitting.py:104` selects timepoint columns with a literal whitelist:

```python
timepoint_cols = [col for col in lfc_sampled.columns if col in ['YES0', 'YES1', 'YES2', 'YES3', 'YES4']]
```

Measured against the LFC tables actually present in `projects/`:

| Project | LFC columns | Whitelist matches | Result |
|---|---|---|---|
| `HD_DIT_HAP` | 5 | 5 | works |
| `Spore2YES6_1328` | 7 | 5 | `ValueError: x and y must be the same size` |
| `LD_haploid` | 8 (incl. `0h`) | 5 | same, plus a separate `KeyError: 'DR'` |
| `Spikein` | 6 (`Spikein0..5`) | 0 | guaranteed failure |

So the figure only renders for the two projects whose timepoints happen to be
named `YES0`–`YES4`. Four of six projects in this repository cannot produce it.

Three further defects, all confirmed by reading the code:

- **Panel-label overflow.** `chr(65 + i)` emits `[`, `\`, `]` for i >= 26.
  `distribution.py` (15 metrics) and `curve_fitting.py` (32 curves by default)
  both use this bare form; `curve_fitting.py` is already past the boundary and
  is emitting non-letter panel labels today. Two incompatible variants of this
  logic exist across five modules — `orientation.py` and `read_counts.py` carry
  the correct overflow-handling version.
- **House style bypassed.** `curve_fitting.py:125,135` hardcode `#1f77b4` and
  `#ff7f0e`, matplotlib's defaults, ignoring the Cell palette that
  `apply_house_style()` installs.
- **Clipped y-axis.** `curve_fitting.py:147` fixes `set_ylim(-1.5, 8.5)`. LFC
  range varies per project; values outside are silently cut off.

**Cross-cutting duplication.** Panel-letter generation is rewritten in five
modules; required-column validation is rewritten in nine loaders; grid sizing
is duplicated between `orientation.py` and `read_counts.py`.

---

## Solution

Reorganize by **plot grammar** rather than by figure. Each `src/figure_render/`
module owns one kind of plot and takes every figure-specific value as an
argument. Column names, axis labels and required-column lists move into the
`workflow/scripts/figures/plot_*.py` entrypoint that owns each figure.

## Target Layout

```
workflow/src/figure_render/
  __init__.py
  _layout.py      panel_labels(), grid_panel_size(), PANEL_DECORATION_PX
  _schema.py      require_columns()
  scatter.py      render_grouped_regression_figure()   <- correlation + orientation
                  render_scatter_grid_figure()          <- density
  histogram.py    render_histogram_grid_figure()        <- distribution + read_counts
  series.py       render_series_scatter_figure()        <- dispersions
  ma.py           render_ma_figure()                    <- ma_plot + ma_plot_replicates
  curves.py       render_fitted_curves_figure()         <- curve_fitting
  composition.py  render_composition_figure()           <- coverage
```

Ten figure modules become six grammar modules plus two shared helpers. **All
ten existing modules are deleted**; each figure's specific constants move to its
entrypoint script. The mapping from old module to superseding grammar module:

| Deleted module | Superseded by |
|---|---|
| `correlation.py`, `orientation.py` | `scatter.py` (`render_grouped_regression_figure`) |
| `density.py` | `scatter.py` (`render_scatter_grid_figure`) |
| `distribution.py`, `read_counts.py` | `histogram.py` |
| `ma_plot.py`, `ma_plot_replicates.py` | `ma.py` |
| `dispersions.py` | `series.py` |
| `curve_fitting.py` | `curves.py` |
| `coverage.py` | `composition.py` |

`correlation.py` + `orientation.py` is a true near-duplicate merge: same plot
shape, same input schema, only the column names and labels differ.

The other two merges are looser and must not be oversold:

- **`ma_plot` + `ma_plot_replicates`** share the MA plot shape but not their
  inputs. `ma_plot` reads two wide tables (baseMean and LFC, one column per
  timepoint) and offers a vertical/horizontal `Orientation`; the replicate
  branch reads one long table with a `padj` column and colours points by
  significance. `ma.py` must expose both: a `point_colors` hook (default: a
  single colour) and the orientation switch. The two branches are selected by
  `has_replicates` in `depletion_scoring.smk` and never run together, so
  neither may regress in service of the other.
- **`distribution` + `read_counts`** share the histogram-grid shape, but
  `read_counts` additionally replays pre-binned counts via `weights` and draws
  a cutoff annotation plus a figure-level retention footer. `histogram.py` keeps
  a pre-binned mode alongside the raw-values mode.

`dispersions` -> `series.py`, `curve_fitting` -> `curves.py`, `coverage` ->
`composition.py` and `density` -> `render_scatter_grid_figure` are one-to-one
moves into grammar-named modules with generalized signatures.

### Loaders move to scripts

Every current loader is `read_csv` + required-column validation + a few derived
columns. The required-column list *is* figure-specific knowledge, so under the
chosen boundary it belongs in `scripts/`. `src` keeps only the generic
`require_columns()` primitive.

Consequence: `load_and_prepare_data`, `load_coverage_data`,
`load_dispersion_data`, `load_density_data`, `load_distribution_data`,
`load_cutoff_stats`, `load_fitting_stats`, `load_ma_data` and
`load_and_sample_data` all disappear from `src`. All ten test files must
repoint their imports.

---

## Module Boundaries

Every module stays a **library module** per `python-script-conventions`:
IMPORTS -> CONSTANTS -> CORE LOGIC, no `main()`, no `parse_args()`, no
`setup_logger()`.

### Renderer signature convention

Illustrated with the merged scatter/regression renderer:

```python
def render_grouped_regression_figure(
    df: pd.DataFrame,
    output_stem: Path,
    *,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    group_by: Sequence[str] = ("sample", "timepoint"),
    row_key: str | None = None,
    col_key: str | None = None,
    scatter_kws: Mapping[str, Any] | None = None,
    panel_decoration_px: int = PANEL_DECORATION_PX,
) -> None:
```

Rules applied to all six renderers:

1. **`x`, `y`, `xlabel`, `ylabel` are keyword-only with no default.** A stale
   call like `render_correlation_figure(df, stem)` must raise `TypeError`, never
   silently draw PBL/PBR labels over someone else's data.
2. **`group_by` keeps the `("sample", "timepoint")` default.** That is the
   long-format convention of the figure-data TSVs, not one figure's private
   knowledge — verified identical in `pbl_pbr_pairs.tsv` and `strand_pairs.tsv`.
3. **`row_key` / `col_key` control the grid.** Extracted from
   `orientation.py`'s behaviour: rows per sample, columns per timepoint, sample
   name in the first column's ylabel only. Default `None` falls back to the
   first two entries of `group_by`.
4. **No business defaults for labels or column names.** Anything naming a
   biological quantity lives in `scripts/`.

### cnsplots constraints the renderers must preserve

Verified against installed cnsplots 0.6.0; these silently produce wrong output
if violated.

- **`regplot` resets the axis labels it manages.** Labels must be applied
  *after* drawing, as `orientation.py:158-161` already does. The generic
  renderer keeps that ordering.
- **`regplot` drops its own `s` argument when `scatter_kws` is supplied.**
  Marker size must live inside `scatter_kws`.
- **`regplot` computes its own `r` and `P` from the columns it is given.** The
  project convention is PCC on log10 values, so callers pass explicit `log10_*`
  columns. Passing raw values reports a different statistic while looking
  equally plausible (PBL/PBR: raw r ~ -0.03 vs log-space r = 0.85).
- **`cns.scatterplot` always passes `edgecolor=None` internally.** Supplying
  `edgecolor` in its kwargs raises on the duplicate. Because `regplot` and
  `scatterplot` need different kwargs, `scatter_kws` defaults are per-renderer:
  `render_grouped_regression_figure` defaults to the regplot set
  (`facecolor="none", edgecolor="gray"`), `render_scatter_grid_figure` to the
  scatterplot set. They do not share a default.
- **Sizes are pixels, not inches** (`inches = px / 72`), and `cns.multipanel`
  is a class, not a function.

---

## Correctness Fixes

Carried as part of generalization, not as separate work.

### 1. Timepoint whitelist (fix)

Replace the `['YES0'..'YES4']` whitelist with `list(lfc.columns)`. The LFC table
columns are the full ordered timepoint set, matching the length of the
`time_points` field in the stats table — verified across all four projects that
have both files.

Add an explicit guard before rendering:

```python
if len(timepoint_cols) != len(time_points):
    raise ValueError(...)  # reports both actual lengths
```

This turns a `ValueError` surfacing from deep inside matplotlib into a
diagnosable error at the boundary.

### 2. `LD_haploid` legacy parameter names (documented limitation, not fixed)

That project's stats file carries `um` / `lam` instead of `A` / `DR` / `DL`, and
currently fails with `KeyError: 'DR'`. This is a historical data-format issue
outside the rendering layer's remit. The loader will surface it through
`require_columns()` with a readable message ("missing A, DR, DL — file may
predate the current curve model") rather than pretending to render. Explicitly
a **known limitation**, not a fix.

### 3. Panel-label overflow (fix)

Single `panel_labels()` in `_layout.py` with the overflow handling already used
by `orientation.py` and `read_counts.py`. All five call sites converge on it.

### 4. Hardcoded colours (fix, changes pixels)

`curve_fitting.py`'s `#1f77b4` / `#ff7f0e` are replaced by house-style palette
colours. This changes that figure's pixels and requires visual sign-off.

### 5. Fixed y-limits (parameterized, default unchanged)

`set_ylim(-1.5, 8.5)` becomes a parameter **defaulting to `(-1.5, 8.5)`**, so
panels stay comparable across projects. Callers may override; the default
preserves current behaviour and therefore current pixels.

---

## Verification Strategy

### Pixel comparison is exact, not tolerance-based

Confirmed empirically: rendering `render_coverage_figure` twice on identical
input produces byte-identical PNGs (1624x1953x3, max abs diff = 0). So
comparison uses `np.array_equal` with no threshold. PDFs embed timestamps and
cannot be compared directly; `pdftoppm` and `gs` are both on PATH for
rasterization if a PDF-level check is wanted. The render env already has PIL
12.3.0 and numpy 2.4.6 — nothing new to install.

### Generalization coverage is uneven — stated plainly

| Figure | Baseline (`HD_DIT_HAP`) | Real multi-project data |
|---|---|---|
| curve_fitting, distribution, ma_plot | yes | `LD_DIT_HAP`, `LD_haploid`, `Spore2YES6_1328`, `Spikein` |
| correlation, orientation, read_counts, density, coverage, ma_plot_replicates, dispersions | yes | **none** |

Only `HD_DIT_HAP` has `18_figure_data/arc/*.tsv`; the other projects have never
run that computation layer. For those seven figures, generalization is verified
with **synthetic TSVs** varying timepoint count, timepoint naming
(`Spikein0..5` style) and sample count, asserting the figure renders and the
panel count is correct. This is weaker than real data and is not claimed to be
equivalent; it does establish that `YES0`–`YES4` is no longer assumed.

### Three verification layers, each gating the next

1. **Pixel baseline.** Before any code change, render all ten figures from
   `HD_DIT_HAP` and store the PNGs. After refactoring, re-render: eight must be
   byte-identical; correlation and curve_fitting go to the user for visual
   sign-off.
2. **Generalization.** Three figures against real multi-project data; seven
   against synthetic TSVs.
3. **Regression.** `pytest workflow/tests -q` fully green,
   `snakemake -n --use-conda` resolves the DAG, and **zero `findfont`
   warnings** — the established pass criterion for figure work in this project.

### Expected pixel changes

Two figures change and need visual confirmation:

- **correlation** — adopts `orientation.py`'s grid layout (chosen deliberately;
  that implementation is more mature than `plt.tight_layout()`).
- **curve_fitting** — palette correction from fix #4 only. Fix #5 keeps the
  default y-limits, so no change from that.

The other eight are strictly byte-identical.

---

## Phased Rollout

Shared layer first, then risk ascending. **Phase 1's pixel baseline must be
captured before any code change**, or the baseline is contaminated.

**Phase 1 — Shared helpers.** Create `_layout.py` and `_schema.py` with unit
tests. Pure addition; no figure touched. Capture the pixel baseline first.

**Phase 2 — `scatter.py`.** Merge `correlation.py` and `orientation.py`. The
originally reported problem, and the one figure needing layout sign-off.

**Phase 3 — `curves.py`.** Timepoint whitelist fix plus palette correction.
Highest risk (the only reproduced crash), but Phase 2 has already proven the
pattern.

**Phase 4 — Remaining grammars.** `histogram.py`, `series.py`, `ma.py`,
`composition.py` — six figures, the most mechanical stage.

**Phase 5 — Wiring.** Ten entrypoint scripts updated, ten test files rewritten,
four superseded modules deleted.

**Invariant:** at the end of every phase all ten figures still render. There is
no intermediate state where a figure is broken.

---

## Risks & Mitigations

| Risk | Mitigation |
|---|---|
| Layout regression invisible to existing tests | Pixel baseline on all ten figures before any change |
| Seven figures have no real multi-project data | Synthetic TSVs; limitation stated rather than papered over |
| A stale call silently renders wrong labels | `x`/`y`/`xlabel`/`ylabel` keyword-only with no default -> `TypeError` |
| `regplot` label reset reintroduced | Post-draw label application encoded in the generic renderer |
| Wrong statistic from raw-value columns | Callers pass explicit `log10_*` columns; documented in the renderer docstring |
| `scatterplot` duplicate-`edgecolor` crash | Per-renderer `scatter_kws` defaults, not shared |
| Ten test files break at once | Test rewrite confined to Phase 5, after renderers are pixel-verified |

---

## Success Criteria

1. No module in `src/figure_render/` names a biological quantity (`pbl`, `pbr`,
   `plus_count`, `FYPOviability`) or a timepoint literal (`YES0`).
2. Eight figures byte-identical to baseline; correlation and curve_fitting
   visually approved.
3. `curve_fitting` renders for `Spore2YES6_1328` and `Spikein`, which it cannot
   today.
4. `pytest workflow/tests -q` green; `snakemake -n --use-conda` resolves; zero
   `findfont` warnings.
5. Panel-letter generation, required-column validation and grid sizing each
   exist exactly once.
6. Entrypoint scripts remain docstring + `parse_args()` + `main()` + figure
   constants.

---

## Out of Scope

- **No Snakefile or rule changes.** Rules keep invoking
  `python workflow/scripts/figures/plot_*.py` with the same CLI flags.
- **No env changes.** No new dependencies; PIL and numpy are already present.
- **No computation-layer changes.** `src/qc/` and `src/depletion/` untouched.
- **No fix for `LD_haploid`'s legacy `um`/`lam` stats format** — surfaced as a
  readable error, documented as a known limitation.
- **No new figures or statistical methods.** Structural change plus the five
  enumerated correctness fixes.

---

## References

- **Skill:** `python-script-conventions` (library vs standalone module layout).
- **Prior design:** `docs/plans/2026-08-19-scripts-src-normalization-design.md`
  (established `src/` as a library; this design completes `figure_render/`).
- **Prior design:** `docs/plans/2026-08-13-cnsplots-figure-refactor-design.md`
  (compute/render env split these modules live in).
- **Memory:** [[cnsplots-api-gotchas]] — the 0.6.0 behaviours encoded above.
