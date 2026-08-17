# cnsplots Figure Refactor — Design

**Date:** 2026-08-13
**Status:** Approved, pending implementation
**Scope:** All 8 figure-producing scripts + shared plotting library

## Goal

Replace the pipeline's matplotlib plotting with `cnsplots` to produce
publication-quality figures, while separating computation from rendering so a
plotting-library dependency downgrade cannot affect analysis results.

## Motivation

Three problems with the current state:

1. **Figures are poster-scale, not journal-scale.** `config/DIT_HAP.mplstyle`
   sets 7×6 inch figures, 16 pt base font, 4 pt line widths. Nothing in
   `reports/` can go to a journal without being redrawn.
2. **Plotting and computation are entangled.** `curve_fitting.py` (633 lines,
   16 threads) and `insertion_density_analysis.py` (1033 lines) both compute
   numbers and draw figures. A figure tweak forces a full recompute.
3. **Dead and broken plotting code.** `distribution_of_curve_fitting_results.py`
   has zero references in any `.smk`. `curve_fitting.py`'s fitted-curve grid is
   commented out. `insertion_orientation_analysis.py` opens `PdfPages` inside a
   per-file loop, so every file overwrites the previous PDF and only the last
   file's figures survive.

## Constraint: dependency conflict

Measured against the real Python 3.12 environment, not assumed:

| Package | `statistics_and_figure_plotting.yml` (current) | cnsplots 0.6.0 requires |
|---|---|---|
| matplotlib | 3.11.1 | 3.10.9 (`<3.11`) |
| pandas | 3.0.3 | 2.3.3 |
| numpy | 2.5.1 | 2.4.6 |

cnsplots is PyPI-only (no conda-forge package). It cannot coexist with the
current pins. Downgrading the shared environment would pull the numerical
scripts back to pandas 2.3, putting analysis correctness at risk for a
cosmetic gain — pandas 3 changed copy-on-write and string dtype behaviour.

**Decision:** a dedicated `workflow/envs/cnsplots.yml` for rendering only. The
computation environment keeps pandas 3.

## Architecture: two layers

**Computation layer** — environment `statistics_and_computation.yml`
(pandas 3). Computes every number a figure needs and writes it to an
intermediate TSV. Does not import matplotlib.

**Rendering layer** — environment `cnsplots.yml` (pandas 2.3 + cnsplots
0.6.0). Reads intermediate TSVs, draws with cnsplots, saves two artifacts.
Performs no statistical inference and recomputes no values.

The practical payoff: a figure that looks wrong is fixed by re-running a
seconds-long plotting script, not a 16-thread curve fit. The pandas downgrade
is confined to read-TSV-and-draw operations, which are insensitive to the
version differences that matter.

### Directory layout

```
workflow/
├── envs/
│   ├── statistics_and_computation.yml   # renamed from statistics_and_figure_plotting.yml
│   └── cnsplots.yml                     # new: rendering only
├── src/
│   ├── plot.py                          # deleted (superseded by figures.py)
│   └── figures.py                       # new: shared rendering library
└── scripts/
    ├── quality_control/                 # computation halves stay here
    └── figures/                         # new: 8 pure rendering scripts
```

### `workflow/src/figures.py`

A **library module** per python-script-conventions §1.1 — IMPORTS, CONSTANTS,
CORE LOGIC only; no `main()`, no `parse_args()`, no `setup_logger()`.

| Function | Responsibility |
|---|---|
| `apply_house_style()` | Ecotyper1 palette, Arial font, `pdf_fonttype=42` for Illustrator editability |
| `save_dual(stem)` | Journal artifact (PDF, 180 mm two-column) + review artifact (PNG, enlarged) in one call |
| `scatter_with_stats()` | Replaces `create_scatter_correlation_plot`; built on `cns.regplot`, keeps the PCC / R² / slope / RMSE box |

**Font note:** cnsplots defaults to Helvetica, which is absent on this machine,
emitting dozens of `findfont` warnings per figure and silently falling back.
Arial is installed and is already third in `cns.settings.font_sans_serif`.
`apply_house_style()` removes Helvetica from that tuple and sets
`panel_label_fontname = "Arial"` — warnings gone, and consistent with the
outgoing mplstyle, which also used Arial.

## Script mapping

| Original script | Computation output | Rendering script | Figure form |
|---|---|---|---|
| PBL_PBR_correlation | `pbl_pbr_pairs.tsv` | `plot_pbl_pbr_correlation.py` | `cns.regplot` per sample, log-log + stats box |
| read_count_distribution | `read_count_distribution.tsv` + `read_count_cutoff_stats.tsv` | `plot_read_count_distribution.py` | Pre-binned histogram + cutoff line |
| insertion_orientation | `strand_pairs.tsv` | `plot_insertion_orientation.py` | `cns.regplot` panel matrix (fixes overwrite bug) |
| insertion_density | *reuses* `insertion_density_analysis.tsv` | `plot_insertion_density.py` | 2×2 scatter (essentiality-coloured) + metric histograms |
| gene_coverage | `gene_coverage_stats.tsv` | `plot_gene_coverage.py` | `cns.barplot` + `cns.donutplot` |
| curve_fitting | *reuses* `*_fitting_statistics.tsv` | `plot_curve_fitting.py` | Sampled fitted-curve panel grid |
| MA plot | `ma_values.tsv` | `plot_ma.py` | `cns.scatterplot` per timepoint |
| distribution_of_curve_fitting | *reuses* `*_fitting_statistics.tsv` | `plot_fitting_distributions.py` | `cns.histplot` metric grid |

Only 4 scripts need a new intermediate table. The other 4 read tables the
pipeline already writes.

`gene_level_fitting_statistics.tsv` already carries observations (`YES0..YES4`),
fitted values (`*_fitted`), sigmoid parameters (`A`, `DR`, `DL`),
`FYPOviability`, and `DeletionLibrary_essentiality` on one row. So
`plot_curve_fitting.py` reads that single table — evaluating the sigmoid
analytically for a smooth curve, which is rendering, not refitting.

## Data contract

The interface between layers is TSV column names. All new tables are **long
(tidy)** format so cnsplots' `x=` / `y=` / `hue=` arguments take column names
with zero reshaping. Names follow python-script-conventions §3.3b snake_case.

```
pbl_pbr_pairs.tsv         sample  timepoint  condition  pbl  pbr
strand_pairs.tsv          sample  timepoint  plus_count  minus_count
read_count_distribution   sample  timepoint  bin_left  bin_right  count
read_count_cutoff_stats   sample  original_rows  rows_kept  pct_rows_kept
                                  original_counts  counts_kept  pct_counts_kept
ma_values.tsv             timepoint  a_value  m_value
```

Written uniformly with `sep="\t"`, `index=False`, `float_format="%.4f"`,
matching existing pipeline conventions.

### Size control

`raw_reads.filtered.tsv` is ~10 MB and `insertion_level_fitting_statistics.tsv`
is ~22 MB. Naive long-format conversion multiplies row count by the number of
timepoints. Two mitigations:

- **`read_count_distribution.tsv` stores binned counts, not raw values.** This
  caps size and keeps the binning decision in the computation layer; the
  rendering script draws pre-binned data rather than re-binning. Bin count is a
  computation-layer CLI parameter, default 50 (the current value).
- **`plot_curve_fitting.py` samples deterministically.** The original grid drew
  8×4 = 32 curves per page. The script exposes `--n-curves` (default 32) and
  `--seed` (default 42). A fixed seed keeps Snakemake re-runs reproducible.
  Sample size and seed are written into the figure title, so the figure is
  self-documenting.

### Error handling

Existing conventions carry over: `@logger.catch` on core functions, `Config`
dataclass `__post_init__` for path validation, `return 1` on failure paths.
The rendering layer adds an **empty-data guard** — a sample with no valid
points after filtering gets a panel with "No valid data" text instead of a
crash. This preserves behaviour already present in
`read_count_distribution_analysis.py`.

## Snakemake wiring

Each of the 8 figures splits into two rules. PBL/PBR shown; the rest are
structurally identical:

```python
rule pbl_pbr_pairs:                    # computation
    input:  expand(rules.merge_strand_insertions.output, ...)
    output: f"projects/{project_name}/results/18_figure_data/pbl_pbr_pairs.tsv"
    conda:  "../envs/statistics_and_computation.yml"

rule plot_pbl_pbr_correlation:         # rendering
    input:  rules.pbl_pbr_pairs.output
    output:
        journal = report(".../PBL_PBR_correlation_analysis.pdf",
                         caption=..., category="Quality Control",
                         labels={"name": "3. PBL-PBR Correlation Analysis",
                                 "type": "Correlation Plot", "format": "PDF"}),
        review  = ".../PBL_PBR_correlation_analysis.review.png"
    conda:  "../envs/cnsplots.yml"
```

- Intermediate tables land in `results/18_figure_data/`, continuing the
  existing `NN_name` numbering.
- `report()` decorates only the journal PDF. The review PNG stays out of the
  Snakemake report so each figure appears once.
- Existing `caption` / `category` / `labels` are preserved verbatim; report
  numbering and names do not change.
- Output paths and filenames are unchanged, so `packaging.smk`'s `RELEASE_MAP`
  and the README directory listing need no edits.

### Environment rename

`statistics_and_figure_plotting.yml` → `statistics_and_computation.yml`. After
the split it performs no plotting, so the old name misleads. 19 references
across 3 `.smk` files (`read_processing.smk` ×6, `depletion_scoring.smk` ×7,
`quality_control.smk` ×6) are updated together with the file itself.

Two historical design documents mention the old filename
(`2026-07-16-pipeline-performance-optimization-design.md`,
`optimization_plan.md`). These record what was true when written and are left
unchanged.

## Verification

The cnsplots skill forbids delivering untested plotting code, so every script
is run against real data.

1. Run each rendering script on real `HD_DIT_HAP` data (94 TSVs available);
   validate artifacts are non-empty via the skill's `validate_output.py`.
2. Render each PNG and **inspect visually**: text clipping, tick overlap,
   legend collision, panel misalignment, colour distinguishability.
3. Compare against the pre-refactor PDFs and confirm statistics (PCC, R²,
   coverage percentages) are numerically identical. Better-looking figures must
   not come with changed numbers.
4. `snakemake -n` to confirm a complete DAG with no orphaned rules.
5. Confirm the computation layer imports cleanly under pandas 3 and the
   rendering layer under pandas 2.3.

## Implementation order

1. `src/figures.py` + `cnsplots.yml` — everything else depends on these.
2. PBL/PBR end-to-end as the reference pattern; verify before proceeding.
3. Remaining 7 scripts, following the established pattern.
4. Rename the computation env across all 19 references.
5. Rewire `.smk` files; run `snakemake -n`.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Scope | All 8 figure-producing scripts | User-selected |
| Environment | Dedicated `cnsplots.yml` | Isolates pandas downgrade from analysis code |
| Layering | Computation and rendering fully separated | Keeps recompute-heavy scripts on pandas 3 |
| Figure size | Dual artifact: journal PDF + review PNG | Submission-ready without losing on-screen QC review |
| Palette | cnsplots default Ecotyper1 | Colour-blind-safe and print-validated; supersedes the mplstyle cycle |
| Env rename | `statistics_and_computation.yml` | Name matches post-split responsibility |

## Known follow-ups

- `config/DIT_HAP.mplstyle` becomes unused once all 8 scripts migrate. Left in
  place for now; removing it is a separate cleanup.
- `distribution_of_curve_fitting_results.py` was never wired into any rule.
  This refactor wires its replacement (`plot_fitting_distributions.py`) in, so
  the metric distributions become a real pipeline output for the first time.
