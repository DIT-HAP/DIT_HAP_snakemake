# Scripts/Src Normalization Design

**Date:** 2026-08-19  
**Status:** Approved  
**Goal:** Normalize `workflow/src/` as an importable library, reduce `workflow/scripts/` to thin CLI entrypoints.

---

## Problem

The current `workflow/scripts/` structure violates separation of concerns:

- **33 scripts hold both domain logic and CLI machinery** (12.6k lines total).
- **32 byte-identical copies of `setup_logger()`** across scripts (the 33rd is `setup_logging()` in `extract_genome_region.py`).
- **Repeated I/O conventions**: 14 scripts use `index_col=[0,1,2,3]` for insertion tables; no shared constant or function.
- **Two near-duplicate `update_sysID` implementations**: `src/utils.py:23` (takes Path, uses `print`) and `reference_data/extract_genome_region.py:141` (takes DataFrame, uses logger). One should be canonical.
- **`workflow/src/` is not a library**: holds `concat_tables.py` (a stray `# %%` cell script, not imported by anything), `utils.py` (imported only by `insertion_orientation_analysis.py:67`), `figures.py` (imported by 8 scripts via `sys.path.append("../../src")`), and one test file.
- **Tests are co-located with scripts** but import via `sys.path.append` — moving core logic into `src/` breaks 12 test files.

This makes scripts hard to test, reuse, and maintain. A bug fix in table-reading logic requires patching 14 files.

---

## Solution

**Restructure into a library/CLI split:**

1. `workflow/src/` becomes a clean library of domain modules (no `main()`, no `parse_args()`, no logger setup).
2. `workflow/scripts/<domain>/<name>.py` becomes docstring + `parse_args()` + `main()` only (~60–90 lines per script).
3. Shared boilerplate (logger setup, table I/O, ID resolution, figure styling) consolidates into `src/` modules.
4. Per-domain core logic (depletion computation, BAM parsing, QC metrics, figure rendering) moves into `src/<domain>/` modules.
5. Tests move to `workflow/tests/` mirroring `src/`, with a root `pyproject.toml` so `pytest` works without `sys.path` hacks.
6. Scripts import from `src/` via the existing `sys.path.append` bootstrap (no Snakefile changes, no env modifications).

---

## Target Layout

```
workflow/
  src/
    __init__.py
    logging_setup.py        setup_logger()              [replaces 33 copies]
    io_tables.py            read/write insertion tables, INSERTION_INDEX constants
    gene_metadata.py        resolve_gene_ids()          [canonical; replaces 2 copies]
    figures.py              apply_house_style, save_dual, JOURNAL_* [unchanged]
    
    depletion/
      __init__.py
      weights.py            Scheme enum, compute_weights, normalize_per_gene_timepoint
      gene_level.py         aggregate_to_gene_level
      insertion_level.py    DeSeq2 wrapper, median normalization, MA computation
      curve_fitting.py      sigmoid fns, fit_single_curve, generate_summary_stats
      imputation.py         filter_insertions, impute_missing_values, calc_statistics
      ctr_insertions.py     get_control_insertions
    
    read_processing/
      __init__.py
      bam.py                extract_read_info, process_read_pair, build_schema
      filtering.py          build_filter_mask, process_chunk, load_config_from_yaml
      insertions.py         calc_insertion_coords, create_validation_mask, count_insertions
      annotation.py         calc_codon_distances, calc_affected_residue, annotate_insertions
      concat.py             concatenate (counts + annotations)
      merge.py              merge_strand_insertions, merge_timepoints
      hard_filtering.py     apply_hard_filtering
      sequence_extraction.py  load_reference_data, extract_target_sequence
    
    qc/
      __init__.py
      density.py            load_insertion_data, filter_in_gene, calc_insertion_stats, analyze_gene_insertions
      coverage.py           load_covered_genes, load_gene_viability, compute_coverage_stats
      orientation.py        extract_strand_pairs
      correlation.py        read_tsv_file, parse_filename (PBL/PBR)
      read_counts.py        load_and_validate_data, calc_cutoff_stats, compute_binned_distribution
      mapping_stats.py      parse_log_file, extract_summary_data, create_dataframe
    
    figure_render/
      __init__.py
      ma_plot.py            load_ma_data, render_ma_figure, Orientation enum
      curve_fitting.py      sigmoid_function, load_and_sample_data, render_curve_fitting_figure
      dispersions.py        load_dispersion_data, render_dispersion_figure
      density.py            load_density_data, render_density_figure
      orientation.py        load_and_prepare_data, render_orientation_figure
      coverage.py           load_coverage_data, render_coverage_figure
      read_counts.py        load_distribution_data, load_cutoff_stats, render_distribution_figure
      correlation.py        load_and_prepare_data (PBL/PBR), render_correlation_figure
      distribution.py       load_fitting_stats, render_distribution_figure
    
    reference_data/
      __init__.py
      genome_region.py      parse_gff_data, build_intergenic_bed, select_primary_transcripts, annotate_intergenic_region_flanks
  
  tests/                   [mirrors src/ tree]
    test_logging_setup.py
    test_io_tables.py
    test_gene_metadata.py
    test_figures.py
    depletion/
      test_weights.py
      test_gene_level.py
      test_insertion_level.py
      test_curve_fitting.py
      ... (1 test file moved from scripts/depletion_scoring/)
    figure_render/
      test_ma_plot.py
      test_curve_fitting.py
      ... (10 test files moved from scripts/figures/)
  
  scripts/
    <domain>/<name>.py     [docstring + parse_args() + main() only, ~60-90 lines each]
```

**Changes to existing files:**

- `workflow/src/concat_tables.py` **deleted** (dead exploration code, not imported by anything, hardcoded `./brite_table` paths).
- `workflow/src/utils.py` **split and merged**:
  - `read_file()` → `src/io_tables.py:read_table()` (generalized).
  - `update_sysIDs()` → `src/gene_metadata.py:resolve_gene_ids()` (canonical implementation, replaces both `utils.py:23` and `extract_genome_region.py:141`).
- `workflow/src/test_figures.py` → `workflow/tests/test_figures.py`.

---

## Module Boundaries

Every `src/` module is a **library module** per the `python-script-conventions` skill:

- **IMPORTS → CONSTANTS → CORE LOGIC** only.
- **No `main()`**, no `parse_args()`, no `setup_logger()`.
- Functions take explicit arguments and return values — no reading `argparse.Namespace`, no writing files as a side effect of computation.

Where a script currently computes and writes in one function, that splits into:
- A pure compute function in `src/` (returns DataFrame or dataclass).
- A `write_*` call in the script's `main()`.

**Config dataclasses** (`PlotConfig`, `AnalysisConfig`, `InputOutputConfig`) **stay in scripts/** — they validate CLI input paths, which is a CLI concern. Tests already exercise them there.

---

## Import Mechanism

**Chosen approach:** Keep `sys.path.append` bootstrap (no PYTHONPATH, no pip install).

Every script gets the identical two-line preamble:

```python
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.append(str((SCRIPT_DIR / "../../src").resolve()))

from logging_setup import setup_logger  # noqa: E402
from io_tables import read_insertion_table  # noqa: E402
from depletion.gene_level import aggregate_to_gene_level  # noqa: E402
```

**Why this instead of PYTHONPATH or pip install:**

- No Snakefile changes, no env yml modifications, no rebuild of 11 conda envs.
- Snakemake env hashes stay unchanged (no rule cache invalidation).
- Works identically in all 11 existing envs without installation.

**Naming constraints:**

- `src/logging_setup.py` (not `logging.py`) to avoid shadowing stdlib `logging` once `src` is on `sys.path`.
- `src/io_tables.py` (not `io.py`) to avoid shadowing stdlib `io`.

---

## Test Layout

The 12 co-located `workflow/scripts/<domain>/test_*.py` files move to `workflow/tests/` mirroring `src/`:

```
workflow/tests/
  test_figures.py                      (already exists, stays)
  depletion/test_insertion_level.py    (was scripts/depletion_scoring/test_insertion_level_depletion_analysis_has_replicates.py)
  figure_render/test_ma_plot.py        (was scripts/figures/test_plot_ma_plot.py)
  ... (10 more from figures/)
```

Assertions on **domain logic** import from `src/`. Assertions on **CLI config validation** import from the `scripts/` module (via `sys.path.append("../../scripts/<domain>")`).

Add root `pyproject.toml`:

```toml
[tool.pytest.ini_options]
pythonpath = ["workflow/src", "workflow/scripts"]
testpaths = ["workflow/tests"]
```

So `pytest` works without per-file `sys.path` juggling.

---

## Verification Strategy

**Baseline → refactor → byte-compare per phase.**

For each script with usable `HD_DIT_HAP` inputs:

1. **Before refactor:** Run the script, save outputs to `tmp/refactor_baseline/<phase>/`.
2. **After refactor:** Re-run → `tmp/refactor_after/<phase>/`.
3. **Compare:**
   - TSV/CSV: `sha256sum` for byte-identity.
   - Where float formatting may differ: `pd.testing.assert_frame_equal(rtol=1e-9, atol=1e-12)`.
   - PDF/PNG: existence + non-zero size (PDFs embed timestamps, so byte-compare won't work). Compare the underlying figure-data TSVs instead where those exist.
4. **Unit tests:** `pytest workflow/tests/<phase> -q`.
5. **DAG:** `snakemake -n --use-conda`.

**Caveat:** `read_processing` scripts `parse_bam_to_tsv`, `filter_aligned_reads`, and `extract_insertion_sites` cannot be byte-verified end-to-end without re-running from FASTQ. `HD_DIT_HAP` results show `1_fastp/` empty and `3_mapped/` has 1 file. For those three:

- Baseline against whatever `6_filtered`/`7_insertions` inputs still exist.
- Where that's impossible, rely on unit tests + dry-run only and document the gap explicitly.

---

## Phased Rollout

**Phase 1: Foundation** (shared modules + all imports repointed + tests moved + dead code deleted)

- Create `src/logging_setup.py`, `src/io_tables.py`, `src/gene_metadata.py`.
- Delete `src/concat_tables.py`, split `src/utils.py` into the above.
- Repoint all 33 scripts' imports to use the new shared modules (no domain logic moved yet).
- Move 12 test files to `workflow/tests/` and add `pyproject.toml`.
- Verify: all 33 scripts still produce identical outputs.

**Phase 2: Figures** (10 scripts, 4030 lines, 10 tests)

- Create `src/figure_render/` modules.
- Extract domain logic from `workflow/scripts/figures/*.py` into `src/figure_render/`.
- Repoint figure script imports.
- Verify: compare figure-data TSVs + PDF/PNG existence.

**Phase 3: Depletion Scoring** (7 scripts, 2582 lines, 1 test)

- Create `src/depletion/` modules.
- Extract domain logic from `workflow/scripts/depletion_scoring/*.py`.
- Verify: byte-compare all outputs (stages 14–17 have usable inputs).

**Phase 4: Read Processing** (9 scripts, 3097 lines, 0 tests)

- Create `src/read_processing/` modules.
- Extract domain logic from `workflow/scripts/read_processing/*.py`.
- Verify: baseline where inputs exist; unit tests + dry-run for `parse_bam`/`filter_aligned_reads`/`extract_insertion_sites`.

**Phase 5: Quality Control** (6 scripts, 1986 lines, 0 tests)

- Create `src/qc/` modules.
- Extract domain logic from `workflow/scripts/quality_control/*.py`.
- Verify: byte-compare all QC outputs.

**Phase 6: Reference Data** (1 script, 690 lines, 0 tests)

- Create `src/reference_data/genome_region.py`.
- Extract domain logic from `workflow/scripts/reference_data/extract_genome_region.py`.
- Verify: dry-run (this script has never run for `HD_DIT_HAP`, so no baseline exists).

Each phase verified before the next.

---

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Silent numerical drift during extraction | Byte-compare outputs per phase before proceeding |
| Float formatting differences mask real bugs | Use `pd.testing.assert_frame_equal` with tight tolerances where needed |
| `parse_bam`/`filter_aligned_reads`/`extract_insertion_sites` unverifiable without full FASTQ re-run | Rely on unit tests + dry-run; document the gap explicitly |
| Refactor breaks existing tests | Move and repoint tests in Phase 1 before touching domain logic |
| Import bootstrap fails in one of the 11 envs | The bootstrap is identical to the existing `figures/` pattern that already works across envs |

---

## Success Criteria

1. All 33 scripts produce byte-identical (or float-tolerant-equal) outputs after refactor.
2. `pytest workflow/tests -q` passes.
3. `snakemake -n --use-conda` resolves the DAG without errors.
4. Every `src/` module is a library module (no `main()`, no CLI machinery).
5. No duplicated `setup_logger()` or table I/O logic.
6. Scripts are 60–90 lines of docstring + CLI + orchestration only.

---

## Out of Scope

- **No Snakefile or env changes** — the refactor is internal to `workflow/scripts` and `workflow/src`.
- **No rule renames or path changes** — rules still invoke `python workflow/scripts/<domain>/<name>.py`.
- **No pip installation** — `src/` remains a runtime-imported library.
- **No new features or algorithms** — purely structural refactor.

---

## References

- **Skill:** `python-script-conventions` (defines library vs standalone script layout).
- **Template:** `workflow/src/figures.py` (already a library module, no `main()`).
- **In-progress refactor:** [[cnsplots-figure-refactor]] — this design complements that work (compute/render split is orthogonal to scripts/src normalization).
