# Unified Snakemake Report — Design

Date: 2026-07-28
Status: Approved, ready for implementation planning

## Goal

Produce one Snakemake HTML/zip report (via the existing `--report
--report-after-run` mechanism) that:

- Has a clear entry point summarizing which pipeline branches ran and why.
- Automatically surfaces the config parameters, software versions, and code
  version behind the run, in both a quick-glance summary and a searchable
  detailed table.
- Orders its sections by pipeline stage instead of alphabetically, so reading
  top-to-bottom follows the DAG.

## Scope

This design covers **only the reporting layer**: a new `reporting.smk`
module, its scripts, and the `rule all` target list. It does **not** change
any computation in `reference_data.smk`, `read_processing.smk`,
`depletion_scoring.smk`, or the metrics those scripts compute. Read
processing's own outputs remain without `report()` annotations, as they are
today — this design does not add any.

## 1. Category structure

Snakemake orders report categories by dictionary-key sort (alphabetical), not
by DAG execution order or declaration order (confirmed against snakemake
source: `render_categories` / `dictsort` in the Jinja template). To make the
report read in pipeline order, every category name gets a numeric prefix.

Final category list, top to bottom:

| Category | Content | Origin |
|---|---|---|
| `0. Overview` | Summary HTML: pipeline stage diagram + branch-decision text + key numbers | New |
| `1. Parameters & Versions` | 3 sub-tables: Config Parameters / Software Versions / Code Version | New |
| `2. Read Processing QC` | MultiQC (fastp/FastQC/samtools/Picard), mapping_filtering_statistics, PBL/PBR correlation, read count distribution, insertion orientation | Existing, renamed + moved out of old "Quality Control" |
| `3. Insertion-level results` | Insertion-level LFC/padj, insertion-level curve fitting statistics | Existing, unchanged |
| `4. Gene-level results` | Gene-level LFC, gene-level curve fitting statistics | Existing, unchanged |
| `5. Depletion QC` | Insertion density analysis, gene coverage analysis | Existing, renamed + moved out of old "Quality Control" (these depend on post-depletion outputs, hence placed after 3/4) |

The old flat "Quality Control" category is split across `2.` and `5.`
because its members depend on different pipeline stages (some need only
`hard_filtering` output, others need `14_insertion_level_depletion_analysis`
output).

## 2. Overview page (`0. Overview`)

Content, top to bottom:

1. **Pipeline stage diagram** — a horizontal row of stage boxes: `Reference
   Data → Read Processing → Depletion Scoring (branch) → Curve Fitting &
   Gene-level (conditional) → QC → Packaging`. Branch points are
   color-coded: the path actually taken is solid, skipped paths (e.g. the
   no-replicate branch, or the gene-level steps when `time_points` is unset)
   are greyed out.
2. **Branch decision text** — 3-4 sentences directly quoting config fields
   and their values, e.g.: "`use_DEseq2_for_biological_replicates = True`,
   this run uses the DESeq2 has-replicates branch: insertion-level depletion
   is preceded by imputation and produces padj rather than only LFC.
   `time_points` is set (5 timepoints), so gene-level curve fitting (steps
   15-17) is included."
3. **Key number cards** — 4-6 small cards, config-derived only (no sample
   manifest details): number of timepoints, `hard_filtering_cutoff` value,
   PomBase reference version, active branch name(s).

**Generation:** new script `workflow/scripts/reporting/generate_overview.py`,
following the repo's standard script conventions (7-section layout, frozen
dataclass config, loguru, `parse_args()` + `sys.exit(main())`). CLI: `-c/--config`
(path to the project's already-validated `config.yaml`) and `-o/--output`
(output HTML path). The script parses config fields, determines the active
branches, and renders a self-contained single-file HTML (inline CSS, no
external resources) — mirroring how `multiqc_preprocessing`'s HTML output is
marked as a `report()` target directly, with no `.rst` caption. (Snakemake's
`report()` caption mechanism only accepts `.rst`, not `.md` — confirmed
against the snakemake source before choosing this HTML-output approach over a
caption-based one.)

## 3. Parameters & Versions detailed tables (`1. Parameters & Versions`)

Three independent datavzrd tables, following the existing
"generate TSV → datavzrd rule" pattern used by `gene_level_LFC` and
`mapping_filtering_statistics`.

**Table 1 — Config Parameters** (`config_parameters.tsv`)
Columns: `Category | Parameter | Value | Description`. `Category` groups
`config.yaml` fields: Reference (`Pombase_release_version`), Preprocessing
(adapter sequences, `chunk_size`), Filtering (mapq/ncigar/nm thresholds,
`hard_filtering_cutoff`), Depletion (`use_DEseq2_for_biological_replicates`,
`initial_time_point`, `final_time_point`, `time_points`), Timepoint merging
(`merge_similar_timepoints` and related fields). Nested fields (e.g.
`aligned_read_filtering.read_1_filtering.mapq_threshold`) are flattened to a
dot-path string in `Parameter`. `Description` comes from a hardcoded
field-description dict in the script — not parsed from config.yaml comments
(unreliable to parse).

**Table 2 — Software Versions** (`software_versions.tsv`)
Columns: `Conda Env | Package | Version`. Script iterates
`workflow/envs/*.yml`, parses each file's `dependencies` list, and includes
only packages with an explicit pinned version (containing `=`); unpinned
packages are skipped (not shown as "unknown", to avoid a misleading entry).

**Table 3 — Code Version** (`code_version.tsv`)
Single-row table: `Git Commit | Branch | Commit Date | Dirty`. `Dirty` flags
whether `git status --porcelain` is non-empty. This rule is a plain `shell:`
block running git commands directly — no Python script needed.

Six new rules total (3× generate + 3× datavzrd), all in `reporting.smk`.

## 4. `rule all` and module wiring

**`Snakefile`**: uncomment and replace `rule all`'s `input:` with a
branch-aware target list ordered by category, using the same
conditional-list pattern already used in `packaging.smk`'s `RELEASE_MAP`:

```python
rule all:
    input:
        # 0. Overview + 1. Parameters & Versions
        f"projects/{project_name}/reports/overview/overview.html",
        f"projects/{project_name}/reports/parameters/datavzrd_config_parameters",
        f"projects/{project_name}/reports/parameters/datavzrd_software_versions",
        f"projects/{project_name}/reports/parameters/datavzrd_code_version",
        # 2. Read Processing QC
        f"projects/{project_name}/reports/multiqc/quality_control_multiqc_report.html",
        f"projects/{project_name}/reports/mapping_filtering_statistics/datavzrd_mapping_filtering_statistics",
        f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
        f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis.pdf",
        f"projects/{project_name}/reports/insertion_orientation_analysis/insertion_orientation_analysis.pdf",
        # 3. Insertion-level + 4. Gene-level (gene-level only when time_points is set)
        f"projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        *([f"projects/{project_name}/results/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
           f"projects/{project_name}/results/16_gene_level_depletion_analysis/datavzrd_gene_level_LFC",
           f"projects/{project_name}/results/17_gene_level_curve_fitting/datavzrd_gene_level_curve_fitting"]
          if config.get("time_points") else []),
        # 5. Depletion QC
        f"projects/{project_name}/reports/insertion_density_analysis/datavzrd_insertion_density_analysis",
        *([f"projects/{project_name}/reports/gene_coverage_analysis/gene_coverage_analysis.pdf"]
          if config.get("time_points") else []),
```

The old commented-out candidate target lines are deleted, not kept as
comments.

**`reporting.smk`**: `include`d in `Snakefile` after `quality_control.smk`
and before `packaging.smk`. Its rules depend only on `config`,
`workflow/envs/*.yml`, and git repo state — never on `results/` — so they run
first and never block any downstream rule.

## 5. Error handling

- `generate_overview.py` / the config-parameters extraction script read
  `config`, which is already validated against
  `workflow/schemas/config.schema.yaml` at `Snakefile` load time. Only
  genuinely optional fields (e.g. `time_points`, unset for QC-only projects)
  need defensive `config.get(...)` handling.
- The software-versions script iterates `workflow/envs/*.yml`; a single
  unparsable file is logged as a warning and skipped rather than aborting the
  whole run. Unpinned packages are skipped by design (see §3), not an error
  case.
- Per this repo's established convention (`docs/optimization_plan.md`):
  computation functions are not wrapped in `@logger.catch` (it silently
  swallows exceptions and returns `None`); only I/O boundaries use it.
- `code_version` is a plain `shell:` rule (`git rev-parse`, `git status
  --porcelain`); failures propagate as a normal non-zero exit, no extra
  fallback needed since the repo is guaranteed to be a git repo.

## 6. Validation plan

Test project: **`LD_DIT_HAP`** (DESeq2 has-replicates branch, `time_points`
set, `results/` and `release/` already fully populated from a prior run).
Cross-branch dry-run comparison uses **`Spikein`**
(`use_DEseq2_for_biological_replicates: False`, `time_points` unset).

1. `snakemake -n --use-conda` dry-run on both `LD_DIT_HAP` and `Spikein` to
   confirm `rule all`'s conditional targets resolve correctly under both
   branch combinations with no wildcard/dependency errors.
2. Run only the reporting-module targets (no dependency on `results/`, so
   this is fast) against `LD_DIT_HAP`'s config; manually verify
   `overview.html`'s text/numbers and the three parameter tables' contents
   match `config.yaml`, `workflow/envs/*.yml`, and `git log` ground truth.
3. Run the full `rule all` + `--report --report-after-run` against
   `LD_DIT_HAP` (already fully computed); open the resulting `report.zip` and
   confirm all 6 categories appear in the intended order, the Overview page
   renders correctly, and the 3 new datavzrd tables open correctly.
4. Confirm `package_release` is unaffected (`reporting.smk` never touches
   `release/` or `RELEASE_MAP`).

## Out of scope (explicitly deferred)

- Any change to `reference_data.smk` / `read_processing.smk` computation
  logic or their lack of `report()` annotations.
- Sample manifest / sample-sheet display in the report.
- Any QC metric, flag, PASS/FAIL gating, or gene/insertion-level QC column
  work (raised in a separate, earlier discussion — not part of this design).
