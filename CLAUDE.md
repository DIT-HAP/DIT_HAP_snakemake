# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project Overview

DIT-HAP (Diploid for Insertional Mutagenesis by Transposon and Haploid for
Analysis of Phenotype) is a Snakemake (v9.0+) workflow for analyzing piggyBac
transposon insertion sequencing data in *Schizosaccharomyces pombe*.

This repository is the **upstream pipeline**: it runs from raw reads through to
**gene-level depletion / curve-fitting tables** plus QC reports, packaged into
a `release/` folder per project. Downstream analysis (enrichment, clustering,
ML, thesis figures) lives in a separate project (`DIT_HAP_pipeline`) and
consumes the packaged `release/` tables produced here. Do not add downstream
notebooks or analysis to this repo.

The repository is **multi-project**: each experiment is self-contained under
`projects/{project_name}/` (config, sample sheet, results, reports, release).
Only one project is active per Snakemake invocation, selected in `Snakefile`.

Python 3.12 is the scripting language.

## Common Commands

```bash
# Dry run
snakemake -n --use-conda

# Full run (active project is set in Snakefile via the project = "..." line)
snakemake --use-conda --cores 16

# Full run + snakemake report, in one command (see "Running a full project +
# report" below for why --report-after-run is required instead of a rule)
snakemake --use-conda --cores 16 \
    --report projects/<project_name>/reports/snakemake_report/report.zip \
    --report-after-run -- all

# Package this project's release/ folder (run after the full run above)
snakemake --use-conda --cores 16 -- package_release

# Create conda environments only
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only

# Lint
snakemake --lint
```

Note: in this environment the `snakemake` CLI is only available inside the mamba
`snakemake` environment, not on the base PATH.

## Configuration

The active project and working directory are set in `Snakefile`:
- `project = "..."` near the top — exactly one line must be uncommented; edit
  to switch experiments (the rest stay commented as a menu of known projects).
- `workdir:` — set to this repository's path.

Each project's own config lives at `projects/{project_name}/config/config.yaml`
(HD/LD density, generationRAW/PLUS1, haploid/diploid, spikein, etc.), validated
against `workflow/schemas/config.schema.yaml` at load time. `project_name` in
the config must equal the `projects/{project_name}/` folder name (asserted in
`Snakefile`). Each config references a sample sheet (TSV with `Sample`,
`Timepoint`, `Condition`, `read1`, `read2`) via `sample_sheet:`, validated
against `workflow/schemas/samples.schema.yaml`. Shared, non-project-specific
config (`multiqc_config.yml`, `DIT_HAP.mplstyle`, templates) lives in the
top-level `config/`.

Two config flags gate which rules run for a given project:
- `use_DEseq2_for_biological_replicates` — DESeq2 vs. no-replicate LFC branch
  for insertion-level depletion (`depletion_scoring.smk`).
- `time_points` — required for curve fitting (steps 15-17) and, transitively,
  gene-level depletion in the no-replicate branch. Leave unset for a QC-only /
  normalization project (e.g. a spike-in run); `packaging.smk` skips
  gene-level release targets accordingly.

## Architecture

### Rule modules (`workflow/rules/`)
- `reference_data.smk` — Download PomBase genome/annotation; index (samtools faidx, BWA-MEM2); extract genome regions.
- `read_processing.smk` — fastp → PBL/PBR junction classification → BWA mapping → BAM→TSV parsing → aligned-read filtering → insertion-site extraction → strand merge → timepoint concatenation → feature annotation → concatenation → hard filtering.
- `depletion_scoring.smk` — Control-insertion selection → insertion-level depletion (DESeq2 when `use_DEseq2_for_biological_replicates`, else no-replicate LFC) → insertion-level curve fitting → gene-level aggregation → gene-level curve fitting. Steps 15-17 require `time_points`.
- `quality_control.smk` — MultiQC, mapping/filtering statistics, PBL/PBR correlation, read-count distribution, insertion orientation, insertion density, gene coverage.
- `packaging.smk` — `package_release`: copies this project's depletion/curve-fitting tables (declared outputs + undeclared side-outputs, see comments in the file) into `release/`, branch-aware on `time_points` like `depletion_scoring.smk`. `package_reference`: separately archives shared PomBase reference data.

Included in `Snakefile` in the order above (`reference_data` → `read_processing`
→ `depletion_scoring` → `quality_control` → `packaging`).
Wildcards `sample`, `timepoint`, `condition` are constrained from the sample
sheet.

### Running a full project + report
```bash
snakemake --use-conda --cores 16 \
    --report projects/{project_name}/reports/snakemake_report/report.zip \
    --report-after-run -- all
```
Report generation cannot be a plain rule (`shell: "snakemake --report ..."`
would deadlock on the directory lock the parent invocation still holds).
`--report-after-run` is Snakemake's own mechanism: it runs the workflow,
releases the lock, then builds the report — all in one process. Run
`package_release` as a separate step afterward to populate `release/`.

### Output structure
- `projects/{project_name}/results/1_fastp/` … `17_gene_level_curve_fitting/` — numbered pipeline steps; 15-17 only exist when `time_points` is set.
- `projects/{project_name}/release/{insertion_level,gene_level}/` — downstream interface (packaged by `package_release`).
- `projects/{project_name}/reports/` — QC reports (PDF/HTML/TSV) plus `snakemake_report/report.zip`.
- `projects/{project_name}/logs/{reference_data,read_processing,depletion_scoring,quality_control}/` — per-rule logs.
- `resources/pombase_data/{release_version}/` — auto-downloaded reference data, shared across projects.

### Scripts (`workflow/scripts/`)
Organized by module: `reference_data/`, `read_processing/`, `depletion_scoring/`,
`quality_control/`. Each Snakemake rule calls its script via `shell:` with an
explicit CLI (`-i/-o/...`); there are no inline `run:` blocks.

## Python Script Conventions

All Python scripts follow the **`python-script-conventions` skill** (Modern
Python 3.12+). Load that skill before creating or editing any `.py` file. Key
points: 7-section layout (IMPORTS → DECORATORS → CONSTANTS/ENUMS → CONFIG →
LOGGING → CORE → MAIN), `@dataclass(kw_only=True, slots=True, frozen=True)` for
config objects, native generics (`list`/`dict`, `X | None`), `StrEnum` for
grouped string markers, loguru with `@logger.catch` (no `print`), `pathlib.Path`,
single-line function docstrings, and `parse_args()` + `sys.exit(main())`.

When modifying a pipeline script, preserve its CLI contract exactly — the
Snakemake rule that calls it depends on the flag names and arities.
