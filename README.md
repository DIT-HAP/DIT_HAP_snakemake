# DIT-HAP snakemake

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

**DIT-HAP snakemake** is a Snakemake workflow for the **DIT-HAP** project
(Diploid for Insertional Mutagenesis by Transposon and Haploid for Analysis of
Phenotype), analyzing piggyBac transposon insertion sequencing data in
*Schizosaccharomyces pombe*.

This repository is the **upstream half** of the DIT-HAP project. It takes raw
paired-end sequencing reads and produces **gene-level depletion / curve-fitting
tables** plus quality-control reports. Downstream analysis (enrichment,
clustering, machine learning, thesis figures) lives in a separate project
(`DIT_HAP_pipeline`) and consumes the packaged `release/` tables produced here
as its interface.

The repository is **multi-project**: each experiment lives self-contained under
`projects/{project_name}/` (own config, sample sheet, results, reports, release),
sharing code and reference data at the repo root. Only one project is active
per Snakemake invocation (selected in `Snakefile`); switch projects and re-run
to process another.

## What this pipeline does

1. Downloads the *S. pombe* reference genome and annotations from PomBase.
2. Preprocesses reads (fastp), classifies piggyBac PBL/PBR junctions (cutadapt),
   maps with BWA-MEM2, and extracts insertion sites.
3. Annotates insertions with genomic features and concatenates timepoints.
4. Scores depletion at the insertion level (DESeq2 with replicates, or a
   no-replicate log-fold-change path), fits depletion curves, and aggregates to
   the **gene level**.
5. Generates QC reports (MultiQC, mapping/filtering statistics, PBL/PBR
   correlation, read-count distribution, insertion orientation, insertion
   density, gene coverage) and a self-contained Snakemake HTML/zip report.
6. Packages the project's depletion/curve-fitting tables into a stable
   `release/` folder for downstream consumption.

Gene-level analysis (steps 15-17) requires `time_points` to be set in the
project's config; projects without it (e.g. a spike-in normalization run) stop
at QC + insertion-level depletion — see [Configuration](#configuration).

The interface to downstream analysis is the packaged release folder:
`projects/{project_name}/release/insertion_level/` and
`projects/{project_name}/release/gene_level/` (see
[Packaging & release](#packaging--release)).

## Requirements

- Snakemake ≥ 9.0
- Conda / Mamba (each tool runs in its own environment under `workflow/envs/`)
- Python 3.12

## Quick start

The active project is selected in `Snakefile` — edit the `project = "..."` line
near the top (all other project lines are commented out; exactly one must be
active) and re-run. Each project's config lives at
`projects/{project_name}/config/config.yaml`, validated against
`workflow/schemas/config.schema.yaml` on load.

```bash
# Dry run
snakemake -n --use-conda

# Full run (active project's default target — see Snakefile's rule all)
snakemake --use-conda --cores 16

# Full run + Snakemake HTML/zip report, in one command (see Snakemake report)
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

Sample sheets (TSV with `Sample`, `Timepoint`, `Condition`, `read1`, `read2`)
are referenced from each project's config via `sample_sheet:`.

### Adding a new project

1. Copy `config/config.template.yaml` to `projects/<name>/config/config.yaml`
   and `config/sample_sheet.template.tsv` to
   `projects/<name>/config/sample_sheet.tsv`; fill in both. `project_name` in
   the config must equal the `<name>` folder (`Snakefile` asserts this).
2. Uncomment `project = "<name>"` in `Snakefile` (and comment out whichever was
   active).
3. Run the full-run + report and `package_release` commands shown above.

## Architecture

Rule modules (`workflow/rules/`), included by `Snakefile` in this order:

| Module                    | Responsibility                                                                                                                                 |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `reference_data.smk`    | Download PomBase genome/annotation, index (samtools/BWA), extract genome regions                                                               |
| `read_processing.smk`   | fastp → junction classification → mapping → BAM parsing → filtering → insertion sites → annotation → concatenation → hard filtering    |
| `depletion_scoring.smk` | Control selection → insertion-level depletion (DESeq2 or no-replicate) → curve fitting → gene-level aggregation → gene-level curve fitting |
| `quality_control.smk`   | MultiQC, mapping/filtering statistics, PBL/PBR correlation, read-count distribution, insertion orientation, insertion density, gene coverage   |
| `packaging.smk`         | Copies each project's depletion/curve-fitting tables into `release/`; separately archives shared PomBase reference data                       |

`depletion_scoring.smk` branches on two config flags:
- `use_DEseq2_for_biological_replicates` — DESeq2 (has replicates) vs. a
  no-replicate log-fold-change path for insertion-level depletion.
- `time_points` (list of numeric hours/generations, one per timepoint column)
  — required for curve fitting and everything downstream of it (steps 15-17,
  including gene-level depletion, since its weights transitively depend on
  curve fitting in the no-replicate branch). Projects without it stop at
  insertion-level depletion (step 14) and QC; `packaging.smk` is branch-aware
  and skips the gene-level release targets accordingly.

Analysis scripts live in `workflow/scripts/{reference_data,read_processing,depletion_scoring,quality_control}/`
and follow the `python-script-conventions` standard (Modern Python 3.12+:
7-section layout, frozen dataclasses, loguru, native generics).

## Output structure

```
projects/{project_name}/
├── config/                                 # project_name's config.yaml + sample_sheet.tsv
├── results/
│   ├── 1_fastp/ ... 9_concatenated/            # read processing (per-sample)
│   ├── 10_annotated/ ... 13_filtered/          # annotation, concatenation, hard filtering
│   ├── 14_insertion_level_depletion_analysis/  # insertion-level LFC
│   ├── 15_insertion_level_curve_fitting/       # requires time_points
│   ├── 16_gene_level_depletion_analysis/       # requires time_points
│   └── 17_gene_level_curve_fitting/            # requires time_points
├── reports/
│   ├── multiqc/  fastqc/  samtools_mapping_statistics/  picard_insert_size/
│   ├── mapping_filtering_statistics/  read_count_distribution_analysis/
│   ├── insertion_orientation_analysis/  insertion_density_analysis/  gene_coverage_analysis/
│   └── snakemake_report/report.zip         # ← self-contained Snakemake report (see below)
├── release/                                # ← downstream interface (see Packaging & release)
│   ├── insertion_level/
│   └── gene_level/
└── logs/{reference_data,read_processing,depletion_scoring,quality_control}/

resources/pombase_data/{release_version}/   # shared across projects, auto-downloaded
```

## Snakemake report

Every QC and depletion/curve-fitting output relevant to the active project is
`report()`-decorated, so a plain full run followed by report generation
captures everything — MultiQC, all QC PDFs/datavzrd tables, insertion-/
gene-level depletion and curve-fitting tables. Generate the report with
`--report` + `--report-after-run` (see Quick start) rather than a bare
`--report` call: a Snakemake rule cannot shell out to `snakemake --report` on
its own workdir since the parent process still holds the directory lock;
`--report-after-run` builds the report after the run completes and the lock
releases, all in one command.

## Packaging & release

`package_release` (in `packaging.smk`) copies each project's depletion and
curve-fitting tables — including side-outputs written next to a rule's
declared output (`baseMean.tsv`, `fitting_*.tsv`, `transformed_weights.tsv`,
`imputation_statistics.tsv`) that aren't independently reachable by Snakemake
— into `projects/{project_name}/release/{insertion_level,gene_level}/` under
stable names. This is the contract downstream projects (e.g. the
`DIT_HAP_streamlit` app) consume; `annotations.tsv` is gzip-compressed for
that consumer. Run it as a separate step after the main pipeline (`all`)
finishes, as shown in Quick start. `package_reference` separately archives the
shared PomBase reference data (`resources/pombase_data/{release_version}/`)
into a `.tar.gz`, since it isn't tied to any single project's run — invoke it
directly when needed.
