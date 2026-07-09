# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project Overview

DIT-HAP (Diploid for Insertional Mutagenesis by Transposon and Haploid for
Analysis of Phenotype) is a Snakemake (v8.0+) workflow for analyzing piggyBac
transposon insertion sequencing data in *Schizosaccharomyces pombe*.

This repository is the **upstream pipeline**: it runs from raw reads through to
**gene-level depletion / curve-fitting tables** plus QC reports. Downstream
analysis (enrichment, clustering, ML, thesis figures) lives in a separate
project (`DIT_HAP_pipeline`) and consumes the gene-level tables produced here.
Do not add downstream notebooks or analysis to this repo.

Python 3.12 is the scripting language.

## Common Commands

```bash
# Dry run
snakemake -n --use-conda

# Full run (active config is set in Snakefile)
snakemake --use-conda --cores 16

# Run with a specific config
snakemake --configfile config/config_HD_generationRAW.yaml --use-conda --cores 16

# Create conda environments only
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only

# Lint
snakemake --lint
```

Note: in this environment the `snakemake` CLI is only available inside the mamba
`snakemake` environment, not on the base PATH.

## Configuration

The active config file and working directory are set in `Snakefile`:
- `snakemake_config_file` near the top — edit to switch experiments.
- `workdir:` — set to this repository's path.

Config files live in `config/` (HD/LD density, generationRAW/PLUS1, haploid/
diploid, spikein, etc.). Each references a sample sheet (TSV with `Sample`,
`Timepoint`, `Condition`, `read1`, `read2`).

## Architecture

### Rule modules (`workflow/rules/`)
- `reference_data.smk` — Download PomBase genome/annotation; index (samtools faidx, BWA-MEM2); extract genome regions.
- `read_processing.smk` — fastp → PBL/PBR junction classification → BWA mapping → BAM→TSV parsing → aligned-read filtering → insertion-site extraction → strand merge → timepoint concatenation → feature annotation → concatenation → hard filtering.
- `depletion_scoring.smk` — Control-insertion selection → insertion-level depletion (DESeq2 when `use_DEseq2_for_biological_replicates`, else no-replicate LFC) → insertion-level curve fitting → gene-level aggregation → gene-level curve fitting.
- `quality_control.smk` — MultiQC, mapping/filtering statistics, PBL/PBR correlation, read-count distribution, insertion orientation, insertion density, gene coverage.

Included in `Snakefile` in the order above. Wildcards `sample`, `timepoint`,
`condition` are constrained from the sample sheet.

### Output structure
- `results/{project_name}/1_fastp/` … `17_gene_level_curve_fitting/` — numbered pipeline steps.
- Downstream interface: `16_gene_level_depletion_analysis/` and `17_gene_level_curve_fitting/`.
- `reports/{project_name}/` — QC reports (PDF/HTML/TSV).
- `logs/{project_name}/{reference_data,read_processing,depletion_scoring,quality_control}/` — per-rule logs.
- `resources/pombase_data/{release_version}/` — auto-downloaded reference data.

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
