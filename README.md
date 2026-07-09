# DIT-HAP (Upstream Pipeline)

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

**DIT-HAP** (Diploid for Insertional Mutagenesis by Transposon and Haploid for
Analysis of Phenotype) is a Snakemake workflow for analyzing piggyBac transposon
insertion sequencing data in *Schizosaccharomyces pombe*.

This repository is the **upstream half** of the DIT-HAP project. It takes raw
paired-end sequencing reads and produces **gene-level depletion / curve-fitting
tables** plus quality-control reports. Downstream analysis (enrichment,
clustering, machine learning, thesis figures) lives in a separate project and
consumes the gene-level tables produced here as its interface.

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
   density, gene coverage).

The interface to downstream analysis is the gene-level output:
`results/{project_name}/16_gene_level_depletion_analysis/` and
`results/{project_name}/17_gene_level_curve_fitting/`.

## Requirements

- Snakemake ≥ 8.0
- Conda / Mamba (each tool runs in its own environment under `workflow/envs/`)
- Python 3.12

## Quick start

The active config and working directory are set in `Snakefile`. Edit
`snakemake_config_file` (near the top) to switch experiments, then:

```bash
# Dry run
snakemake -n --use-conda

# Full run
snakemake --use-conda --cores 16

# Override the hardcoded config
snakemake --configfile config/config_HD_generationRAW.yaml --use-conda --cores 16

# Create conda environments only
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only

# Lint
snakemake --lint
```

Sample sheets (TSV with `Sample`, `Timepoint`, `Condition`, `read1`, `read2`)
are referenced from the config file.

## Architecture

Rule modules (`workflow/rules/`), included in this order:

| Module | Responsibility |
|--------|----------------|
| `reference_data.smk` | Download PomBase genome/annotation, index (samtools/BWA), extract genome regions |
| `read_processing.smk` | fastp → junction classification → mapping → BAM parsing → filtering → insertion sites → annotation → concatenation → hard filtering |
| `depletion_scoring.smk` | Control selection → insertion-level depletion (DESeq2 or no-replicate) → curve fitting → gene-level aggregation → gene-level curve fitting |
| `quality_control.smk` | MultiQC, mapping/filtering statistics, PBL/PBR correlation, read-count distribution, insertion orientation, insertion density, gene coverage |

Analysis scripts live in `workflow/scripts/{reference_data,read_processing,depletion_scoring,quality_control}/`
and follow the `python-script-conventions` standard (Modern Python 3.12+:
7-section layout, frozen dataclasses, loguru, native generics).

## Output structure

```
results/{project_name}/
├── 1_fastp/ ... 9_concatenated/            # read processing (per-sample)
├── 10_annotated/ ... 13_filtered/          # annotation, concatenation, hard filtering
├── 14_insertion_level_depletion_analysis/  # insertion-level LFC
├── 15_insertion_level_curve_fitting/
├── 16_gene_level_depletion_analysis/       # ← gene-level DR/DL tables (downstream interface)
└── 17_gene_level_curve_fitting/            # ← gene-level curve fitting (downstream interface)

reports/{project_name}/
├── multiqc/  fastp/  fastqc/  junction_classification/
├── mapping_filtering_statistics/  samtools_mapping_statistics/  picard_insert_size/
├── read_count_distribution_analysis/  insertion_orientation_analysis/
├── insertion_density_analysis/  gene_coverage_analysis/
```

## License

MIT — see [LICENSE](LICENSE).
