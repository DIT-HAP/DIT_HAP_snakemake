# =============================================================================
# packaging.smk — Package a project's results into a clean, distributable
# `release/` folder for downstream consumers (e.g. the DIT_HAP_streamlit app).
#
# release/ holds ONLY this project's own depletion results, under stable names
# that form the upstream<->downstream contract. Shared reference data
# (pombase_data) is packaged separately by package_reference below.
#
# package_release is run as a separate step after the main pipeline (`all`)
# and, if desired, the Snakemake HTML/zip report — it is not folded into any
# aggregate target. package_reference stays standalone too, since it archives
# reference data shared across projects, not any single project's results.
#
# IMPORTANT — side-outputs: several source files (baseMean.tsv, normed_counts.tsv,
# insertion_level_statistics.tsv, fitting_*.tsv, transformed_weights.tsv,
# imputation_statistics.tsv) are written by scripts NEXT TO a rule's declared
# output but are NOT themselves declared in any rule's `output:`. Snakemake
# cannot build an undeclared file on demand, so each release target depends on
# its rule's DECLARED sibling ("anchor") to drive the DAG, then copies the
# real source ("src") in the shell.
# =============================================================================

import re

_RESULTS = f"projects/{project_name}/results"
_RELEASE = f"projects/{project_name}/release"

# release-relative target -> (declared anchor, actual source), both results-
# relative. anchor == src when the source is itself a declared rule output.
# Always available (insertion-level depletion doesn't need time_points).
RELEASE_MAP = {
    "insertion_level/raw_reads.tsv": (
        "12_concatenated/raw_reads.tsv",
        "12_concatenated/raw_reads.tsv",
    ),
    "insertion_level/LFC.tsv": (
        "14_insertion_level_depletion_analysis/LFC.tsv",
        "14_insertion_level_depletion_analysis/LFC.tsv",
    ),
    "insertion_level/baseMean.tsv": (
        "14_insertion_level_depletion_analysis/LFC.tsv",
        "14_insertion_level_depletion_analysis/baseMean.tsv",
    ),
    "insertion_level/normed_counts.tsv": (
        "14_insertion_level_depletion_analysis/LFC.tsv",
        "14_insertion_level_depletion_analysis/normed_counts.tsv",
    ),
}

# DESeq2-replicates branch only: imputed raw reads (declared output of
# impute_missing_values_using_FR) and the padj / combined-statistics
# side-outputs that only insertion_level_depletion_analysis_has_replicates
# writes (the no-replicates branch never produces them).
if config.get("use_DEseq2_for_biological_replicates", False):
    RELEASE_MAP |= {
        "insertion_level/imputed_raw_reads.tsv": (
            "13_filtered/imputed_raw_reads.tsv",
            "13_filtered/imputed_raw_reads.tsv",
        ),
        "insertion_level/padj.tsv": (
            "14_insertion_level_depletion_analysis/padj.tsv",
            "14_insertion_level_depletion_analysis/padj.tsv",
        ),
        "insertion_level/insertion_level_statistics.tsv": (
            "14_insertion_level_depletion_analysis/padj.tsv",
            "14_insertion_level_depletion_analysis/insertion_level_statistics.tsv",
        ),
    }

# Curve-fitting / gene-level outputs (15-17) need config["time_points"] — see
# insertion_level_curve_fitting and gene_level_depletion_analysis's weights_path
# branch in depletion_scoring.smk. Skip these entirely when time_points is
# absent (e.g. Spikein, run QC-only).
if config.get("time_points"):
    RELEASE_MAP |= {
        "insertion_level/insertion_level_fitting_statistics.tsv": (
            "15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
            "15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
        ),
        "insertion_level/fitting_LFCs.tsv": (
            "15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
            "15_insertion_level_curve_fitting/fitting_LFCs.tsv",
        ),
        "insertion_level/fitting_results.tsv": (
            "15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
            "15_insertion_level_curve_fitting/fitting_results.tsv",
        ),
        "insertion_level/transformed_weights.tsv": (
            "16_gene_level_depletion_analysis/gene_level_statistics.tsv",
            "16_gene_level_depletion_analysis/transformed_weights.tsv",
        ),
        "gene_level/LFC.tsv": (
            "16_gene_level_depletion_analysis/LFC.tsv",
            "16_gene_level_depletion_analysis/LFC.tsv",
        ),
        "gene_level/gene_level_fitting_statistics.tsv": (
            "17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
            "17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
        ),
        "gene_level/fitting_LFCs.tsv": (
            "17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
            "17_gene_level_curve_fitting/fitting_LFCs.tsv",
        ),
        "gene_level/fitting_results.tsv": (
            "17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
            "17_gene_level_curve_fitting/fitting_results.tsv",
        ),
    }

# annotations.tsv is produced uncompressed upstream but the app expects .gz.
_ANNOTATIONS_SRC = "12_concatenated/annotations.tsv"
_ANNOTATIONS_DST = "insertion_level/annotations.tsv.gz"

# imputation_statistics.tsv is a side-output of the DESeq2-replicates branch,
# anchored on that branch's declared imputed_raw_reads.tsv.
_HAS_IMPUTATION = config.get("use_DEseq2_for_biological_replicates", False)
_IMPUTATION_ANCHOR = "13_filtered/imputed_raw_reads.tsv"
_IMPUTATION_SRC = "13_filtered/imputation_statistics.tsv"
_IMPUTATION_DST = "insertion_level/imputation_statistics.tsv"


def _release_targets() -> list[str]:
    """All release files this project should produce (branch-aware)."""
    targets = [f"{_RELEASE}/{dst}" for dst in RELEASE_MAP]
    targets.append(f"{_RELEASE}/{_ANNOTATIONS_DST}")
    if _HAS_IMPUTATION:
        targets.append(f"{_RELEASE}/{_IMPUTATION_DST}")
    return targets


# Copy rule: input = declared anchor (drives the DAG), params.src = real file
# to copy (anchor or its sibling side-output). Constrained to mapped paths so
# the slash-containing {rel} wildcard stays unambiguous.
rule package_release_copy:
    input:
        lambda wildcards: f"{_RESULTS}/{RELEASE_MAP[wildcards.rel][0]}",
    output:
        f"{_RELEASE}/{{rel}}",
    params:
        src=lambda wildcards: f"{_RESULTS}/{RELEASE_MAP[wildcards.rel][1]}",
    wildcard_constraints:
        rel="|".join(re.escape(k) for k in RELEASE_MAP),
    message:
        "*** Packaging release file: {output}"
    shell:
        "mkdir -p $(dirname {output}) && cp {params.src} {output}"


# annotations.tsv -> annotations.tsv.gz (app expects gzip-compressed).
rule package_release_annotations:
    input:
        f"{_RESULTS}/{_ANNOTATIONS_SRC}",
    output:
        f"{_RELEASE}/{_ANNOTATIONS_DST}",
    message:
        "*** Packaging release file (gzip): {output}"
    shell:
        "mkdir -p $(dirname {output}) && gzip -c {input} > {output}"


# imputation_statistics.tsv (DESeq2-replicates branch only; side-output).
rule package_release_imputation:
    input:
        f"{_RESULTS}/{_IMPUTATION_ANCHOR}",
    output:
        f"{_RELEASE}/{_IMPUTATION_DST}",
    params:
        src=f"{_RESULTS}/{_IMPUTATION_SRC}",
    message:
        "*** Packaging release file: {output}"
    shell:
        "mkdir -p $(dirname {output}) && cp {params.src} {output}"


rule package_release:
    """Aggregate target: build the full release/ folder for this project."""
    input:
        _release_targets(),
    message:
        "*** Release packaged under projects/{project_name}/release"


# -----------------------------------------------------------------------------
# Reference data packaging (shared across projects, versioned separately).
# -----------------------------------------------------------------------------
_POMBASE_VERSION = config["Pombase_release_version"]
_POMBASE_DIR = f"resources/pombase_data/{_POMBASE_VERSION}"


rule package_reference:
    """Archive the versioned PomBase reference data for downstream use."""
    input:
        _POMBASE_DIR,
    output:
        f"resources/pombase_data/pombase_data_{_POMBASE_VERSION}.tar.gz",
    message:
        "*** Packaging reference data: {output}"
    shell:
        "tar -czf {output} -C resources/pombase_data {_POMBASE_VERSION}"
