# =============================================================================
# Snakefile — DIT-HAP pipeline entry point
# =============================================================================

from snakemake.utils import min_version, validate
from pathlib import Path
import pandas as pd

min_version("9.0")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Select the active project (edit this line to switch experiments).
# Each project lives under projects/{project}/ with its own config/ and outputs.
project = "LD_DIT_HAP"
# project = "HD_DIT_HAP"
# project = "HD_diploid"
# project = "LD_haploid"
# project = "Spikein"
# project = "Spore2YES6_1328"
config_file = f"projects/{project}/config/config.yaml"
configfile: config_file
validate(config, "workflow/schemas/config.schema.yaml")
workdir: "/data/c/yangyusheng_optimized/DIT_HAP_snakemake"

# ---------------------------------------------------------------------------
# Workflow lifecycle hooks
# ---------------------------------------------------------------------------
onstart:
    print("\n--- DIT-HAP analysis started ---\n")

onsuccess:
    print("\n--- Workflow finished successfully! ---\n")

onerror:
    print("\n--- An error occurred! ---\n")

# ---------------------------------------------------------------------------
# Project metadata
# ---------------------------------------------------------------------------
project_name = config["project_name"]
assert project_name == project, (
    f"Project directory '{project}' does not match config project_name "
    f"'{project_name}'. The projects/{{project}}/ folder name must equal "
    f"the project_name in its config.yaml."
)
snakemake_wrapper_version = config["snakemake_wrapper_version"]

# ---------------------------------------------------------------------------
# Sample sheet
# ---------------------------------------------------------------------------
sample_sheet = pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
validate(sample_sheet, "workflow/schemas/samples.schema.yaml")
samples    = sample_sheet["Sample"].unique().tolist()
timepoints = sample_sheet["Timepoint"].unique().tolist()
conditions = sample_sheet["Condition"].unique().tolist()

sample_sheet_dict = {
    s: {t: {c: {"fq1": None, "fq2": None} for c in conditions} for t in timepoints}
    for s in samples
}

for _, row in sample_sheet.iterrows():
    sample_sheet_dict[row["Sample"]][row["Timepoint"]][row["Condition"]]["fq1"] = Path(row["read1"])
    sample_sheet_dict[row["Sample"]][row["Timepoint"]][row["Condition"]]["fq2"] = Path(row["read2"])

# ---------------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------------
wildcard_constraints:
    sample    = "|".join(samples),
    timepoint = "|".join(timepoints),
    condition = "|".join(conditions),
    fragment  = "PBL|PBR",

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
# `all` lists exactly the report()-annotated deliverables, grouped by the
# categories they appear under in the Snakemake report. Intermediates (raw
# read tables, filtered counts, gene_level_statistics.tsv, …) are pulled in
# automatically as dependencies and are deliberately not listed here.
#
# Not part of `all` by design:
#   - reference data (resources/pombase_data/...): downloaded on demand by the
#     rules that consume it, only needed when bootstrapping a new project.
#   - release/ packaging: run separately via
#       snakemake --use-conda --cores <N> package_release
#   - the report itself:
#       snakemake --use-conda --cores <N> \
#           --report projects/{project_name}/reports/snakemake_report/report.zip \
#           --report-after-run -- all
rule all:
    input:
        # --- report category: Quality Control ---
        f"projects/{project_name}/reports/multiqc/quality_control_multiqc_report.html",
        f"projects/{project_name}/reports/mapping_filtering_statistics/datavzrd_mapping_filtering_statistics",
        f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
        f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis.pdf",
        f"projects/{project_name}/reports/insertion_orientation_analysis/insertion_orientation_analysis.pdf",
        f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis.pdf",
        f"projects/{project_name}/reports/insertion_density_analysis/datavzrd_insertion_density_analysis",
        f"projects/{project_name}/reports/gene_coverage_analysis/gene_coverage_analysis.pdf",

        # --- report category: Insertion-level results ---
        f"projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        f"projects/{project_name}/reports/ma_plot/ma_plot.pdf",
        f"projects/{project_name}/reports/ma_plot/ma_plot_horizontal.pdf",
        # DESeq2-replicates branch only (no p-values or dispersions without replicates):
        *(
            [
                f"projects/{project_name}/results/14_insertion_level_depletion_analysis/padj.tsv",
                f"projects/{project_name}/reports/dispersion_analysis/dispersion_analysis.pdf",
            ]
            if config.get("use_DEseq2_for_biological_replicates", False)
            else []
        ),
        # Curve fitting & gene level require a time course; QC-only projects
        # leave config["time_points"] unset and skip these.
        *(
            [
                f"projects/{project_name}/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
                f"projects/{project_name}/reports/insertion_level_curve_fitting/insertion_level_fitted_curves_sampled.pdf",
                f"projects/{project_name}/reports/distribution_of_curve_fitting/distribution_of_curve_fitting.pdf",
                f"projects/{project_name}/results/16_gene_level_depletion_analysis/insertion_weights.tsv",

                # --- report category: Gene-level results ---
                f"projects/{project_name}/results/16_gene_level_depletion_analysis/LFC.tsv",
                f"projects/{project_name}/results/16_gene_level_depletion_analysis/datavzrd_gene_level_LFC",
                f"projects/{project_name}/results/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
                f"projects/{project_name}/results/17_gene_level_curve_fitting/datavzrd_gene_level_curve_fitting",
                f"projects/{project_name}/reports/gene_level_curve_fitting/gene_level_fitted_curves_sampled.pdf",
            ]
            if config.get("time_points")
            else []
        ),


# ---------------------------------------------------------------------------
# Rule modules
# ---------------------------------------------------------------------------
include: "workflow/rules/reference_data.smk"
include: "workflow/rules/read_processing.smk"
include: "workflow/rules/depletion_scoring.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/packaging.smk"
include: "workflow/rules/figures.smk"
include: "workflow/rules/datavzrd.smk"
