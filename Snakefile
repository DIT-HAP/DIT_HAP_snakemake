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
rule all:
    input:
        # --- smoke-test target (uncomment one to run) ---
        f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf"
        # --- reference data ---
        # f"resources/pombase_data/{config['Pombase_release_version']}/genome_region/coding_gene_primary_transcripts.bed",

        # --- read processing ---
        # expand(f"projects/{project_name}/results/10_annotated/{{sample}}_{{timepoint}}_{{condition}}.annotated.tsv", sample=samples, timepoint=timepoints, condition=conditions),
        # expand(f"projects/{project_name}/results/11_concat_timepoints/{{sample}}_{{condition}}.counts.tsv", sample=samples, condition=conditions),
        
        # --- depletion scoring ---
        # f"projects/{project_name}/results/13_filtered/raw_reads.filtered.tsv",
        # f"projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        # f"projects/{project_name}/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
        # f"projects/{project_name}/results/16_gene_level_depletion_analysis/gene_level_statistics.tsv",
        # f"projects/{project_name}/results/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
        # --- quality control ---
        # f"projects/{project_name}/reports/multiqc/quality_control_multiqc_report.html",
        # f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
        # f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis_histograms.pdf",
        # f"projects/{project_name}/reports/gene_coverage_analysis",
        
        # --- packaging ---
        # see workflow/rules/packaging.smk's package_release for the release/
        # folder target; generate the Snakemake HTML/zip report separately with
        # `snakemake --use-conda --cores <N> \
        #     --report projects/{project_name}/reports/snakemake_report/report.zip \
        #     --report-after-run -- all`
        
        
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
