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
snakemake_config_file = "config/config_HD_generationRAW.yaml"
# snakemake_config_file = "config/config_HD_generationPLUS1.yaml"
# snakemake_config_file = "config/config_LD_generationRAW.yaml"
# snakemake_config_file = "config/config_LD_generationPLUS1.yaml"
# snakemake_config_file = "config/config_HD_diploid.yaml"
# snakemake_config_file = "config/config_LD_haploid.yaml"
# snakemake_config_file = "config/config_spikein.yaml"
# snakemake_config_file = "config/config_1328_spore2YES6.yaml"
configfile: snakemake_config_file
validate(config, "workflow/schemas/config.schema.yaml")
workdir: "/data/c/yangyusheng_optimized/DIT_HAP"

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

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        # --- reference data ---
        # f"resources/pombase_data/{config['Pombase_release_version']}/genome_region/coding_gene_primary_transcripts.bed",
        # --- read processing ---
        # expand(f"results/{project_name}/10_annotated/{{sample}}_{{timepoint}}_{{condition}}.annotated.tsv", sample=samples, timepoint=timepoints, condition=conditions),
        # expand(f"results/{project_name}/11_concat_timepoints/{{sample}}_{{condition}}.counts.tsv", sample=samples, condition=conditions),
        # --- depletion scoring ---
        # f"results/{project_name}/13_filtered/raw_reads.filtered.tsv",
        # f"results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
        # f"results/{project_name}/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
        # f"results/{project_name}/16_gene_level_depletion_analysis/gene_level_statistics.tsv",
        # f"results/{project_name}/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
        # --- quality control ---
        # f"reports/{project_name}/multiqc/quality_control_multiqc_report.html",
        # f"reports/{project_name}/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
        # f"reports/{project_name}/insertion_density_analysis/insertion_density_analysis_histograms.pdf",
        # f"reports/{project_name}/gene_coverage_analysis",
        # --- smoke-test target (uncomment one to run) ---
        f"resources/pombase_data/{config['Pombase_release_version']}/genome_region/coding_gene_primary_transcripts.bed",

# ---------------------------------------------------------------------------
# Rule modules
# ---------------------------------------------------------------------------
include: "workflow/rules/reference_data.smk"
include: "workflow/rules/read_processing.smk"
include: "workflow/rules/depletion_scoring.smk"
include: "workflow/rules/quality_control.smk"
