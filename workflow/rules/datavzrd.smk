# Datavzrd report for mapping filtering statistics
# -----------------------------------------------------
rule datavzrd_mapping_filtering_statistics:
    input:
        config="workflow/reports/datavzrd/mapping_filtering_statistics.yaml",
        table=rules.mapping_filtering_statistics.output[0],
    params:
        extra="",
    output:
        report(
            directory(f"projects/{project_name}/reports/mapping_filtering_statistics/datavzrd_mapping_filtering_statistics"),
            htmlindex="index.html",
            category="1. Quality Control",
            labels={
                "name": "2. Mapping Filtering Statistics",
                "type": "Datavzrd Report",
                "format": "Datavzrd HTML",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/mapping_filtering_statistics_datavzrd.log",
    wrapper:
        f"{snakemake_wrapper_version}/utils/datavzrd"

# Datavzrd report for insertion density statistics
# -----------------------------------------------------
rule datavzrd_insertion_density_analysis:
    input:
        config="workflow/reports/datavzrd/insertion_density_analysis.yaml",
        table=rules.insertion_density_data.output[0],
    params:
        extra="",
    output:
        report(
            directory(f"projects/{project_name}/reports/insertion_density_analysis/datavzrd_insertion_density_analysis"),
            htmlindex="index.html",
            category="1. Quality Control",
            labels={
                "name": "6b. Insertion Density (Table)",
                "type": "Datavzrd Report",
                "format": "Datavzrd HTML",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/insertion_density_analysis_datavzrd.log",
    wrapper:
        f"{snakemake_wrapper_version}/utils/datavzrd"

# Datavzrd report for gene-level curve fitting statistics
# -----------------------------------------------------
rule datavzrd_gene_level_curve_fitting:
    input:
        config="workflow/reports/datavzrd/gene_level_curve_fitting.yaml",
        table=rules.gene_level_curve_fitting.output[0],
    params:
        extra="",
    output:
        report(
            directory("projects/{project_name}/results/17_gene_level_curve_fitting/datavzrd_gene_level_curve_fitting"),
            htmlindex="index.html",
            category="3. Gene-level results",
            labels={
                "name": "Gene-level Curve Fitting Statistics (Table)",
                "type": "Datavzrd Report",
                "format": "Datavzrd HTML",
            },
        ),
    log:
        "projects/{project_name}/logs/depletion_scoring/gene_level_curve_fitting_datavzrd.log",
    wrapper:
        f"{snakemake_wrapper_version}/utils/datavzrd"
