# =============================================================================
# depletion_scoring.smk — Insertion filtering, depletion analysis, curve fitting
# =============================================================================


# Control insertion selection
# -----------------------------------------------------
rule control_insertion_selection:
    input:
        counts_df=rules.hard_filtering.output,
        annotation_df=rules.concat_counts_and_annotations.output.annotations,
    output:
        "projects/{project_name}/results/13_filtered/control_insertions.tsv",
    log:
        "projects/{project_name}/logs/depletion_scoring/control_insertion_selection.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Selecting control insertions..."
    shell:
        """
        python workflow/scripts/depletion_scoring/def_ctr_insertions.py \
            -i {input.counts_df} \
            -a {input.annotation_df} \
            -o {output} &> {log}
        """


# Distinguish replicates vs. no-replicates branches
# -----------------------------------------------------
if config.get("use_DEseq2_for_biological_replicates", False):

    # Imputation of missing values using Forward/Reverse insertions
    rule impute_missing_values_using_FR:
        input:
            filtered_reads=rules.hard_filtering.output,
            annotation=rules.concat_counts_and_annotations.output.annotations,
        output:
            "projects/{project_name}/results/13_filtered/imputed_raw_reads.tsv",
        log:
            "projects/{project_name}/logs/depletion_scoring/impute_missing_values_using_FR.log",
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Imputing missing values using FR..."
        shell:
            """
            python workflow/scripts/depletion_scoring/impute_missing_values_using_FR.py \
                -i {input.filtered_reads} \
                -a {input.annotation} \
                -o {output} &> {log}
            """

    # Insertion-level depletion analysis (DESeq2 — has replicates)
    rule insertion_level_depletion_analysis_has_replicates:
        input:
            counts_df=rules.impute_missing_values_using_FR.output,
            control_insertions_df=rules.control_insertion_selection.output,
        output:
            LFC=report(
                "projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level LFC",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
            padj=report(
                "projects/{project_name}/results/14_insertion_level_depletion_analysis/padj.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level adjusted p-values",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
        log:
            "projects/{project_name}/logs/depletion_scoring/insertion_level_depletion_analysis_has_replicates.log",
        params:
            initial_time_point=config["initial_time_point"],
        conda:
            "../envs/pydeseq2.yml"
        message:
            "*** Running insertion-level depletion analysis (DESeq2)..."
        shell:
            """
            python workflow/scripts/depletion_scoring/insertion_level_depletion_analysis_has_replicates.py \
                -i {input.counts_df} \
                -c {input.control_insertions_df} \
                -t {params.initial_time_point} \
                -o {output.LFC} &> {log}
            """

else:

    # Insertion-level depletion analysis (no replicates)
    rule insertion_level_depletion_analysis_no_replicates:
        input:
            counts_df=rules.hard_filtering.output,
            control_insertions_df=rules.control_insertion_selection.output,
        output:
            LFC=report(
                "projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level LFC",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
        log:
            "projects/{project_name}/logs/depletion_scoring/insertion_level_depletion_analysis_no_replicates.log",
        params:
            initial_time_point=config["initial_time_point"],
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Running insertion-level depletion analysis (no replicates)..."
        shell:
            """
            python workflow/scripts/depletion_scoring/insertion_level_depletion_analysis_no_replicates.py \
                -i {input.counts_df} \
                -c {input.control_insertions_df} \
                -t {params.initial_time_point} \
                -o {output.LFC} &> {log}
            """


# Insertion-level curve fitting
# -----------------------------------------------------
rule insertion_level_curve_fitting:
    input:
        "projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
    output:
        report(
            "projects/{project_name}/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
            category="Insertion-level results",
            labels={
                "name": "Insertion-level Curve Fitting Statistics",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "projects/{project_name}/logs/depletion_scoring/insertion_level_curve_fitting.log",
    params:
        time_points=lambda wildcards: " ".join(map(str, config["time_points"])),
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Running insertion-level curve fitting..."
    shell:
        """
        python workflow/scripts/depletion_scoring/curve_fitting.py \
            -i {input} \
            -t {params.time_points} \
            -o {output} &> {log}
        """


# R-square as weights (no-replicate branch only)
# -----------------------------------------------------
if not config.get("use_DEseq2_for_biological_replicates", False):

    rule r_square_as_weights:
        input:
            rules.insertion_level_curve_fitting.output,
        output:
            report(
                "projects/{project_name}/results/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level LFC Fitted with R-square as Weights",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
        log:
            "projects/{project_name}/logs/depletion_scoring/r_square_as_weights.log",
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Computing R-square as weights..."
        shell:
            """
            python workflow/scripts/depletion_scoring/compute_r2_weights.py \
                -i {input} -o {output} &> {log}
            """


# Gene-level depletion analysis
# -----------------------------------------------------
rule gene_level_depletion_analysis:
    input:
        lfc_path="projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        weights_path=branch(
            config.get("use_DEseq2_for_biological_replicates", False),
            "projects/{project_name}/results/14_insertion_level_depletion_analysis/padj.tsv",
            "projects/{project_name}/results/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv",
        ),
        annotations_path=rules.concat_counts_and_annotations.output.annotations,
    output:
        all_statistics="projects/{project_name}/results/16_gene_level_depletion_analysis/gene_level_statistics.tsv",
        LFC=report(
            "projects/{project_name}/results/16_gene_level_depletion_analysis/LFC.tsv",
            category="Gene-level results",
            labels={
                "name": "Gene-level LFC",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "projects/{project_name}/logs/depletion_scoring/gene_level_depletion_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Running gene-level depletion analysis..."
    shell:
        """
        python workflow/scripts/depletion_scoring/gene_level_depletion_analysis.py \
            -l {input.lfc_path} \
            -a {input.annotations_path} \
            -w {input.weights_path} \
            -o {output.all_statistics} &> {log}
        """


# Datavzrd report for gene-level LFC
# -----------------------------------------------------
rule datavzrd_gene_level_LFC:
    input:
        config="workflow/reports/datavzrd/gene_level_LFC.yaml",
        table=rules.gene_level_depletion_analysis.output.LFC,
    params:
        extra="",
    output:
        report(
            directory("projects/{project_name}/results/16_gene_level_depletion_analysis/datavzrd_gene_level_LFC"),
            htmlindex="index.html",
            category="Gene-level results",
            labels={
                "name": "Gene-level LFC (Table)",
                "type": "Datavzrd Report",
                "format": "Datavzrd HTML",
            },
        ),
    log:
        "projects/{project_name}/logs/depletion_scoring/gene_level_LFC_datavzrd.log",
    wrapper:
        f"{snakemake_wrapper_version}/utils/datavzrd"


# Gene-level curve fitting
# -----------------------------------------------------
rule gene_level_curve_fitting:
    input:
        LFC=rules.gene_level_depletion_analysis.output.LFC,
    output:
        report(
            "projects/{project_name}/results/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
            category="Gene-level results",
            labels={
                "name": "Gene-level Curve Fitting Statistics",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "projects/{project_name}/logs/depletion_scoring/gene_level_curve_fitting.log",
    params:
        time_points=lambda wildcards: " ".join(map(str, config["time_points"])),
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Running gene-level curve fitting..."
    shell:
        """
        python workflow/scripts/depletion_scoring/curve_fitting.py \
            -i {input.LFC} \
            -t {params.time_points} \
            -o {output} &> {log}
        """


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
            category="Gene-level results",
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
