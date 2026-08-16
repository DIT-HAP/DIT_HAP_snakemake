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
        "../envs/statistics_and_computation.yml"
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
            "../envs/statistics_and_computation.yml"
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
            dispersion_data="projects/{project_name}/results/18_figure_data/dispersion_data.tsv",
            ma_values="projects/{project_name}/results/18_figure_data/ma_values_replicates.tsv",
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
            ma_values="projects/{project_name}/results/18_figure_data/ma_values.tsv",
        log:
            "projects/{project_name}/logs/depletion_scoring/insertion_level_depletion_analysis_no_replicates.log",
        params:
            initial_time_point=config["initial_time_point"],
        conda:
            "../envs/statistics_and_computation.yml"
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


# Plot MA plot
# -----------------------------------------------------
if not config.get("use_DEseq2_for_biological_replicates", False):

    rule plot_ma_plot:
        input:
            rules.insertion_level_depletion_analysis_no_replicates.output.ma_values,
        output:
            report(
                "projects/{project_name}/results/18_figure_data/ma_plot.pdf",
                category="Insertion-level results",
                labels={
                    "name": "MA Plot",
                    "type": "Figure",
                    "format": "PDF",
                },
            ),
        log:
            "projects/{project_name}/logs/figures/plot_ma_plot.log",
        params:
            stem=lambda wildcards, output: Path(output[0]).with_suffix(''),
        conda:
            "../envs/cnsplots.yml"
        message:
            "*** Plotting MA plot..."
        shell:
            """
            python workflow/scripts/figures/plot_ma_plot.py \
                -i {input} \
                -o {params.stem} &> {log}
            """


# Insertion-level curve fitting (computation only)
# -----------------------------------------------------
rule insertion_level_curve_fitting:
    input:
        "projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
    output:
        stats=report(
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
    threads: 16
    conda:
        "../envs/statistics_and_computation.yml"
    message:
        "*** Running insertion-level curve fitting..."
    shell:
        """
        python workflow/scripts/depletion_scoring/curve_fitting.py \
            -i {input} \
            -t {params.time_points} \
            -j {threads} \
            -o {output.stats} &> {log}
        """


# Plot insertion-level curve fitting
# -----------------------------------------------------
rule plot_insertion_level_curve_fitting:
    input:
        stats=rules.insertion_level_curve_fitting.output.stats,
        lfc="projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
    output:
        report(
            "projects/{project_name}/results/18_figure_data/insertion_level_fitted_curves.pdf",
            category="Insertion-level results",
            labels={
                "name": "Insertion-level Fitted Curves",
                "type": "Figure",
                "format": "PDF",
            },
        ),
    log:
        "projects/{project_name}/logs/figures/plot_insertion_level_curve_fitting.log",
    params:
        stem=lambda wildcards, output: Path(output[0]).with_suffix(''),
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Plotting insertion-level fitted curves..."
    shell:
        """
        python workflow/scripts/figures/plot_curve_fitting.py \
            -s {input.stats} \
            -l {input.lfc} \
            -o {params.stem} &> {log}
        """


# R-square as weights (no-replicate branch only)
# -----------------------------------------------------
if not config.get("use_DEseq2_for_biological_replicates", False):

    rule r_square_as_weights:
        input:
            rules.insertion_level_curve_fitting.output.stats,
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
            "../envs/statistics_and_computation.yml"
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
        "../envs/statistics_and_computation.yml"
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
    threads: 16
    conda:
        "../envs/statistics_and_computation.yml"
    message:
        "*** Running gene-level curve fitting..."
    shell:
        """
        python workflow/scripts/depletion_scoring/curve_fitting.py \
            -i {input.LFC} \
            -t {params.time_points} \
            -j {threads} \
            -o {output} &> {log}
        """


# Plot distribution of curve fitting results
# -----------------------------------------------------
rule plot_distribution_of_curve_fitting:
    input:
        rules.insertion_level_curve_fitting.output.stats,
    output:
        report(
            "projects/{project_name}/results/18_figure_data/distribution_of_curve_fitting.pdf",
            category="Insertion-level results",
            labels={
                "name": "Distribution of Curve Fitting Results",
                "type": "Figure",
                "format": "PDF",
            },
        ),
    log:
        "projects/{project_name}/logs/figures/plot_distribution_of_curve_fitting.log",
    params:
        stem=lambda wildcards, output: Path(output[0]).with_suffix(''),
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Plotting distribution of curve fitting results..."
    shell:
        """
        python workflow/scripts/figures/plot_distribution_of_curve_fitting.py \
            -i {input} \
            -o {params.stem} &> {log}
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
