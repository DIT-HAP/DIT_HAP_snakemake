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
        "results/{project_name}/13_filtered/control_insertions.tsv",
    log:
        "logs/{project_name}/depletion_scoring/control_insertion_selection.log",
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
            "results/{project_name}/13_filtered/imputed_raw_reads.tsv",
        log:
            "logs/{project_name}/depletion_scoring/impute_missing_values_using_FR.log",
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
                "results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level LFC",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
            padj=report(
                "results/{project_name}/14_insertion_level_depletion_analysis/padj.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level adjusted p-values",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
        log:
            "logs/{project_name}/depletion_scoring/insertion_level_depletion_analysis_has_replicates.log",
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
                "results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level LFC",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
        log:
            "logs/{project_name}/depletion_scoring/insertion_level_depletion_analysis_no_replicates.log",
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
        "results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
    output:
        report(
            "results/{project_name}/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
            category="Insertion-level results",
            labels={
                "name": "Insertion-level Curve Fitting Statistics",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "logs/{project_name}/depletion_scoring/insertion_level_curve_fitting.log",
    params:
        time_points=" ".join(map(str, config["time_points"])),
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
                "results/{project_name}/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv",
                category="Insertion-level results",
                labels={
                    "name": "Insertion-level LFC Fitted with R-square as Weights",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
        log:
            "logs/{project_name}/depletion_scoring/r_square_as_weights.log",
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Computing R-square as weights..."
        run:
            import pandas as pd

            fitting_res = pd.read_csv(input[0], sep="\t", index_col=[0, 1, 2, 3])
            fitting_res["confidence"] = fitting_res["R2"].clip(lower=1e-6, upper=1 - 1e-6)

            timepoint_columns = fitting_res.filter(regex=r".*_fitted$").columns.tolist()
            timepoint_columns = [col.rstrip("_fitted") for col in timepoint_columns]

            weights = pd.DataFrame(index=fitting_res.index, columns=timepoint_columns)
            for col in timepoint_columns:
                weights[col] = 1 - fitting_res["confidence"]

            weights.to_csv(output[0], sep="\t")


# Gene-level depletion analysis
# -----------------------------------------------------
rule gene_level_depletion_analysis:
    input:
        lfc_path="results/{project_name}/14_insertion_level_depletion_analysis/LFC.tsv",
        weights_path=branch(
            config.get("use_DEseq2_for_biological_replicates", False),
            "results/{project_name}/14_insertion_level_depletion_analysis/padj.tsv",
            "results/{project_name}/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv",
        ),
        annotations_path=rules.concat_counts_and_annotations.output.annotations,
    output:
        all_statistics="results/{project_name}/16_gene_level_depletion_analysis/gene_level_statistics.tsv",
        LFC=report(
            "results/{project_name}/16_gene_level_depletion_analysis/LFC.tsv",
            category="Gene-level results",
            labels={
                "name": "Gene-level LFC",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "logs/{project_name}/depletion_scoring/gene_level_depletion_analysis.log",
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


# Gene-level curve fitting
# -----------------------------------------------------
rule gene_level_curve_fitting:
    input:
        LFC=rules.gene_level_depletion_analysis.output.LFC,
    output:
        report(
            "results/{project_name}/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
            category="Gene-level results",
            labels={
                "name": "Gene-level Curve Fitting Statistics",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "logs/{project_name}/depletion_scoring/gene_level_curve_fitting.log",
    params:
        time_points=" ".join(map(str, config["time_points"])),
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
