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


# =============================================================================
# Insertion-level analysis
# =============================================================================

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
                category="2. Insertion-level results",
                labels={
                    "name": "3. Insertion-level LFC",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
            padj="projects/{project_name}/results/14_insertion_level_depletion_analysis/padj.tsv",
            dispersion_data="projects/{project_name}/reports/dispersion_analysis/dispersion_data.tsv",
            basemean="projects/{project_name}/results/14_insertion_level_depletion_analysis/baseMean.tsv"
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
                -o {output.LFC} \
                --dispersion_data {output.dispersion_data} &> {log}
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
                category="2. Insertion-level results",
                labels={
                    "name": "3. Insertion-level LFC",
                    "type": "Statistics Table",
                    "format": "TSV",
                },
            ),
            basemean="projects/{project_name}/results/14_insertion_level_depletion_analysis/baseMean.tsv"
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

# Insertion-level curve fitting (computation only)
# -----------------------------------------------------
rule insertion_level_curve_fitting:
    input:
        "projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
    output:
        stats=report(
            "projects/{project_name}/results/15_insertion_level_curve_fitting/insertion_level_fitting_statistics.tsv",
            category="2. Insertion-level results",
            labels={
                "name": "4. Insertion-level Curve Fitting Statistics",
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


# Insertion-level aggregation weights
# -----------------------------------------------------
# The branch decision lives here, not in gene_level_depletion_analysis: with
# replicates, DESeq2's padj is the weighting signal; without replicates, there
# are no p-values at all, so rule 15's curve-fitting R2 is the only option. The
# `stats` input follows that same switch, so only one rule is needed.
rule compute_insertion_weights:
    input:
        stats=branch(
            config.get("use_DEseq2_for_biological_replicates", False),
            "projects/{project_name}/results/14_insertion_level_depletion_analysis/padj.tsv",
            rules.insertion_level_curve_fitting.output.stats,
        ),
        lfc="projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        annotations=rules.concat_counts_and_annotations.output.annotations,
    output:
        report(
            "projects/{project_name}/results/16_gene_level_depletion_analysis/insertion_weights.tsv",
            category="2. Insertion-level results",
            labels={
                "name": "5. Insertion-level Aggregation Weights",
                "type": "Statistics Table",
                "format": "TSV",
            },
        ),
    log:
        "projects/{project_name}/logs/depletion_scoring/compute_insertion_weights.log",
    params:
        has_replicates=(
            "--has-replicates"
            if config.get("use_DEseq2_for_biological_replicates", False)
            else "--no-has-replicates"
        ),
    conda:
        "../envs/statistics_and_computation.yml"
    message:
        "*** Computing insertion-level weights..."
    shell:
        """
        python workflow/scripts/depletion_scoring/compute_insertion_weights.py \
            {params.has_replicates} \
            -s {input.stats} -l {input.lfc} -a {input.annotations} \
            -o {output} &> {log}
        """


# =============================================================================
# Gene-level analysis
# =============================================================================

# Gene-level depletion analysis
# -----------------------------------------------------
rule gene_level_depletion_analysis:
    input:
        lfc_path="projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        weights_path=rules.compute_insertion_weights.output[0],
        annotations_path=rules.concat_counts_and_annotations.output.annotations,
    output:
        all_statistics="projects/{project_name}/results/16_gene_level_depletion_analysis/gene_level_statistics.tsv",
        LFC=report(
            "projects/{project_name}/results/16_gene_level_depletion_analysis/LFC.tsv",
            category="3. Gene-level results",
            labels={
                "name": "1. Gene-level LFC",
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

# Gene-level curve fitting
# -----------------------------------------------------
rule gene_level_curve_fitting:
    input:
        LFC=rules.gene_level_depletion_analysis.output.LFC,
    output:
        report(
            "projects/{project_name}/results/17_gene_level_curve_fitting/gene_level_fitting_statistics.tsv",
            category="3. Gene-level results",
            labels={
                "name": "3. Gene-level Curve Fitting Statistics",
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

