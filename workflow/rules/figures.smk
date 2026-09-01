# Plot MA plot (unified for both replicate and no-replicate branches)
# -----------------------------------------------------
rule plot_ma_plot:
    input:
        basemean="projects/{project_name}/results/14_insertion_level_depletion_analysis/baseMean.tsv",
        lfc="projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
    output:
        report(
            "projects/{project_name}/reports/ma_plot/ma_plot.pdf",
            category="2. Insertion-level results",
            labels={
                "name": "MA Plot",
                "type": "Figure",
                "format": "PDF",
            },
        ),
        report(
            "projects/{project_name}/reports/ma_plot/ma_plot_horizontal.pdf",
            category="2. Insertion-level results",
            labels={
                "name": "MA Plot (horizontal)",
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
            -b {input.basemean} \
            -l {input.lfc} \
            -o {params.stem} &> {log}
        """

# Only the DESeq2 branch estimates dispersions, so only it produces the figure data.
if config.get("use_DEseq2_for_biological_replicates", False):

    rule plot_dispersions:
        input:
            f"projects/{project_name}/reports/dispersion_analysis/dispersion_data.tsv",
        output:
            journal=report(
                f"projects/{project_name}/reports/dispersion_analysis/dispersion_analysis.pdf",
                caption="../reports/captions/dispersion_analysis.rst",
                category="2. Insertion-level results",
                labels={
                    "name": "DESeq2 Dispersion Estimates",
                    "type": "Figure",
                    "format": "PDF",
                },
            ),
            review=f"projects/{project_name}/reports/dispersion_analysis/dispersion_analysis.review.png",
        log:
            f"projects/{project_name}/logs/figures/plot_dispersions.log",
        params:
            stem=f"projects/{project_name}/reports/dispersion_analysis/dispersion_analysis",
        conda:
            "../envs/cnsplots.yml"
        message:
            "*** Rendering DESeq2 dispersion figure..."
        shell:
            """
            python workflow/scripts/figures/plot_dispersions.py \
                -i {input} \
                -o {params.stem} &> {log}
            """

rule plot_pbl_pbr_correlation:
    input:
        expand(rules.merge_strand_insertions.output, sample=samples, timepoint=timepoints, condition=conditions),
    output:
        journal=report(
            f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
            caption="../reports/captions/PBL_PBR_correlation_analysis.rst",
            category="1. Quality Control",
            labels={
                "name": "3. PBL-PBR Correlation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
            },
        ),
        review=f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.review.png",
    log:
        f"projects/{project_name}/logs/figures/plot_pbl_pbr_correlation.log",
    params:
        stem=f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Rendering PBL-PBR correlation figure..."
    shell:
        """
        python workflow/scripts/figures/plot_pbl_pbr_correlation.py \
            -i {input} \
            -o {params.stem} &> {log}
        """

rule plot_read_count_distribution:
    input:
        branch(
            config["merge_similar_timepoints"],
            expand(f"projects/{project_name}/results/11_merged/{{sample}}_{{condition}}.merged.tsv", sample=samples, condition=conditions),
            expand(rules.concat_timepoints.output.Reads, sample=samples, condition=conditions),
        ),
    output:
        journal=report(
            f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis.pdf",
            caption="../reports/captions/read_count_distribution_analysis.rst",
            category="1. Quality Control",
            labels={
                "name": "4. Read Count Distribution Analysis",
                "type": "Distribution Plot",
                "format": "PDF",
            },
        ),
        review=f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis.review.png",
        retention=f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis.retention.txt",
    log:
        f"projects/{project_name}/logs/figures/plot_read_count_distribution.log",
    params:
        stem=f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis",
        initial_time_point=config["initial_time_point"],
        hard_filtering_cutoff=config["hard_filtering_cutoff"],
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Rendering read count distribution figure..."
    shell:
        """
        python workflow/scripts/figures/plot_read_count_distribution.py \
            -i {input} \
            -q 0.999 \
            -t {params.initial_time_point} \
            -c {params.hard_filtering_cutoff} \
            -o {params.stem} &> {log}
        """

rule plot_insertion_orientation:
    input:
        rules.hard_filtering.output,
    output:
        journal=report(
            f"projects/{project_name}/reports/insertion_orientation_analysis/insertion_orientation_analysis.pdf",
            caption="../reports/captions/insertion_orientation_analysis.rst",
            category="1. Quality Control",
            labels={
                "name": "5. Insertion Orientation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
            },
        ),
        review=f"projects/{project_name}/reports/insertion_orientation_analysis/insertion_orientation_analysis.review.png",
    log:
        f"projects/{project_name}/logs/figures/plot_insertion_orientation.log",
    params:
        stem=f"projects/{project_name}/reports/insertion_orientation_analysis/insertion_orientation_analysis",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Rendering insertion orientation figure..."
    shell:
        """
        python workflow/scripts/figures/plot_insertion_orientation.py \
            -i {input} \
            -o {params.stem} &> {log}
        """

rule plot_insertion_density:
    input:
        rules.insertion_density_data.output,
    output:
        journal=report(
            f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis.pdf",
            caption="../reports/captions/insertion_density_analysis.rst",
            category="1. Quality Control",
            labels={
                "name": "6a. Insertion Density (Distributions)",
                "type": "Distribution Plot",
                "format": "PDF",
            },
        ),
        review=f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis.review.png",
    log:
        f"projects/{project_name}/logs/figures/plot_insertion_density.log",
    params:
        stem=f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis",
        initial_time_point=config["initial_time_point"],
        final_time_point=config["final_time_point"],
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Rendering insertion density figure..."
    shell:
        """
        python workflow/scripts/figures/plot_insertion_density.py \
            -i {input} \
            -t {params.initial_time_point} \
            -f {params.final_time_point} \
            -o {params.stem} &> {log}
        """

rule plot_gene_coverage:
    input:
        lfc=f"projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        annotation=rules.concat_counts_and_annotations.output.annotations,
        gene_viability=(
            f"resources/pombase_data/{config['Pombase_release_version']}/Gene_metadata/gene_viability.tsv"
        ),
    output:
        journal=report(
            f"projects/{project_name}/reports/gene_coverage_analysis/gene_coverage_analysis.pdf",
            caption="../reports/captions/gene_coverage_analysis.rst",
            category="1. Quality Control",
            labels={
                "name": "7. Gene Coverage",
                "type": "Coverage Plot",
                "format": "PDF",
            },
        ),
        review=f"projects/{project_name}/reports/gene_coverage_analysis/gene_coverage_analysis.review.png",
    log:
        f"projects/{project_name}/logs/figures/plot_gene_coverage.log",
    params:
        stem=f"projects/{project_name}/reports/gene_coverage_analysis/gene_coverage_analysis",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** Rendering gene coverage figure..."
    shell:
        """
        python workflow/scripts/figures/plot_gene_coverage.py \
            -i {input.lfc} \
            -a {input.annotation} \
            -v {input.gene_viability} \
            -o {params.stem} &> {log}
        """