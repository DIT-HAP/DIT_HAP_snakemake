# Plot MA plot (unified for both replicate and no-replicate branches)
# -----------------------------------------------------
rule plot_ma_plot:
    input:
        basemean="projects/{project_name}/results/14_insertion_level_depletion_analysis/baseMean.tsv",
        lfc="projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
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
        report(
            "projects/{project_name}/results/18_figure_data/ma_plot_horizontal.pdf",
            category="Insertion-level results",
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