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