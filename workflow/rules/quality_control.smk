# =============================================================================
# quality_control.smk — QC reports for the DIT-HAP pipeline
# =============================================================================


# FastQC on junction-classified reads
# -----------------------------------------------------
rule fastqc_junction_classification:
    input:
        PBL_r1=rules.junction_classification.output.PBL_r1,
        PBL_r2=rules.junction_classification.output.PBL_r2,
        PBR_r1=rules.junction_classification.output.PBR_r1,
        PBR_r2=rules.junction_classification.output.PBR_r2,
    output:
        PBL_r1_html=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_1_fastqc.html",
        PBL_r1_zip=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_1_fastqc.zip",
        PBL_r2_html=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_2_fastqc.html",
        PBL_r2_zip=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBL_2_fastqc.zip",
        PBR_r1_html=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_1_fastqc.html",
        PBR_r1_zip=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_1_fastqc.zip",
        PBR_r2_html=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_2_fastqc.html",
        PBR_r2_zip=f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.PBR_2_fastqc.zip",
    log:
        f"projects/{project_name}/logs/quality_control/fastqc/{{sample}}_{{timepoint}}_{{condition}}_fastqc_demultiplexed.log",
    conda:
        "../envs/fastqc.yml"
    params:
        output_dir=f"projects/{project_name}/reports/fastqc",
    threads: 4
    message:
        "*** Running FastQC for junction classified paired-end reads for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        # One fastqc invocation over all four files: --threads N means "process N
        # files concurrently", so a single call with 4 files uses all 4 threads,
        # whereas the previous four separate calls left threads idle (each ran one
        # file). Same outputs, ~4x faster.
        """
        fastqc --threads {threads} --noextract -o {params.output_dir} \
            {input.PBL_r1} {input.PBL_r2} {input.PBR_r1} {input.PBR_r2} &> {log}
        """


# Coordinate-sort and index BAMs (QC-only: feeds samtools/picard stats below,
# not the core bam_to_tsv path, which reads the name-sorted BAM directly)
# -----------------------------------------------------
rule samtools_sorting_and_indexing:
    input:
        PBL=rules.bwa_mem_mapping.output.PBL,
        PBR=rules.bwa_mem_mapping.output.PBR,
        ref_index=rules.samtools_faidx.output[0].format(release_version=config["Pombase_release_version"]),
    output:
        PBL_sorted=f"projects/{project_name}/results/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBL.sorted.bam",
        PBR_sorted=f"projects/{project_name}/results/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBR.sorted.bam",
        PBL_index=f"projects/{project_name}/results/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBL.sorted.bam.bai",
        PBR_index=f"projects/{project_name}/results/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.PBR.sorted.bam.bai",
    log:
        f"projects/{project_name}/logs/quality_control/samtools_sorting_and_indexing/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/bwa_mapping.yml"
    threads: 2
    message:
        "*** Sorting and indexing {input.PBL} and {input.PBR}..."
    shell:
        """
        echo "*** Sorting PBL reads..." > {log}
        samtools sort -@ {threads} {input.PBL} -O BAM -o {output.PBL_sorted} &>> {log}
        echo "*** Indexing PBL reads..." >> {log}
        samtools index -@ {threads} {output.PBL_sorted} &>> {log}
        echo "*** Sorting PBR reads..." >> {log}
        samtools sort -@ {threads} {input.PBR} -O BAM -o {output.PBR_sorted} &>> {log}
        echo "*** Indexing PBR reads..." >> {log}
        samtools index -@ {threads} {output.PBR_sorted} &>> {log}
        """


# Samtools mapping statistics
# -----------------------------------------------------
rule samtools_mapping_statistics:
    input:
        PBL=rules.samtools_sorting_and_indexing.output.PBL_sorted,
        PBR=rules.samtools_sorting_and_indexing.output.PBR_sorted,
    output:
        PBL_stats=f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBL.stats.txt",
        PBR_stats=f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBR.stats.txt",
        PBL_flagstat=f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBL.flagstat.txt",
        PBR_flagstat=f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBR.flagstat.txt",
        PBL_idxstats=f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBL.idxstats.txt",
        PBR_idxstats=f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.PBR.idxstats.txt",
    log:
        f"projects/{project_name}/logs/quality_control/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/bwa_mapping.yml"
    threads: 2
    message:
        "*** Running Samtools mapping statistics for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Samtools stats/flagstat/idxstats for PBL..." > {log}
        samtools stats    {input.PBL} > {output.PBL_stats}    2>> {log}
        samtools flagstat {input.PBL} > {output.PBL_flagstat} 2>> {log}
        samtools idxstats {input.PBL} > {output.PBL_idxstats} 2>> {log}
        echo "*** Samtools stats/flagstat/idxstats for PBR..." >> {log}
        samtools stats    {input.PBR} > {output.PBR_stats}    2>> {log}
        samtools flagstat {input.PBR} > {output.PBR_flagstat} 2>> {log}
        samtools idxstats {input.PBR} > {output.PBR_idxstats} 2>> {log}
        """


# Picard insert-size metrics
# -----------------------------------------------------
rule insert_size:
    input:
        f"projects/{project_name}/results/4_sorted/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.sorted.bam",
    output:
        txt=f"projects/{project_name}/reports/picard_insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.txt",
        pdf=f"projects/{project_name}/reports/picard_insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.pdf",
    log:
        f"projects/{project_name}/logs/quality_control/insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.log",
    params:
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE",
    resources:
        mem_mb=1024,
    wrapper:
        f"{snakemake_wrapper_version}/bio/picard/collectinsertsizemetrics"


# MultiQC for preprocessing reports
# -----------------------------------------------------
rule multiqc_preprocessing:
    input:
        fastp_json=expand(
            rules.fastp_preprocessing.output.json,
            sample=samples, timepoint=timepoints, condition=conditions
        ),
        junction_classification_json=expand(
            rules.junction_classification.output.json,
            sample=samples, timepoint=timepoints, condition=conditions
        ),
        fastqc=expand(
            f"projects/{project_name}/reports/fastqc/{{sample}}_{{timepoint}}_{{condition}}.{{read}}_fastqc.zip",
            sample=samples, timepoint=timepoints, condition=conditions,
            read=["PBL_1", "PBL_2", "PBR_1", "PBR_2"]
        ),
        bam_stats=expand(
            f"projects/{project_name}/reports/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.{{type}}.txt",
            sample=samples, timepoint=timepoints, condition=conditions,
            fragment=["PBL", "PBR"], type=["stats", "flagstat", "idxstats"]
        ),
        picard_metrics=expand(
            f"projects/{project_name}/reports/picard_insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.txt",
            sample=samples, timepoint=timepoints, condition=conditions,
            fragment=["PBL", "PBR"]
        ),
    output:
        report(
            f"projects/{project_name}/reports/multiqc/quality_control_multiqc_report.html",
            category="Quality Control",
            labels={
                "name": "1. Preprocessing MultiQC Report",
                "type": "MultiQC Report",
                "format": "HTML",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/multiqc_quality_control.log",
    conda:
        "../envs/multiqc.yml"
    params:
        outdir=f"projects/{project_name}/reports/multiqc",
        title=f"Quality Control Report for {project_name}",
        multiqc_config=config["multiqc_config"],
    message:
        "*** Generating MultiQC report for quality control..."
    shell:
        """
        multiqc \
            --title "{params.title}" \
            --filename quality_control_multiqc_report.html \
            --outdir {params.outdir} \
            --force \
            --verbose \
            --fn_as_s_name \
            --config {params.multiqc_config} \
            {input.fastp_json} \
            {input.junction_classification_json} \
            {input.fastqc} \
            {input.bam_stats} \
            {input.picard_metrics} \
            &> {log}
        """


# Merge per-fragment filter logs back into one combined log
# -----------------------------------------------------
# filter_aligned_reads was split into independent PBL/PBR jobs (each writes its
# own log). The QC parser (extract_mapping_filtering_statistics.py) keys the
# sample name off the log file stem and expects a single log containing both the
# PBL and PBR "FILTERING SUMMARY" blocks. This rule concatenates the two
# per-fragment logs into that combined form so the QC contract is unchanged.
rule merge_filter_logs:
    input:
        PBL=f"projects/{project_name}/logs/read_processing/filter_aligned_reads/{{sample}}_{{timepoint}}_{{condition}}.PBL.log",
        PBR=f"projects/{project_name}/logs/read_processing/filter_aligned_reads/{{sample}}_{{timepoint}}_{{condition}}.PBR.log",
    output:
        f"projects/{project_name}/logs/read_processing/filter_aligned_reads_combined/{{sample}}_{{timepoint}}_{{condition}}.log",
    message:
        "*** Merging PBL/PBR filter logs for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        cat {input.PBL} {input.PBR} > {output}
        """


# Mapping filtering statistics
# -----------------------------------------------------
rule mapping_filtering_statistics:
    input:
        expand(
            f"projects/{project_name}/logs/read_processing/filter_aligned_reads_combined/{{sample}}_{{timepoint}}_{{condition}}.log",
            sample=samples, timepoint=timepoints, condition=conditions,
        ),
    output:
        f"projects/{project_name}/reports/mapping_filtering_statistics/mapping_filtering_statistics.tsv",
    log:
        f"projects/{project_name}/logs/quality_control/mapping_filtering_statistics.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Generating mapping filtering statistics..."
    shell:
        """
        python workflow/scripts/quality_control/extract_mapping_filtering_statistics.py \
        -i {input} \
        -o {output} &> {log}
        """


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
            category="Quality Control",
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


# PBL-PBR correlation analysis
# -----------------------------------------------------
rule PBL_PBR_correlation_analysis:
    input:
        expand(rules.merge_strand_insertions.output, sample=samples, timepoint=timepoints, condition=conditions),
    output:
        report(
            f"projects/{project_name}/reports/PBL_PBR_correlation_analysis/PBL_PBR_correlation_analysis.pdf",
            caption="../reports/captions/PBL_PBR_correlation_analysis.rst",
            category="Quality Control",
            labels={
                "name": "3. PBL-PBR Correlation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/PBL_PBR_correlation_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing PBL-PBR correlation analysis..."
    shell:
        """
        python workflow/scripts/quality_control/PBL_PBR_correlation_analysis.py -i {input} -o {output} &> {log}
        """


# Read count distribution analysis
# -----------------------------------------------------
rule read_count_distribution_analysis:
    input:
        branch(
            config["merge_similar_timepoints"],
            expand(f"projects/{project_name}/results/11_merged/{{sample}}_{{condition}}.merged.tsv", sample=samples, condition=conditions),
            expand(rules.concat_timepoints.output.Reads, sample=samples, condition=conditions),
        ),
    output:
        report(
            f"projects/{project_name}/reports/read_count_distribution_analysis/read_count_distribution_analysis.pdf",
            caption="../reports/captions/read_count_distribution_analysis.rst",
            category="Quality Control",
            labels={
                "name": "4. Read Count Distribution Analysis",
                "type": "Distribution Plot",
                "format": "PDF",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/read_count_distribution_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing read count distribution analysis..."
    params:
        initial_time_point=config["initial_time_point"],
        hard_filtering_cutoff=config["hard_filtering_cutoff"],
    shell:
        """
        python workflow/scripts/quality_control/read_count_distribution_analysis.py \
            -i {input} \
            -t {params.initial_time_point} \
            -c {params.hard_filtering_cutoff} \
            -o {output} &> {log}
        """


# Insertion orientation analysis
# -----------------------------------------------------
rule insertion_orientation_analysis:
    input:
        rules.hard_filtering.output,
    output:
        report(
            f"projects/{project_name}/reports/insertion_orientation_analysis/insertion_orientation_analysis.pdf",
            caption="../reports/captions/insertion_orientation_analysis.rst",
            category="Quality Control",
            labels={
                "name": "5. Insertion Orientation Analysis",
                "type": "Correlation Plot",
                "format": "PDF",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/insertion_orientation_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing insertion orientation analysis..."
    shell:
        """
        python workflow/scripts/quality_control/insertion_orientation_analysis.py -i {input} -o {output} &> {log}
        """


# Insertion density analysis
# -----------------------------------------------------
rule insertion_density_analysis:
    input:
        insertion_data=rules.hard_filtering.output,
        annotation=rules.concat_counts_and_annotations.output.annotations,
    output:
        table=f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis.tsv",
        plot=report(
            f"projects/{project_name}/reports/insertion_density_analysis/insertion_density_analysis_histograms.pdf",
            caption="../reports/captions/insertion_density_analysis.rst",
            category="Quality Control",
            labels={
                "name": "6a. Insertion Density (Distributions)",
                "type": "Distribution Plot",
                "format": "PDF",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/insertion_density_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing insertion density analysis..."
    params:
        initial_time_point=config["initial_time_point"],
        final_time_point=config["final_time_point"],
    shell:
        """
        python workflow/scripts/quality_control/insertion_density_analysis.py \
            -i {input.insertion_data} \
            -a {input.annotation} \
            -t {params.initial_time_point} \
            -f {params.final_time_point} \
            -o {output.table} &> {log}
        """


# Datavzrd report for insertion density statistics
# -----------------------------------------------------
rule datavzrd_insertion_density_analysis:
    input:
        config="workflow/reports/datavzrd/insertion_density_analysis.yaml",
        table=rules.insertion_density_analysis.output.table,
    params:
        extra="",
    output:
        report(
            directory(f"projects/{project_name}/reports/insertion_density_analysis/datavzrd_insertion_density_analysis"),
            htmlindex="index.html",
            category="Quality Control",
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


# Gene coverage analysis
# -----------------------------------------------------
rule gene_coverage_analysis:
    input:
        insertion_data=f"projects/{project_name}/results/14_insertion_level_depletion_analysis/LFC.tsv",
        annotation=rules.concat_counts_and_annotations.output.annotations,
    output:
        report(
            f"projects/{project_name}/reports/gene_coverage_analysis/gene_coverage_analysis.pdf",
            caption="../reports/captions/gene_coverage_analysis.rst",
            category="Quality Control",
            labels={
                "name": "7. Gene Coverage",
                "type": "Coverage Plot",
                "format": "PDF",
            },
        ),
    log:
        f"projects/{project_name}/logs/quality_control/gene_coverage_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Performing gene coverage analysis..."
    params:
        gene_viability=rules.download_pombase_data.output.gene_viability.format(
            release_version=config["Pombase_release_version"]
        ),
    shell:
        """
        python workflow/scripts/quality_control/gene_coverage_analysis.py \
            -i {input.insertion_data} \
            -a {input.annotation} \
            -v {params.gene_viability} \
            -o {output} &> {log}
        """
