# =============================================================================
# read_processing.smk — Read QC, mapping, insertion extraction, and annotation
# =============================================================================


# Fastp preprocessing for QC and adapter trimming
# -----------------------------------------------------
rule fastp_preprocessing:
    input:
        fq1=lambda wildcards: sample_sheet_dict[wildcards.sample][wildcards.timepoint][wildcards.condition]["fq1"],
        fq2=lambda wildcards: sample_sheet_dict[wildcards.sample][wildcards.timepoint][wildcards.condition]["fq2"],
    output:
        fq1=temp(f"projects/{project_name}/results/1_fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp_1.fq.gz"),
        fq2=temp(f"projects/{project_name}/results/1_fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp_2.fq.gz"),
        html=f"projects/{project_name}/reports/fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp.html",
        json=f"projects/{project_name}/reports/fastp/{{sample}}_{{timepoint}}_{{condition}}.fastp.json",
    log:
        f"projects/{project_name}/logs/read_processing/fastp/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/fastp.yml"
    params:
        adapter_sequence=config["adapter_sequence"],
        adapter_sequence_r2=config["adapter_sequence_r2"],
    threads: 6
    message:
        "*** Preprocessing fastp for {input.fq1} and {input.fq2}..."
    shell:
        """
        fastp --adapter_sequence {params.adapter_sequence} \
              --adapter_sequence_r2 {params.adapter_sequence_r2} \
              --disable_quality_filtering \
              --disable_length_filtering \
              --overrepresentation_analysis \
              --correction \
              -j {output.json} \
              -h {output.html} \
              --thread {threads} \
              --in1 {input.fq1} \
              --in2 {input.fq2} \
              --out1 {output.fq1} \
              --out2 {output.fq2} &> {log}
        """


# Cutadapt junction classification into PBL and PBR
# -----------------------------------------------------
rule junction_classification:
    input:
        fq1=rules.fastp_preprocessing.output.fq1,
        fq2=rules.fastp_preprocessing.output.fq2,
    output:
        PBL_r1=temp(f"projects/{project_name}/results/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBL_1.fq.gz"),
        PBL_r2=temp(f"projects/{project_name}/results/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBL_2.fq.gz"),
        PBR_r1=temp(f"projects/{project_name}/results/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBR_1.fq.gz"),
        PBR_r2=temp(f"projects/{project_name}/results/2_junction_classification/{{sample}}_{{timepoint}}_{{condition}}.PBR_2.fq.gz"),
        json=f"projects/{project_name}/reports/junction_classification/{{sample}}_{{timepoint}}_{{condition}}.json",
    log:
        f"projects/{project_name}/logs/read_processing/junction_classification/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/cutadapt.yml"
    params:
        PBL_adapter=config["PBL_adapter"],
        PBR_adapter=config["PBR_adapter"],
        PBL_reverseComplement_adapter=config["PBL_reverseComplement_adapter"],
        PBR_reverseComplement_adapter=config["PBR_reverseComplement_adapter"],
        output_folder=f"projects/{project_name}/results/2_junction_classification",
    threads: 6
    message:
        "*** Junction classification {input.fq1} and {input.fq2}..."
    shell:
        """
        cutadapt --cores {threads} \
                 -q 15 \
                 --overlap 15 \
                 -g PBL={params.PBL_adapter} \
                 -g PBR={params.PBR_adapter} \
                 -A PBL={params.PBL_reverseComplement_adapter} \
                 -A PBR={params.PBR_reverseComplement_adapter} \
                 -o {params.output_folder}/{wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}.{{name}}_1.fq.gz \
                 -p {params.output_folder}/{wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}.{{name}}_2.fq.gz \
                 --json {output.json} \
                 {input.fq1} {input.fq2} &> {log}
        """


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
        f"projects/{project_name}/logs/read_processing/fastqc/{{sample}}_{{timepoint}}_{{condition}}_fastqc_demultiplexed.log",
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


# BWA-MEM2 mapping, name-sorted BAM
# -----------------------------------------------------
rule bwa_mem_mapping:
    input:
        ref=rules.download_pombase_data.output.fasta.format(release_version=config["Pombase_release_version"]),
        ref_index=expand(rules.bwa_index.output, release_version=config["Pombase_release_version"]),
        PBL_fq1=rules.junction_classification.output.PBL_r1,
        PBL_fq2=rules.junction_classification.output.PBL_r2,
        PBR_fq1=rules.junction_classification.output.PBR_r1,
        PBR_fq2=rules.junction_classification.output.PBR_r2,
        # FastQC was previously listed here as an artificial input, forcing QC to
        # complete before mapping could start. Mapping does not need QC results, so
        # the dependency is removed to let FastQC and bwa run in parallel. FastQC
        # remains reachable via the `multiqc_preprocessing` rule (quality_control.smk),
        # which consumes the fastqc zips; request that target (or its report) to run QC.
    output:
        PBL=temp(f"projects/{project_name}/results/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.PBL.name_sorted.bam"),
        PBR=temp(f"projects/{project_name}/results/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.PBR.name_sorted.bam"),
    log:
        f"projects/{project_name}/logs/read_processing/bwa_mem_mapping/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/bwa_mapping.yml"
    threads: 8
    message:
        "*** Mapping {input.PBL_fq1} and {input.PBL_fq2} to {input.ref}..."
    shell:
        """
        echo "*** Mapping PBL reads..." > {log}
        bwa-mem2 mem -t {threads} {input.ref} {input.PBL_fq1} {input.PBL_fq2} 2>> {log} | \
            samtools sort -n -@ {threads} -O BAM -o {output.PBL} &>> {log}
        echo "*** Mapping PBR reads..." >> {log}
        bwa-mem2 mem -t {threads} {input.ref} {input.PBR_fq1} {input.PBR_fq2} 2>> {log} | \
            samtools sort -n -@ {threads} -O BAM -o {output.PBR} &>> {log}
        """


# Coordinate-sort and index BAMs
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
        f"projects/{project_name}/logs/read_processing/samtools_sorting_and_indexing/{{sample}}_{{timepoint}}_{{condition}}.log",
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
        f"projects/{project_name}/logs/read_processing/samtools_mapping_statistics/{{sample}}_{{timepoint}}_{{condition}}.log",
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
        f"projects/{project_name}/logs/read_processing/insert_size/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.log",
    params:
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE",
    resources:
        mem_mb=1024,
    wrapper:
        f"{snakemake_wrapper_version}/bio/picard/collectinsertsizemetrics"


# BAM → TSV (name-sorted BAM, no coordinate requirement)
# -----------------------------------------------------
rule bam_to_tsv:
    input:
        bam=f"projects/{project_name}/results/3_mapped/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.name_sorted.bam",
    output:
        tsv=temp(f"projects/{project_name}/results/5_tabulated/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.parquet"),
    log:
        f"projects/{project_name}/logs/read_processing/bam_to_tsv/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.log",
    conda:
        "../envs/pysam.yml"
    # parse_bam_to_tsv.py is a single-core Python streaming loop; pysam's
    # decompression threads give ~no throughput gain (measured: threads 1/4/8
    # all ~99-102% CPU, ~175s, byte-identical output). threads: 1 maximizes
    # scheduling concurrency (floor(cores/1)). PBL and PBR are now independent
    # jobs (fragment wildcard) so they run in parallel instead of serially.
    threads: 1
    resources:
        # Flat-memory streaming parser; real RSS is a few hundred MB. 4 GB is
        # ample. (mem_mb only constrains scheduling under an explicit
        # --resources mem_mb budget, which this project does not use.)
        mem_mb=4000,
    message:
        "*** Transforming {wildcards.fragment} BAM to Parquet for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/read_processing/parse_bam_to_tsv.py \
            -i {input.bam} -o {output.tsv} -t {threads} &> {log}
        """


# Filter aligned read pairs
# -----------------------------------------------------
rule filter_aligned_reads:
    input:
        tsv=rules.bam_to_tsv.output.tsv,
    output:
        filtered=f"projects/{project_name}/results/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.filtered.parquet",
    log:
        f"projects/{project_name}/logs/read_processing/filter_aligned_reads/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        snakemake_config_file=config_file,
        chunk_size=config["chunk_size"],
    message:
        "*** Filtering {wildcards.fragment} aligned read pairs for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/read_processing/filter_aligned_reads.py \
            -i {input.tsv} -o {output.filtered} \
            -c {params.chunk_size} --config {params.snakemake_config_file} &> {log}
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


# Extract insertion sites
# -----------------------------------------------------
rule extract_insertion_sites:
    input:
        filtered=rules.filter_aligned_reads.output.filtered,
    output:
        insertions=f"projects/{project_name}/results/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/extract_insertion_sites/{{sample}}_{{timepoint}}_{{condition}}.{{fragment}}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        chunk_size=config["chunk_size"],
    message:
        "*** Extracting {wildcards.fragment} insertion sites for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/read_processing/extract_insertion_sites.py \
            -i {input.filtered} -o {output.insertions} -c {params.chunk_size} &> {log}
        """


# Merge PBL + PBR strand insertions per sample/timepoint/condition
# -----------------------------------------------------
rule merge_strand_insertions:
    input:
        PBL_insertions=f"projects/{project_name}/results/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBL.tsv",
        PBR_insertions=f"projects/{project_name}/results/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBR.tsv",
    output:
        f"projects/{project_name}/results/8_merged/{{sample}}_{{timepoint}}_{{condition}}.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/merge_strand_insertions/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Merging strand insertions for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/read_processing/merge_strand_insertions.py \
            -i {input.PBL_insertions} \
            -j {input.PBR_insertions} \
            -o {output} &> {log}
        """


# Concatenate timepoints per sample/condition
# -----------------------------------------------------
rule concat_timepoints:
    input:
        counts=lambda wildcards: expand(
            rules.merge_strand_insertions.output,
            sample=wildcards.sample,
            timepoint=timepoints,
            condition=wildcards.condition,
        ),
        ref=rules.download_pombase_data.output.fasta.format(release_version=config["Pombase_release_version"]),
    output:
        PBL=f"projects/{project_name}/results/9_concatenated/{{sample}}_{{condition}}.PBL.tsv",
        PBR=f"projects/{project_name}/results/9_concatenated/{{sample}}_{{condition}}.PBR.tsv",
        Reads=f"projects/{project_name}/results/9_concatenated/{{sample}}_{{condition}}.Reads.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/concat_timepoints/{{sample}}_{{condition}}.log",
    conda:
        "../envs/biopython.yml"
    params:
        timepoints=" ".join(timepoints),
    shell:
        """
        python workflow/scripts/read_processing/concatenate_timepoint_data.py \
            -s {wildcards.sample}_{wildcards.condition} \
            -i {input.counts} \
            -tp {params.timepoints} \
            -g {input.ref} \
            -ol {output.PBL} \
            -or {output.PBR} \
            -o {output.Reads} &> {log}
        """


# Annotate insertions with genomic feature intervals
# -----------------------------------------------------
rule annotate_insertions:
    input:
        insertions=rules.concat_timepoints.output.Reads,
        genome_region=rules.extract_genome_region.output.genome_intervals_bed.format(
            release_version=config["Pombase_release_version"]
        ),
    output:
        f"projects/{project_name}/results/10_annotated/{{sample}}_{{condition}}.annotated.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/annotate_insertions/{{sample}}_{{condition}}.log",
    conda:
        "../envs/pybedtools.yml"
    message:
        "*** Annotating insertions for {wildcards.sample}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/read_processing/annotate_genomic_features.py \
            -i {input.insertions} -g {input.genome_region} -o {output} &> {log}
        """


# Optional: merge similar timepoints
# -----------------------------------------------------
if config["merge_similar_timepoints"]:
    rule merge_similar_timepoints:
        input:
            rules.concat_timepoints.output.Reads,
        output:
            f"projects/{project_name}/results/11_merged/{{sample}}_{{condition}}.merged.tsv",
        log:
            f"projects/{project_name}/logs/read_processing/merge_similar_timepoints/{{sample}}_{{condition}}.log",
        params:
            similar_timepoints=config["similar_timepoints"],
            merged_timepoint=config["merged_timepoint"],
            drop_columns=config["drop_columns"],
        conda:
            "../envs/statistics_and_figure_plotting.yml"
        message:
            "*** Merging similar time points for {wildcards.sample}_{wildcards.condition}..."
        shell:
            """
            python workflow/scripts/read_processing/merge_similar_timepoints.py \
                -i {input} -o {output} \
                -s {params.similar_timepoints} \
                -m {params.merged_timepoint} \
                -d {params.drop_columns} &> {log}
            """


# Concatenate all sample counts and annotations
# -----------------------------------------------------
rule concat_counts_and_annotations:
    input:
        counts=branch(
            config["merge_similar_timepoints"],
            expand(
                f"projects/{project_name}/results/11_merged/{{sample}}_{{condition}}.merged.tsv",
                sample=samples,
                condition=conditions,
            ),
            expand(rules.concat_timepoints.output.Reads, sample=samples, condition=conditions),
        ),
        annotations=expand(rules.annotate_insertions.output, sample=samples, condition=conditions),
    output:
        counts=f"projects/{project_name}/results/12_concatenated/raw_reads.tsv",
        annotations=f"projects/{project_name}/results/12_concatenated/annotations.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/concat_counts_and_annotations.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Concatenating counts and annotations..."
    shell:
        """
        python workflow/scripts/read_processing/concat_counts_and_annotations.py \
            -i {input.counts} \
            -a {input.annotations} \
            -oc {output.counts} \
            -oa {output.annotations} &> {log}
        """


# Hard-filter insertions by read count at initial timepoint
# -----------------------------------------------------
rule hard_filtering:
    input:
        rules.concat_counts_and_annotations.output.counts,
    output:
        f"projects/{project_name}/results/13_filtered/raw_reads.filtered.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/hard_filtering.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        cutoff=config["hard_filtering_cutoff"],
        init_timepoint=config["initial_time_point"],
    message:
        "*** Hard filtering insertions..."
    shell:
        """
        python workflow/scripts/read_processing/reads_hard_filtering.py \
            -i {input} -o {output} \
            -c {params.cutoff} -itp {params.init_timepoint} &> {log}
        """
