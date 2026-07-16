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
        """
        echo "*** Running FastQC for PBL R1 reads..." > {log}
        fastqc --threads {threads} --noextract {input.PBL_r1} -o {params.output_dir} &>> {log}
        echo "*** Running FastQC for PBL R2 reads..." >> {log}
        fastqc --threads {threads} --noextract {input.PBL_r2} -o {params.output_dir} &>> {log}
        echo "*** Running FastQC for PBR R1 reads..." >> {log}
        fastqc --threads {threads} --noextract {input.PBR_r1} -o {params.output_dir} &>> {log}
        echo "*** Running FastQC for PBR R2 reads..." >> {log}
        fastqc --threads {threads} --noextract {input.PBR_r2} -o {params.output_dir} &>> {log}
        echo "*** FastQC analysis completed for all reads" >> {log}
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
        PBL_fastqc_r1=rules.fastqc_junction_classification.output.PBL_r1_html,
        PBL_fastqc_r2=rules.fastqc_junction_classification.output.PBL_r2_html,
        PBR_fastqc_r1=rules.fastqc_junction_classification.output.PBR_r1_html,
        PBR_fastqc_r2=rules.fastqc_junction_classification.output.PBR_r2_html,
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
        PBL=rules.bwa_mem_mapping.output.PBL,
        PBR=rules.bwa_mem_mapping.output.PBR,
    output:
        PBL_tsv=temp(f"projects/{project_name}/results/5_tabulated/{{sample}}_{{timepoint}}_{{condition}}.PBL.tsv"),
        PBR_tsv=temp(f"projects/{project_name}/results/5_tabulated/{{sample}}_{{timepoint}}_{{condition}}.PBR.tsv"),
    log:
        f"projects/{project_name}/logs/read_processing/bam_to_tsv/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/pysam.yml"
    # parse_bam_to_tsv.py is a single-core Python streaming loop; pysam's
    # decompression threads give ~no throughput gain (measured: threads 1/4/8
    # all ~99-102% CPU, ~175s, byte-identical output). Declaring 8 threads only
    # throttled scheduling to floor(cores/8) concurrent jobs. threads: 1 lets all
    # per-timepoint jobs run concurrently (floor(16/1)=16), which is the real win.
    threads: 1
    resources:
        # Flat-memory streaming parser; real RSS is a few hundred MB. 4 GB is
        # ample. (The previous 200000/200GB value only affects scheduling when
        # snakemake is run with an explicit --resources mem_mb budget, which this
        # project does not use, so it was never the concurrency limiter.)
        mem_mb=4000,
    message:
        "*** Transforming BAM to TSV for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Transforming PBL BAM to TSV..." > {log}
        python workflow/scripts/read_processing/parse_bam_to_tsv.py -i {input.PBL} -o {output.PBL_tsv} -t {threads} &>> {log}
        echo "*** Transforming PBR BAM to TSV..." >> {log}
        python workflow/scripts/read_processing/parse_bam_to_tsv.py -i {input.PBR} -o {output.PBR_tsv} -t {threads} &>> {log}
        """


# Filter aligned read pairs
# -----------------------------------------------------
rule filter_aligned_reads:
    input:
        PBL_tsv=rules.bam_to_tsv.output.PBL_tsv,
        PBR_tsv=rules.bam_to_tsv.output.PBR_tsv,
    output:
        PBL_filtered=f"projects/{project_name}/results/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.PBL.filtered.tsv",
        PBR_filtered=f"projects/{project_name}/results/6_filtered/{{sample}}_{{timepoint}}_{{condition}}.PBR.filtered.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/filter_aligned_reads/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        snakemake_config_file=config_file,
        chunk_size=config["chunk_size"],
    message:
        "*** Filtering aligned read pairs for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        python workflow/scripts/read_processing/filter_aligned_reads.py \
            -i {input.PBL_tsv} -o {output.PBL_filtered} \
            -c {params.chunk_size} --config {params.snakemake_config_file} &> {log}
        python workflow/scripts/read_processing/filter_aligned_reads.py \
            -i {input.PBR_tsv} -o {output.PBR_filtered} \
            -c {params.chunk_size} --config {params.snakemake_config_file} &>> {log}
        """


# Extract insertion sites
# -----------------------------------------------------
rule extract_insertion_sites:
    input:
        PBL_filtered=rules.filter_aligned_reads.output.PBL_filtered,
        PBR_filtered=rules.filter_aligned_reads.output.PBR_filtered,
    output:
        PBL_insertions=f"projects/{project_name}/results/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBL.tsv",
        PBR_insertions=f"projects/{project_name}/results/7_insertions/{{sample}}_{{timepoint}}_{{condition}}.PBR.tsv",
    log:
        f"projects/{project_name}/logs/read_processing/extract_insertion_sites/{{sample}}_{{timepoint}}_{{condition}}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    params:
        chunk_size=config["chunk_size"],
    message:
        "*** Extracting insertion sites for {wildcards.sample}_{wildcards.timepoint}_{wildcards.condition}..."
    shell:
        """
        echo "*** Extracting PBL insertion sites..." > {log}
        python workflow/scripts/read_processing/extract_insertion_sites.py \
            -i {input.PBL_filtered} -o {output.PBL_insertions} -c {params.chunk_size} &>> {log}
        echo "*** Extracting PBR insertion sites..." >> {log}
        python workflow/scripts/read_processing/extract_insertion_sites.py \
            -i {input.PBR_filtered} -o {output.PBR_insertions} -c {params.chunk_size} &>> {log}
        """


# Merge PBL + PBR strand insertions per sample/timepoint/condition
# -----------------------------------------------------
rule merge_strand_insertions:
    input:
        PBL_insertions=rules.extract_insertion_sites.output.PBL_insertions,
        PBR_insertions=rules.extract_insertion_sites.output.PBR_insertions,
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
