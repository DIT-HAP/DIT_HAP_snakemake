# =============================================================================
# reference_data.smk — Download PomBase reference data and extract genome regions
# =============================================================================


# download pombase data from ftp server
# -----------------------------------------------------
rule download_pombase_data:
    output:
        fasta = "resources/pombase_data/{release_version}/genome_sequence_and_features/Schizosaccharomyces_pombe_all_chromosomes.fa",
        gff = "resources/pombase_data/{release_version}/genome_sequence_and_features/Schizosaccharomyces_pombe_all_chromosomes.gff3",
        peptide_stats = "resources/pombase_data/{release_version}/Protein_features/peptide_stats.tsv",
        gene_IDs_names_products = "resources/pombase_data/{release_version}/Gene_metadata/gene_IDs_names_products.tsv",
        gene_viability = "resources/pombase_data/{release_version}/Gene_metadata/gene_viability.tsv",
    params:
        release_version = "{release_version}",
        download_dir = "resources/pombase_data/{release_version}",
    log:
        f"logs/{project_name}/reference_data/download_pombase_data_{{release_version}}.log",
    message:
        "*** Downloading PomBase data (release {wildcards.release_version}) from FTP server"
    shell:
        "workflow/scripts/reference_data/fetch_pombase_datasets.sh"
        " {params.release_version}"
        " {params.download_dir}"
        " &> {log}"


# index genome FASTA with samtools faidx
# -----------------------------------------------------
rule samtools_faidx:
    input:
        rules.download_pombase_data.output.fasta,
    output:
        "resources/pombase_data/{release_version}/genome_sequence_and_features/Schizosaccharomyces_pombe_all_chromosomes.fa.fai",
    log:
        f"logs/{project_name}/reference_data/samtools_faidx_{{release_version}}.log",
    message:
        "*** Indexing genome FASTA with samtools faidx"
    wrapper:
        "v5.8.3/bio/samtools/faidx"


# index genome FASTA with bwa-mem2
# -----------------------------------------------------
rule bwa_index:
    input:
        rules.download_pombase_data.output.fasta,
    output:
        rules.download_pombase_data.output.fasta + ".0123",
        rules.download_pombase_data.output.fasta + ".amb",
        rules.download_pombase_data.output.fasta + ".ann",
        rules.download_pombase_data.output.fasta + ".bwt.2bit.64",
        rules.download_pombase_data.output.fasta + ".pac",
    log:
        f"logs/{project_name}/reference_data/bwa_index_{{release_version}}.log",
    message:
        "*** Indexing genome FASTA with bwa-mem2"
    wrapper:
        "v5.8.3/bio/bwa-mem2/index"


# extract genome regions from GFF3 annotation
# -----------------------------------------------------
rule extract_genome_region:
    input:
        gff = rules.download_pombase_data.output.gff,
        fai = rules.samtools_faidx.output[0],
        peptide_stats = rules.download_pombase_data.output.peptide_stats,
        gene_ids = rules.download_pombase_data.output.gene_IDs_names_products,
        fypo = rules.download_pombase_data.output.gene_viability,
    output:
        primary_transcripts_bed = "resources/pombase_data/{release_version}/genome_region/coding_gene_primary_transcripts.bed",
        intergenic_regions_bed = "resources/pombase_data/{release_version}/genome_region/intergenic_regions.bed",
        non_coding_rna_bed = "resources/pombase_data/{release_version}/genome_region/non_coding_rna.bed",
        genome_intervals_bed = "resources/pombase_data/{release_version}/genome_region/genome_intervals.bed",
        overlapped_region_bed = "resources/pombase_data/{release_version}/genome_region/overlapped_region.bed",
    params:
        hayles = "resources/Literature/Hayles_2013_OB_merged_categories.xlsx",
    log:
        f"logs/{project_name}/reference_data/extract_genome_region_{{release_version}}.log",
    message:
        "*** Extracting genome regions from GFF3 annotation (release {wildcards.release_version})"
    conda:
        "../envs/pybedtools.yml"
    script:
        "../scripts/reference_data/extract_genome_region.py"
