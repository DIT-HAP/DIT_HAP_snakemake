#!/usr/bin/env bash

#===============================================================================
# PomBase Data Fetching Script
#===============================================================================
#
# OVERVIEW:
#   Downloads comprehensive datasets from PomBase (Schizosaccharomyces pombe
#   database) for specified release versions. This script automates the
#   acquisition of genome sequences, annotations, ontologies, and associated
#   metadata with automatic version detection and error recovery.
#
# DATA SOURCES:
#   - PomBase Monthly Releases: genome sequences, annotations, metadata
#   - Gene Ontology Consortium: GO terms and slim subsets
#   - GitHub API: FYPO and Mondo ontology version detection
#   - Monarch Initiative: Mondo disease ontology
#
# USAGE:
#   ./fetch_pombase_datasets.sh <release_version> [download_directory]
#
# PARAMETERS:
#   release_version    - PomBase release version in YYYY-MM-DD format
#                       (e.g., "2024-06-05")
#   download_directory - Target directory for downloaded files
#                       (optional, default: ./release/<version>)
#
# EXAMPLES:
#   # Download default release to default location
#   ./fetch_pombase_datasets.sh "2024-06-05"
#
#   # Download to specific directory
#   ./fetch_pombase_datasets.sh "2024-06-05" "/data/pombase"
#
# OUTPUT STRUCTURE:
#   <download_dir>/
#   ├── genome_sequence_and_features/    # FASTA sequences and GFF3 annotations
#   │   ├── Schizosaccharomyces_pombe_all_chromosomes.fa
#   │   ├── peptide.fa
#   │   └── Schizosaccharomyces_pombe_all_chromosomes.gff3
#   ├── Gene_metadata/                   # Gene information and viability data
#   │   ├── gene_IDs_names_products.tsv
#   │   └── gene_viability.tsv
#   ├── RNA_metadata/                    # Expression datasets
#   │   ├── qualitative_gene_expression.tsv
#   │   └── quantitative_gene_expression.tsv
#   ├── Protein_features/                # Protein annotations
#   │   ├── peptide_stats.tsv
#   │   ├── protein_families_and_domains.tsv
#   │   ├── disordered_regions.tsv
#   │   └── protein_modifications.tsv
#   └── ontologies_and_associations/     # Ontology files and annotations
#       ├── go-basic.obo                 # Gene Ontology core
#       ├── goslim_pombe.obo             # GO slim subset
#       ├── fypo-simple-pombase.obo      # Fission Yeast Phenotype Ontology
#       ├── mondo-simple.obo             # Mondo Disease Ontology
#       ├── gene_ontology_annotation.gaf.tsv    # GO annotations
#       ├── pombase_phenotype_annotation.phaf.tsv  # Phenotype annotations
#       └── human_disease_association.tsv        # Disease associations
#
# ERROR HANDLING:
#   - Validates release version format (YYYY-MM-DD)
#   - Verifies file downloads and non-zero file sizes
#   - Falls back to default ontology versions if API calls fail
#   - Creates directories with proper error checking
#   - Implements retry logic for failed downloads
#
# REQUIREMENTS:
#   - wget: for file downloads
#   - curl: for API calls to GitHub
#   - Internet connection to PomBase and GitHub
#
# VERSION: 2.0
# UPDATED: $(date +%Y-%m-%d)
#===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

#===============================================================================
# CONFIGURATION AND CONSTANTS
#===============================================================================

# Script identification
readonly SCRIPT_NAME="$(basename "${0}")"

# Default fallback ontology versions when GitHub API fails
readonly DEFAULT_FALLBACK_FYPO_VERSION="v2025-08-13"
readonly DEFAULT_FALLBACK_MONDO_VERSION="v2025-09-02"

# Base URLs for different data sources
readonly POMBASE_BASE_URL="https://www.pombase.org/monthly_releases"           # PomBase monthly releases
readonly GO_OBO_URL="https://purl.obolibrary.org/obo/go/go-basic.obo"           # Gene Ontology core
readonly GO_SLIM_URL="https://current.geneontology.org/ontology/subsets/goslim_pombe.obo"  # GO PomBe slim
readonly FYPO_GITHUB_API="https://api.github.com/repos/pombase/fypo/releases/latest"     # FYPO version API
readonly FYPO_RELEASE_URL="https://github.com/pombase/fypo/releases/download"            # FYPO download base
readonly MONDO_GITHUB_API="https://api.github.com/repos/monarch-initiative/mondo/releases/latest"  # Mondo version API
readonly MONDO_RELEASE_URL="https://github.com/monarch-initiative/mondo/releases/download"         # Mondo download base

# Wget configuration options
# -nc: no-clobber (don't overwrite existing files)
# --progress=bar:force: show progress bar always
# --timeout=30: timeout after 30 seconds
# --tries=3: retry up to 3 times
readonly WGET_OPTS="-nc --progress=bar:force --timeout=30 --tries=3"

#===============================================================================
# LOGGING AND OUTPUT FUNCTIONS
#===============================================================================

# Print formatted informational messages with timestamps
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*" >&2
}

# Print formatted warning messages with timestamps
log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*" >&2
}

# Print formatted error messages with timestamps
log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

#===============================================================================
# ARGUMENT PARSING AND VALIDATION
#===============================================================================

# Display usage information with detailed examples
show_usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} <release_version> [download_directory]

Download comprehensive PomBase datasets for a specific release version.

Arguments:
    release_version     PomBase release version (format: YYYY-MM-DD)
                       Example: "2024-06-05" for June 2024 release
    download_directory  Target directory (optional)
                       Default: ./release/<version>

Examples:
    # Download June 2024 release to default location
    ${SCRIPT_NAME} "2024-06-05"

    # Download to custom directory
    ${SCRIPT_NAME} "2024-06-05" "/data/pombase"

    # For help information
    ${SCRIPT_NAME} --help
    ${SCRIPT_NAME} -h

For detailed documentation about data sources and output structure,
see the script header documentation above.
EOF
}

# Validate release version format (YYYY-MM-DD)
validate_release_version() {
    local version="$1"

    # Use regex to validate YYYY-MM-DD format
    if [[ ! "${version}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
        log_error "Invalid release version format: '${version}'"
        log_error "Expected format: YYYY-MM-DD (e.g., 2025-08-01)"
        return 1
    fi

    # Additional validation for reasonable date ranges
    local year month day
    year=$(echo "${version}" | cut -d'-' -f1)
    month=$(echo "${version}" | cut -d'-' -f2)
    day=$(echo "${version}" | cut -d'-' -f3)

    # Basic validation (year >= 2000, month 1-12, day 1-31)
    if (( year < 2000 || year > 2030 )); then
        log_error "Invalid year in release version: ${year}"
        return 1
    fi

    if (( month < 1 || month > 12 )); then
        log_error "Invalid month in release version: ${month}"
        return 1
    fi

    if (( day < 1 || day > 31 )); then
        log_error "Invalid day in release version: ${day}"
        return 1
    fi

    return 0
}

# Create directory with comprehensive error handling
create_directory() {
    local dir_path="$1"
    local description="$2"

    log_info "Creating ${description} directory: ${dir_path}"

    # Check if directory already exists
    if [[ -d "${dir_path}" ]]; then
        log_warn "Directory already exists: ${dir_path}"
        return 0
    fi

    # Create directory with parent directories if needed
    if ! mkdir -p "${dir_path}"; then
        log_error "Failed to create ${description} directory: ${dir_path}"
        log_error "Check permissions and available disk space"
        return 1
    fi

    log_info "Successfully created directory: ${dir_path}"
    return 0
}

# Download file with comprehensive error handling and validation
download_file() {
    local url="$1"
    local target_dir="$2"
    local description="${3:-file}"

    log_info "Downloading ${description}: $(basename "${url}")"

    # Download with retry logic and progress tracking
    if ! wget ${WGET_OPTS} -P "${target_dir}" "${url}"; then
        log_error "Failed to download ${description}: ${url}"
        log_error "Network error, timeout, or file not available"
        return 1
    fi

    # Verify file was downloaded and has non-zero size
    local filename
    filename="$(basename "${url}")"
    local filepath="${target_dir}/${filename}"

    if [[ ! -f "${filepath}" ]]; then
        log_error "Downloaded file not found: ${filepath}"
        return 1
    fi

    if [[ ! -s "${filepath}" ]]; then
        log_error "Downloaded file is empty: ${filepath}"
        return 1
    fi

    # Log success with file size
    local file_size
    file_size=$(du -h "${filepath}" | cut -f1)
    log_info "Successfully downloaded: ${filename} (${file_size})"
    return 0
}

# Get latest ontology version from GitHub API with fallback
get_latest_ontology_version() {
    log_info "Fetching latest ontology version from GitHub API..."
    local ontology_github_api="$1"
    local default_fallback_ontology_version="$2"
    local ontology_version

    # Query GitHub API with timeout and error handling
    if ! ontology_version=$(curl -s --max-time 10 "${ontology_github_api}" 2>/dev/null); then
        log_warn "Failed to connect to GitHub API"
        ontology_version="${default_fallback_ontology_version}"
    else
        # Extract version from JSON response
        ontology_version=$(echo "${ontology_version}" | grep '"tag_name"' | cut -d'"' -f4 2>/dev/null)

        if [[ -z "${ontology_version}" ]]; then
            log_warn "Could not parse version from GitHub API response"
            ontology_version="${default_fallback_ontology_version}"
        else
            log_info "Latest ontology version detected: ${ontology_version}"
        fi
    fi

    if [[ "${ontology_version}" == "${default_fallback_ontology_version}" ]]; then
        log_warn "Using fallback ontology version: ${default_fallback_ontology_version}"
    fi

    echo "${ontology_version}"
}

#===============================================================================
# GENOME DATA DOWNLOAD FUNCTIONS
#===============================================================================

# Download genome sequence and annotation files
download_genome_files() {
    local base_url="$1"
    local target_dir="$2"

    log_info "Starting genome sequence and annotation downloads..."

    local genome_dir="${target_dir}/genome_sequence_and_features"
    create_directory "${genome_dir}" "genome sequence and features" || return 1

    # Download chromosome FASTA file (all chromosomes)
    log_info "Downloading chromosome sequences..."
    download_file \
        "${base_url}/genome_sequence_and_features/fasta_format/chromosomes/Schizosaccharomyces_pombe_all_chromosomes.fa" \
        "${genome_dir}" \
        "chromosome sequences" || return 1

    # Download peptide sequences FASTA file
    log_info "Downloading peptide sequences..."
    download_file \
        "${base_url}/genome_sequence_and_features/fasta_format/feature_sequences/peptide.fa" \
        "${genome_dir}" \
        "peptide sequences" || return 1

    # Download genome annotation file in GFF3 format
    log_info "Downloading genome annotation (GFF3)..."
    download_file \
        "${base_url}/genome_sequence_and_features/gff_format/Schizosaccharomyces_pombe_all_chromosomes.gff3" \
        "${genome_dir}" \
        "genome annotation (GFF3)" || return 1

    log_info "Genome files download completed successfully"
    return 0
}

#===============================================================================
# GENE METADATA DOWNLOAD FUNCTIONS
#===============================================================================

# Download gene metadata files
download_gene_metadata() {
    local base_url="$1"
    local target_dir="$2"

    log_info "Starting gene metadata downloads..."

    local metadata_dir="${target_dir}/Gene_metadata"
    create_directory "${metadata_dir}" "gene metadata" || return 1

    # Download gene identifiers, names, and product descriptions
    log_info "Downloading gene identifiers and names..."
    download_file \
        "${base_url}/gene_names_and_identifiers/gene_IDs_names_products.tsv" \
        "${metadata_dir}" \
        "gene IDs and names" || return 1

    # Download gene viability phenotype data
    log_info "Downloading gene viability data..."
    download_file \
        "${base_url}/phenotypes_and_genotypes/gene_viability.tsv" \
        "${metadata_dir}" \
        "gene viability data" || return 1

    log_info "Gene metadata download completed successfully"
    return 0
}

#===============================================================================
# RNA EXPRESSION DATA DOWNLOAD FUNCTIONS
#===============================================================================

# Download RNA expression data files
download_rna_metadata() {
    local base_url="$1"
    local target_dir="$2"

    log_info "Starting RNA metadata downloads..."

    local rna_dir="${target_dir}/RNA_metadata"
    create_directory "${rna_dir}" "RNA metadata" || return 1

    # Download qualitative gene expression data
    log_info "Downloading qualitative gene expression data..."
    download_file \
        "${base_url}/gene_expression/qualitative_gene_expression.tsv" \
        "${rna_dir}" \
        "qualitative gene expression" || return 1

    # Download quantitative gene expression data
    log_info "Downloading quantitative gene expression data..."
    download_file \
        "${base_url}/gene_expression/quantitative_gene_expression.tsv" \
        "${rna_dir}" \
        "quantitative gene expression" || return 1

    log_info "RNA metadata download completed successfully"
    return 0
}

#===============================================================================
# PROtein FEATURE DOWNLOAD FUNCTIONS
#===============================================================================

# Download protein feature and annotation files
download_protein_features() {
    local base_url="$1"
    local target_dir="$2"

    log_info "Starting protein features downloads..."

    local protein_dir="${target_dir}/Protein_features"
    create_directory "${protein_dir}" "protein features" || return 1

    # Download peptide length and composition statistics
    log_info "Downloading peptide statistics..."
    download_file \
        "${base_url}/protein_features/peptide_stats.tsv" \
        "${protein_dir}" \
        "peptide statistics" || return 1

    # Download protein family and domain annotations
    log_info "Downloading protein families and domains..."
    download_file \
        "${base_url}/protein_features/protein_families_and_domains.tsv" \
        "${protein_dir}" \
        "protein families and domains" || return 1

    # Download intrinsically disordered region annotations
    log_info "Downloading disordered regions..."
    download_file \
        "${base_url}/protein_features/disordered_regions.tsv" \
        "${protein_dir}" \
        "disordered regions" || return 1

    # Download post-translational modification annotations
    log_info "Downloading protein modifications..."
    download_file \
        "${base_url}/protein_features/protein_modifications.tsv" \
        "${protein_dir}" \
        "protein modifications" || return 1

    log_info "Protein features download completed successfully"
    return 0
}

#===============================================================================
# CURATED ORTHOLOGS DOWNLOAD FUNCTIONS
# ===============================================================================
download_curated_orthologs() {
    local base_url="$1"
    local target_dir="$2"

    log_info "Starting curated orthologs download..."

    local orthologs_dir="${target_dir}/curated_orthologs"
    create_directory "${orthologs_dir}" "curated orthologs" || return 1

    # Download curated pombe_japonicus_orthologs file
    log_info "Downloading curated pombe_japonicus_orthologs file..."
    download_file \
        "${base_url}/curated_orthologs/pombe_japonicus_orthologs.txt" \
        "${orthologs_dir}" \
        "curated orthologs" || return 1

    # Download curated pombe_cerevisiae_orthologs file
    log_info "Downloading curated pombe_cerevisiae_orthologs file..."
    download_file \
        "${base_url}/curated_orthologs/pombe_cerevisiae_orthologs.txt" \
        "${orthologs_dir}" \
        "curated orthologs" || return 1

    # Download curated pombe_human_orthologs file
    log_info "Downloading curated pombe_human_orthologs file..."
    download_file \
        "${base_url}/curated_orthologs/pombe_human_orthologs.txt" \
        "${orthologs_dir}" \
        "curated orthologs" || return 1

    log_info "Curated orthologs download completed successfully"
    return 0
}

#===============================================================================
# ONTOLOGY AND ANNOTATION DOWNLOAD FUNCTIONS
#===============================================================================

# Download ontology files and term associations
download_ontologies() {
    local base_url="$1"
    local target_dir="$2"
    local fypo_version="$3"
    local mondo_version="$4"

    log_info "Starting ontology and association downloads..."

    local onto_dir="${target_dir}/ontologies_and_associations"
    create_directory "${onto_dir}" "ontologies and associations" || return 1

    # Download core ontology files from external sources
    log_info "Downloading core ontology files..."
    download_file "${GO_OBO_URL}" "${onto_dir}" "GO basic ontology" || return 1
    download_file "${MONDO_RELEASE_URL}/${mondo_version}/mondo-simple.obo" "${onto_dir}" "Mondo ontology" || return 1
    download_file "${FYPO_RELEASE_URL}/${fypo_version}/fypo-simple-pombase.obo" "${onto_dir}" "FYPO ontology" || return 1
    download_file "${GO_SLIM_URL}" "${onto_dir}" "GO slim pombe" || return 1

    # Download slim ontology metadata files from PomBase
    log_info "Downloading slim ontology metadata..."
    local slim_files=(
        "gene_ontology/bp_go_slim_terms.tsv"                 # Biological Process slim
        "gene_ontology/cc_go_slim_terms.tsv"                 # Cellular Component slim
        "gene_ontology/mf_go_slim_terms.tsv"                 # Molecular Function slim
        "phenotypes_and_genotypes/fypo_slim_ids_and_names.tsv"  # FYPO slim terms
        "human_disease_annotation/pombe_mondo_disease_slim_terms.tsv"  # Disease slim terms
    )

    for file in "${slim_files[@]}"; do
        log_info "Downloading $(basename "${file}" .tsv) metadata..."
        download_file \
            "${base_url}/${file}" \
            "${onto_dir}" \
            "$(basename "${file}" .tsv) metadata" || return 1
    done

    # Download term association files
    log_info "Downloading term associations..."
    local association_files=(
        "gene_ontology/gene_ontology_annotation.gaf.tsv"              # GO annotations (GAF format)
        "macromolecular_complexes/macromolecular_complex_annotation.tsv"  # Complex annotations
        "phenotypes_and_genotypes/pombase_phenotype_annotation.phaf.tsv"  # Phenotype annotations (PHAF format)
        "phenotypes_and_genotypes/pombase_phenotype_annotation.eco.phaf.tsv"  # Phenotype annotations with ECO evidence
        "human_disease_annotation/human_disease_association.tsv"       # Human disease associations
    )

    for file in "${association_files[@]}"; do
        log_info "Downloading $(basename "${file}" .tsv) associations..."
        download_file \
            "${base_url}/${file}" \
            "${onto_dir}" \
            "$(basename "${file}" .tsv) associations" || return 1
    done

    log_info "Ontology and association downloads completed successfully"
    return 0
}

#===============================================================================
# MAIN EXECUTION
#===============================================================================

main() {
    # Parse and validate command-line arguments
    local release_version="$1"
    local download_dir="${2:-./release/${release_version}}"

    log_info "Starting PomBase data download for release: ${release_version}"
    log_info "Target directory: ${download_dir}"

    # Validate release version format
    validate_release_version "${release_version}" || {
        show_usage
        exit 1
    }

    # Extract year from release version for URL construction
    local year
    year=$(echo "${release_version}" | cut -d'-' -f1)
    local base_url="${POMBASE_BASE_URL}/${year}/pombase-${release_version}"
    log_info "Constructed PomBase URL: ${base_url}"

    # Create main download directory
    create_directory "${download_dir}" "main download" || exit 1

    # Get latest FYPO version with automatic fallback
    local fypo_version
    fypo_version=$(get_latest_ontology_version "${FYPO_GITHUB_API}" "${DEFAULT_FALLBACK_FYPO_VERSION}")

    # Get latest Mondo version with automatic fallback
    local mondo_version
    mondo_version=$(get_latest_ontology_version "${MONDO_GITHUB_API}" "${DEFAULT_FALLBACK_MONDO_VERSION}")

    # Execute download phases in sequence
    log_info "Beginning multi-phase download process..."

    # Phase 1: Download genome data
    download_genome_files "${base_url}" "${download_dir}" || exit 1

    # Phase 2: Download gene metadata
    download_gene_metadata "${base_url}" "${download_dir}" || exit 1

    # Phase 3: Download RNA expression data
    download_rna_metadata "${base_url}" "${download_dir}" || exit 1

    # Phase 4: Download protein features
    download_protein_features "${base_url}" "${download_dir}" || exit 1

    # Phase 5: Download curated orthologs
    download_curated_orthologs "${base_url}" "${download_dir}" || exit 1

    # Phase 6: Download ontologies and associations
    download_ontologies "${base_url}" "${download_dir}" "${fypo_version}" "${mondo_version}" || exit 1

    # Generate comprehensive completion summary
    log_info "All downloads completed successfully!"
    echo ""
    echo "==================== DOWNLOAD SUMMARY ===================="
    echo "PomBase Release: ${release_version}"
    echo "Download Location: ${download_dir}"
    echo "FYPO Version: ${fypo_version}"
    echo "Mondo Version: ${mondo_version}"
    echo ""
    echo "Directory Structure:"
    echo "  ${download_dir}/"
    echo "  ├── genome_sequence_and_features/    # FASTA and GFF3 files"
    echo "  ├── Gene_metadata/                   # Gene information and viability"
    echo "  ├── RNA_metadata/                    # Expression data (qualitative/quantitative)"
    echo "  ├── Protein_features/                # Protein annotations and modifications"
    echo "  ├── curated_orthologs/               # Curated ortholog datasets"
    echo "  └── ontologies_and_associations/     # Ontologies (GO, FYPO, Mondo) and associations"
    echo ""
    echo "Total downloaded files: $(find "${download_dir}" -type f | wc -l)"
    echo "=========================================================="
}

#===============================================================================
# SCRIPT EXECUTION
#===============================================================================

# Script entry point with comprehensive argument validation
if [[ $# -lt 1 || $# -gt 2 ]]; then
    log_error "Invalid number of arguments (expected 1-2, got $#)"
    show_usage
    exit 1
fi

# Handle help flags
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

# Validate required external tools before starting download
log_info "Checking required tools..."
for tool in wget curl; do
    if ! command -v "${tool}" &> /dev/null; then
        log_error "Required tool '${tool}' is not installed or not in PATH"
        log_error "Please install ${tool} and ensure it's available in your PATH"
        exit 1
    fi
    log_info "✓ ${tool} is available"
done

# Execute main function with all arguments
main "$@"
