#!/bin/bash
# Description: Script for processing and analyzing nanopore adaptive sampling data.
# It splits reads by treatment, maps them to a reference, and generates summaries.

# --- Script Behavior ---
# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u
# Cause a pipeline to return the exit status of the last command in the pipe that failed.
set -o pipefail

# --- 0. Configuration & User-Defined Variables ---
echo "--- Configuration ---"

# User-specified SAMPLEID - now taken from command-line argument
if [ "$#" -eq 0 ]; then # Check if no arguments were provided
    echo "Usage: $0 <SAMPLEID>"
    echo "Error: SAMPLEID must be provided as the first command-line argument."
    exit 1
fi
SAMPLEID="$1" # Now it's safe to access $1
echo "SAMPLEID set from command line: ${SAMPLEID}"


# Pipeline Variables
THREADS=8
MAXCHAN=256 # Channel threshold for "AS" treatment

# Input Data Location (e.g., from a Google Cloud Bucket or local path)
# Using ${HOME} to make it more portable if user's home directory is standard.
# Adjust BUCKET_DIR if your data is elsewhere.
BUCKET_DIR="${HOME}/dsc-nanopore-data/as_data"

# Reference Genome
REFERENCE_BASENAME="D6322.custom.reference.fasta" # Just the filename
REFERENCE="${BUCKET_DIR}/${REFERENCE_BASENAME}"   # Full path to the reference

# Main Output Directory (All outputs for this sample will go under here)
# This structure means SampleID is part of the OUTPUT_DIR path itself.
MAIN_OUTPUT_ROOT="${HOME}/Output_AS" # Root for all sample outputs
OUTPUT_DIR="${MAIN_OUTPUT_ROOT}/${SAMPLEID}"

# Define organism reference names for mapping (names of sequences within your REFERENCE fasta)
# These are used in Part 2 for mapping to specific organisms.
# Using direct variable assignments as in the original script.
Bacillus='Bacillus_subtilis_genome'
Enterococcus='Enterococcus_faecalis_genome'
Escherichia='Escherichia_coli_plasmid Escherichia_coli_chromosome'
Listeria='Listeria_monocytogenes_genome'
Pseudomonas='Pseudomonas_aeruginosa_genome'
Salmonella='Salmonella_enterica_genome'
Staphylococcus='Staphylococcus_aureus_chromosome Staphylococcus_aureus_plasmid1 Staphylococcus_aureus_plasmid2 Staphylococcus_aureus_plasmid3'

# Lists for looping
TRMT_LIST=("AS" "Control")
ISO_LIST=(
    "Bacillus"
    "Enterococcus"
    "Escherichia"
    "Listeria"
    "Pseudomonas"
    "Salmonella"
    "Staphylococcus"
)

# --- Derived Paths & Output Subdirectories ---
# These paths are relative to the OUTPUT_DIR.
FASTA_OUT_DIR="${OUTPUT_DIR}/fasta"
SUMMARY_DIR="${OUTPUT_DIR}/summary_data"
FASTQ_TRMT_DIR="${OUTPUT_DIR}/fastq_treatment"           # Stores filtered FASTQ per treatment
FASTQ_TRMT_MAPPED_DIR="${OUTPUT_DIR}/fastq_treatment_mapped" # Stores FASTQ of reads mapped to ALL community, per treatment
FASTQ_TRMT_ISO_DIR="${OUTPUT_DIR}/fastq_treatment_iso"       # Stores FASTQ of reads mapped to specific ISO, per treatment
MAPPED_DIR="${OUTPUT_DIR}/mapped"                           # Stores BAM files and mapping lists

# Path for the initial concatenated reads (output from cat, input to seqtk)
CONCAT_READS_PASS="${FASTA_OUT_DIR}/${SAMPLEID}.reads.pass.fastq.gz"

# Path to the original sequencing summary from the basecaller
BASECALLED_SEQ_SUMMARY="${BUCKET_DIR}/${SAMPLEID}/basecalled/sequencing_summary.txt"
BASECALLED_FASTQ_PASS_DIR="${BUCKET_DIR}/${SAMPLEID}/basecalled/fastq_pass"

echo "SAMPLEID: ${SAMPLEID}"
echo "THREADS: ${THREADS}"
echo "MAXCHAN: ${MAXCHAN}"
echo "BUCKET_DIR (Input Data): ${BUCKET_DIR}"
echo "REFERENCE: ${REFERENCE}"
echo "MAIN_OUTPUT_ROOT: ${MAIN_OUTPUT_ROOT}"
echo "OUTPUT_DIR (Sample Specific): ${OUTPUT_DIR}"
echo "--- End Configuration ---"
echo ""

# --- Create Base Output Directories ---
echo "Creating base output directories..."
mkdir -p "${OUTPUT_DIR}" # Main output directory for the sample
mkdir -p "${FASTA_OUT_DIR}"
mkdir -p "${SUMMARY_DIR}"
mkdir -p "${FASTQ_TRMT_DIR}"
mkdir -p "${FASTQ_TRMT_MAPPED_DIR}"
mkdir -p "${FASTQ_TRMT_ISO_DIR}"
mkdir -p "${MAPPED_DIR}"
echo "Base output directories ensured."
echo ""

# --- Verify Input Files ---
if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference file not found at ${REFERENCE}"
    exit 1
fi
if [ ! -f "${BASECALLED_SEQ_SUMMARY}" ]; then
    echo "ERROR: Basecalled sequencing summary not found at ${BASECALLED_SEQ_SUMMARY}"
    exit 1
fi
if [ ! -d "${BASECALLED_FASTQ_PASS_DIR}" ]; then
    echo "ERROR: Basecalled fastq_pass directory not found at ${BASECALLED_FASTQ_PASS_DIR}"
    exit 1
fi


# --- Part 1: Splitting the output by Treatment ---
echo "--- Part 1: Splitting output by Treatment ---"

# Collect all the pass reads into a single container
echo "Concatenating pass reads into: ${CONCAT_READS_PASS}"
if [ -z "$(ls -A ${BASECALLED_FASTQ_PASS_DIR}/*.fastq.gz 2>/dev/null)" ]; then
   echo "WARNING: No .fastq.gz files found in ${BASECALLED_FASTQ_PASS_DIR}/"
   # Decide if this is a fatal error or if the script can continue with an empty CONCAT_READS_PASS
   # For now, let's assume it's fatal if no reads are found.
   echo "ERROR: No input FASTQ files to process. Exiting."
   exit 1
else
   cat "${BASECALLED_FASTQ_PASS_DIR}/"*.fastq.gz > "${CONCAT_READS_PASS}"
   echo "Successfully concatenated reads to ${CONCAT_READS_PASS}"
fi
echo ""

# Process each treatment (AS, Control)
for TRMT in "${TRMT_LIST[@]}"; do
    echo "Processing Treatment: ${TRMT}"

    # Designate output files for this treatment
    TRMT_SEQSUM_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.reads.${TRMT}.seqsum.txt"
    TRMT_FASTQ_OUT="${FASTQ_TRMT_DIR}/${TRMT}/${SAMPLEID}.reads.${TRMT}.fastq.gz" # Output path for treatment-specific FASTQ
    TRMT_LST_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.reads.${TRMT}.lst"

    # Create specific subdirectory for this treatment's FASTQ
    mkdir -p "${FASTQ_TRMT_DIR}/${TRMT}"

    echo "  Output sequencing summary: ${TRMT_SEQSUM_OUT}"
    echo "  Output FASTQ (gzipped): ${TRMT_FASTQ_OUT}"
    echo "  Output read list: ${TRMT_LST_OUT}"

    # Generate Sequencing Summary for Treatment
    # Copy header
    echo "  Copying header to ${TRMT_SEQSUM_OUT}"
    head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${TRMT_SEQSUM_OUT}"

    # Define AWK conditional logic (this string will be part of the full awk script)
    # This is the part of the awk script that specifies the conditions for filtering.
    AWK_CONDITIONAL_LOGIC=""
    if [ "${TRMT}" == "AS" ]; then
        # MAXCHAN is expanded by the shell. $5, $10, $14 are for awk. "TRUE" is an awk string.
        AWK_CONDITIONAL_LOGIC='$5 <= '"${MAXCHAN}"' && $10 == "TRUE" && $14 >= 1000'
    elif [ "${TRMT}" == "Control" ]; then
        AWK_CONDITIONAL_LOGIC='$5 > '"${MAXCHAN}"' && $10 == "TRUE" && $14 >= 1000'
    else
        echo "ERROR: Unknown TRMT value: ${TRMT}"
        exit 1
    fi

    # Construct the full awk script for generating the summary (action: {print})
    AWK_SCRIPT_FOR_SUM="${AWK_CONDITIONAL_LOGIC} {print}"
    echo "  Filtering sequencing summary and appending to ${TRMT_SEQSUM_OUT}"
    awk -F'\t' -v OFS='\t' "${AWK_SCRIPT_FOR_SUM}" "${BASECALLED_SEQ_SUMMARY}" >> "${TRMT_SEQSUM_OUT}"

    # Construct the full awk script for generating the read list (action: {print $2})
    # The $2 needs to be escaped for the shell if the action string itself was in double quotes,
    # but here we are concatenating strings where $2 is meant for awk.
    AWK_SCRIPT_FOR_LST="${AWK_CONDITIONAL_LOGIC} {print \$2}" # Escaping $ for $2 to ensure awk sees it literally

    echo "  Generating read ID list: ${TRMT_LST_OUT}"
    awk -F'\t' "${AWK_SCRIPT_FOR_LST}" "${BASECALLED_SEQ_SUMMARY}" > "${TRMT_LST_OUT}"

    # Generate fastq of treatment using seqtk and gzip
    echo "  Generating FASTQ and gzipping: ${TRMT_FASTQ_OUT}"
    if [ -s "${CONCAT_READS_PASS}" ] && [ -s "${TRMT_LST_OUT}" ]; then
        seqtk subseq "${CONCAT_READS_PASS}" "${TRMT_LST_OUT}" | gzip -c > "${TRMT_FASTQ_OUT}"
        echo "  FASTQ for ${TRMT} generated and gzipped."
    else
        echo "  WARNING: Either ${CONCAT_READS_PASS} does not exist/is empty or ${TRMT_LST_OUT} is empty. Skipping FASTQ generation for ${TRMT}."
        # Create an empty gzipped file to prevent downstream errors if files are expected
        gzip -c < /dev/null > "${TRMT_FASTQ_OUT}"
        echo "  Created empty ${TRMT_FASTQ_OUT}."
    fi
    echo ""
done
echo "--- End of Part 1 ---"
echo ""


# --- Part 1b: Combine Treatment-Specific Sequencing Summaries ---
echo "--- Part 1b: Combining Treatment-Specific Sequencing Summaries ---"
COMBINED_READS_TRMT_SEQSUM_FILE="${SUMMARY_DIR}/${SAMPLEID}.reads.all_trmts.combined.seqsum.txt"
HEADER_WRITTEN_TRMT_COMBINE=false

echo "Output will be: ${COMBINED_READS_TRMT_SEQSUM_FILE}"

for TRMT in "${TRMT_LIST[@]}"; do
    INPUT_TRMT_SEQSUM_FILE="${FASTA_OUT_DIR}/${SAMPLEID}.reads.${TRMT}.seqsum.txt"
    echo "Processing for combination: ${INPUT_TRMT_SEQSUM_FILE}"

    if [ ! -f "${INPUT_TRMT_SEQSUM_FILE}" ]; then
        echo "  --> Warning: Input file for TRMT combine not found, skipping: ${INPUT_TRMT_SEQSUM_FILE}"
        continue
    fi

    if [ "${HEADER_WRITTEN_TRMT_COMBINE}" = false ]; then
        awk -v trmt_val="$TRMT" 'BEGIN{OFS="\t"} FNR==1 {print "TRMT", $0} FNR>1 {print trmt_val, $0}' \
            "${INPUT_TRMT_SEQSUM_FILE}" > "${COMBINED_READS_TRMT_SEQSUM_FILE}"
        HEADER_WRITTEN_TRMT_COMBINE=true
        echo "  --> Header written and first TRMT seqsum processed."
    else
        awk -v trmt_val="$TRMT" 'BEGIN{OFS="\t"} FNR>1 {print trmt_val, $0}' \
            "${INPUT_TRMT_SEQSUM_FILE}" >> "${COMBINED_READS_TRMT_SEQSUM_FILE}"
        echo "  --> Appended TRMT seqsum to combined file."
    fi
done
echo "Treatment-specific sequencing summaries combined."
echo ""


# --- Part 2: Analyzing the output (Mapping and Further Processing) ---
echo "--- Part 2: Analyzing the output ---"

# Loop through treatment conditions for mapping
for TRMT in "${TRMT_LIST[@]}"; do
    echo "--- Analyzing Treatment for Mapping: ${TRMT} ---"

    # Input FASTQ for this treatment (generated in Part 1)
    TRMT_FASTQ_INPUT="${FASTQ_TRMT_DIR}/${TRMT}/${SAMPLEID}.reads.${TRMT}.fastq.gz"

    # Check if the treatment FASTQ file exists and is not empty
    if [ ! -s "${TRMT_FASTQ_INPUT}" ]; then
        echo "  WARNING: Input FASTQ for mapping not found or empty: ${TRMT_FASTQ_INPUT}. Skipping mapping for ${TRMT}."
        continue
    fi

    # Output files for "ALL" community analysis for this treatment
    ANALYSIS_ALL_FASTQ_OUT="${FASTQ_TRMT_MAPPED_DIR}/${TRMT}/${SAMPLEID}.mapped.${TRMT}.ALL.fastq.gz"
    ANALYSIS_ALL_BAM_OUT="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.ALL.bam"
    ANALYSIS_ALL_LIST_OUT="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.ALL.lst"
    ANALYSIS_ALL_SEQSUM_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.mapped.${TRMT}.ALL.seqsum.txt" # Changed from OUTPUT_DIR/fasta to FASTA_OUT_DIR

    # Create specific subdirectory for this treatment's mapped FASTQ
    mkdir -p "${FASTQ_TRMT_MAPPED_DIR}/${TRMT}"

    echo "  Input FASTQ for mapping: ${TRMT_FASTQ_INPUT}"
    echo "  Output Mapped FASTQ (ALL, gzipped): ${ANALYSIS_ALL_FASTQ_OUT}"
    echo "  Output BAM (ALL mapped): ${ANALYSIS_ALL_BAM_OUT}"

    # Map reads to entire community reference
    # Using a more specific temporary file name for samtools sort
    SAMTOOLS_SORT_TEMP="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.reads.tmp"
    echo "  Mapping reads to entire community..."
    minimap2 -ax map-ont -t "${THREADS}" "${REFERENCE}" "${TRMT_FASTQ_INPUT}" | \
        samtools sort -@ "${THREADS}" -T "${SAMTOOLS_SORT_TEMP}" -o - | \
        samtools view -F 2308 -b -o "${ANALYSIS_ALL_BAM_OUT}" -
    echo "  Reads mapped to ${ANALYSIS_ALL_BAM_OUT}."

    # Index BAM and Extract list of all mapped reads
    echo "  Indexing BAM file: ${ANALYSIS_ALL_BAM_OUT}"
    samtools index "${ANALYSIS_ALL_BAM_OUT}"

    echo "  Extracting list of all mapped read IDs: ${ANALYSIS_ALL_LIST_OUT}"
    samtools view "${ANALYSIS_ALL_BAM_OUT}" | cut -f1 | sort -u > "${ANALYSIS_ALL_LIST_OUT}" # -u for unique

    # Generate FASTQ for all mapped reads
    echo "  Generating FASTQ of all mapped reads: ${ANALYSIS_ALL_FASTQ_OUT}"
    if [ -s "${ANALYSIS_ALL_LIST_OUT}" ]; then
        seqtk subseq "${TRMT_FASTQ_INPUT}" "${ANALYSIS_ALL_LIST_OUT}" | gzip -c > "${ANALYSIS_ALL_FASTQ_OUT}"
        echo "  FASTQ of all mapped reads generated."
    else
        echo "  WARNING: No reads mapped for ${TRMT} (list is empty). Skipping FASTQ generation for ALL mapped."
        gzip -c < /dev/null > "${ANALYSIS_ALL_FASTQ_OUT}"
        echo "  Created empty ${ANALYSIS_ALL_FASTQ_OUT}."
    fi

    # Generate sequencing summary for all mapped reads
    echo "  Generating sequencing summary for all mapped reads: ${ANALYSIS_ALL_SEQSUM_OUT}"
    head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${ANALYSIS_ALL_SEQSUM_OUT}"
    if [ -s "${ANALYSIS_ALL_LIST_OUT}" ]; then
        grep -F -f "${ANALYSIS_ALL_LIST_OUT}" "${BASECALLED_SEQ_SUMMARY}" >> "${ANALYSIS_ALL_SEQSUM_OUT}"
        echo "  Sequencing summary for all mapped reads generated."
    else
        echo "  WARNING: No reads mapped for ${TRMT}. Mapped sequencing summary will only contain header."
    fi
    echo ""

    # Pulls out mapped reads for each community member in each treatment
    echo "  --- Processing individual organisms for Treatment: ${TRMT} ---"
    for ISO in "${ISO_LIST[@]}"; do
        CURRENT_ISO_REF_NAME="${!ISO}" # Indirect expansion: if ISO="Bacillus", this becomes value of $Bacillus
        echo "    Processing Isolate: ${ISO} (Reference target(s): ${CURRENT_ISO_REF_NAME})"

        ANALYSIS_ISO_BAM_OUT="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.${ISO}.bam"
        ANALYSIS_ISO_LIST_OUT="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.${ISO}.lst"
        ANALYSIS_ISO_FASTQ_OUT="${FASTQ_TRMT_ISO_DIR}/${TRMT}_${ISO}/${SAMPLEID}.mapped.${TRMT}.${ISO}.fastq.gz"
        ANALYSIS_ISO_SEQSUM_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.mapped.${TRMT}.${ISO}.seqsum.txt" # Changed from OUTPUT_DIR/fasta

        # Create specific subdirectory for this treatment_iso FASTQ
        mkdir -p "${FASTQ_TRMT_ISO_DIR}/${TRMT}_${ISO}"

        echo "      Output BAM (${ISO}): ${ANALYSIS_ISO_BAM_OUT}"
        echo "      Output FASTQ (${ISO}, gzipped): ${ANALYSIS_ISO_FASTQ_OUT}"

        # Filter BAM for specific organism references
        # Ensure ANALYSIS_ALL_BAM_OUT exists and is not empty before proceeding
        if [ ! -s "${ANALYSIS_ALL_BAM_OUT}" ]; then
            echo "      WARNING: Parent BAM ${ANALYSIS_ALL_BAM_OUT} is empty or missing. Skipping isolate ${ISO}."
            # Create empty files to avoid downstream errors if files are expected
            samtools view -b -o "${ANALYSIS_ISO_BAM_OUT}" /dev/null # Creates an empty valid BAM
            touch "${ANALYSIS_ISO_LIST_OUT}"
            gzip -c < /dev/null > "${ANALYSIS_ISO_FASTQ_OUT}"
            head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${ANALYSIS_ISO_SEQSUM_OUT}"
            continue
        fi
        
        echo "      Filtering BAM for ${ISO}..."
        samtools view -@ "${THREADS}" -b "${ANALYSIS_ALL_BAM_OUT}" ${CURRENT_ISO_REF_NAME} -o "${ANALYSIS_ISO_BAM_OUT}"
        # Note: Space-separated CURRENT_ISO_REF_NAME is correctly handled by samtools view.
        echo "      BAM filtered for ${ISO}."

        echo "      Indexing BAM file for ${ISO}: ${ANALYSIS_ISO_BAM_OUT}"
        samtools index "${ANALYSIS_ISO_BAM_OUT}"

        # Extract read ID list for the organism (only mapped reads to this organism)
        echo "      Extracting read ID list for ${ISO}: ${ANALYSIS_ISO_LIST_OUT}"
        samtools view -F 0x04 "${ANALYSIS_ISO_BAM_OUT}" | cut -f1 | sort -u > "${ANALYSIS_ISO_LIST_OUT}"

        # Generate FASTQ for the organism using the ALL mapped FASTQ as input
        echo "      Generating FASTQ for ${ISO}: ${ANALYSIS_ISO_FASTQ_OUT}"
        if [ -s "${ANALYSIS_ALL_FASTQ_OUT}" ] && [ -s "${ANALYSIS_ISO_LIST_OUT}" ]; then
            seqtk subseq "${ANALYSIS_ALL_FASTQ_OUT}" "${ANALYSIS_ISO_LIST_OUT}" | gzip -c > "${ANALYSIS_ISO_FASTQ_OUT}"
            echo "      FASTQ for ${ISO} generated."
        else
            echo "      WARNING: Input FASTQ (${ANALYSIS_ALL_FASTQ_OUT}) or list (${ANALYSIS_ISO_LIST_OUT}) is empty/missing. Skipping FASTQ generation for ${ISO}."
            gzip -c < /dev/null > "${ANALYSIS_ISO_FASTQ_OUT}"
            echo "      Created empty ${ANALYSIS_ISO_FASTQ_OUT}."
        fi
        
        # Generate sequencing summary for this isolate's mapped reads
        echo "      Generating sequencing summary for ${ISO}: ${ANALYSIS_ISO_SEQSUM_OUT}"
        head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${ANALYSIS_ISO_SEQSUM_OUT}"
        if [ -s "${ANALYSIS_ISO_LIST_OUT}" ]; then
             grep -F -f "${ANALYSIS_ISO_LIST_OUT}" "${BASECALLED_SEQ_SUMMARY}" >> "${ANALYSIS_ISO_SEQSUM_OUT}"
             echo "      Sequencing summary for ${ISO} generated."
        else
            echo "      WARNING: Read list for ${ISO} is empty. Seqsum for ${ISO} will only contain header."
        fi
        echo ""
    done
    echo ""
done
echo "--- End of Part 2 ---"
echo ""


# --- Part 3: Combine Mapped Isolate-Specific Sequencing Summaries ---
echo "--- Part 3: Combining Mapped Isolate-Specific Sequencing Summaries ---"
COMBINED_MAPPED_ISO_SEQSUM_FILE="${SUMMARY_DIR}/${SAMPLEID}.mapped.all_trmts.all_isos.combined.seqsum.txt"
HEADER_WRITTEN_ISO_COMBINE=false

echo "Output will be: ${COMBINED_MAPPED_ISO_SEQSUM_FILE}"

for TRMT in "${TRMT_LIST[@]}"; do
    for ISO in "${ISO_LIST[@]}"; do
        INPUT_ISO_SEQSUM_FILE="${FASTA_OUT_DIR}/${SAMPLEID}.mapped.${TRMT}.${ISO}.seqsum.txt" # Path from where it was created
        echo "Processing for combination: ${INPUT_ISO_SEQSUM_FILE}"

        if [ ! -f "${INPUT_ISO_SEQSUM_FILE}" ]; then
            echo "  --> Warning: Input file for ISO combine not found, skipping: ${INPUT_ISO_SEQSUM_FILE}"
            continue
        fi

        if [ "${HEADER_WRITTEN_ISO_COMBINE}" = false ]; then
            awk -v trmt_val="$TRMT" -v iso_val="$ISO" 'BEGIN{OFS="\t"} FNR==1 {print "TRMT", "ISO", $0} FNR>1 {print trmt_val, iso_val, $0}' \
                "${INPUT_ISO_SEQSUM_FILE}" > "${COMBINED_MAPPED_ISO_SEQSUM_FILE}"
            HEADER_WRITTEN_ISO_COMBINE=true
            echo "  --> Header written and first ISO seqsum processed."
        else
            awk -v trmt_val="$TRMT" -v iso_val="$ISO" 'BEGIN{OFS="\t"} FNR>1 {print trmt_val, iso_val, $0}' \
                "${INPUT_ISO_SEQSUM_FILE}" >> "${COMBINED_MAPPED_ISO_SEQSUM_FILE}"
            echo "  --> Appended ISO seqsum to combined file."
        fi
    done
done
echo "Mapped isolate-specific sequencing summaries combined."
echo ""

echo "--- Adaptive Sampling Analysis Script Finished ---"
