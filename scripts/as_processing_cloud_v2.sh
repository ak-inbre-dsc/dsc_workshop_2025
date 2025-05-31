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

# Directory containing the original sequencing summary from the basecaller
BASECALLED_SUMMARY_DIR="${BUCKET_DIR}/${SAMPLEID}/basecalled"
BASECALLED_FASTQ_PASS_DIR="${BUCKET_DIR}/${SAMPLEID}/basecalled/fastq_pass"

# Dynamically find the sequencing summary file
# It should start with "sequencing_summary" and end with ".txt"
FOUND_SUMMARIES=()
if [ -d "${BASECALLED_SUMMARY_DIR}" ]; then
    # Use find to populate the array. -print0 and mapfile are safer for filenames with special characters,
    # but this approach is common and works if filenames don't contain newlines.
    mapfile -t FOUND_SUMMARIES < <(find "${BASECALLED_SUMMARY_DIR}" -maxdepth 1 -name "sequencing_summary*.txt" -print)
else
    echo "ERROR: Basecalled summary directory not found: ${BASECALLED_SUMMARY_DIR}"
    exit 1
fi

NUM_FOUND_SUMMARIES=${#FOUND_SUMMARIES[@]}
BASECALLED_SEQ_SUMMARY=""

if [ "${NUM_FOUND_SUMMARIES}" -eq 0 ]; then
    echo "ERROR: No sequencing summary file found in ${BASECALLED_SUMMARY_DIR} matching 'sequencing_summary*.txt'"
    exit 1
elif [ "${NUM_FOUND_SUMMARIES}" -gt 1 ]; then
    echo "ERROR: Multiple sequencing summary files found in ${BASECALLED_SUMMARY_DIR} matching 'sequencing_summary*.txt'. Please ensure only one exists:"
    printf "  %s\n" "${FOUND_SUMMARIES[@]}"
    exit 1
else
    BASECALLED_SEQ_SUMMARY="${FOUND_SUMMARIES[0]}"
    echo "Using sequencing summary file: ${BASECALLED_SEQ_SUMMARY}"
fi


echo "SAMPLEID: ${SAMPLEID}"
echo "THREADS: ${THREADS}"
echo "MAXCHAN: ${MAXCHAN}"
echo "BUCKET_DIR (Input Data): ${BUCKET_DIR}"
echo "REFERENCE: ${REFERENCE}"
echo "MAIN_OUTPUT_ROOT: ${MAIN_OUTPUT_ROOT}"
echo "OUTPUT_DIR (Sample Specific): ${OUTPUT_DIR}"
echo "--- End Configuration ---"
echo ""

# --- Dynamically Determine Column Numbers from Sequencing Summary Header ---
echo "Determining column numbers from header of ${BASECALLED_SEQ_SUMMARY}..."
HEADER_LINE=$(head -n1 "${BASECALLED_SEQ_SUMMARY}")
if [ -z "${HEADER_LINE}" ]; then
    echo "ERROR: Sequencing summary file is empty or header could not be read: ${BASECALLED_SEQ_SUMMARY}"
    exit 1
fi

IFS=$'\t' read -r -a HEADER_FIELDS <<< "$HEADER_LINE"

# Initialize column numbers to -1 (or an invalid value) to check if found
COL_CHANNEL_NUM=-1
COL_PASSES_FILTERING_NUM=-1
COL_SEQUENCE_LENGTH_TEMPLATE_NUM=-1
COL_READ_ID_NUM=-1 # For extracting read IDs for .lst files

for i in "${!HEADER_FIELDS[@]}"; do
    FIELD_NAME="${HEADER_FIELDS[$i]}"
    # Trim potential leading/trailing whitespace from header field (though unlikely with tab-separated)
    # FIELD_NAME=$(echo "$FIELD_NAME" | awk '{$1=$1;print}') # More robust way if spaces are an issue
    case "${FIELD_NAME}" in
        "channel") COL_CHANNEL_NUM=$((i+1)) ;;
        "passes_filtering") COL_PASSES_FILTERING_NUM=$((i+1)) ;;
        "sequence_length_template") COL_SEQUENCE_LENGTH_TEMPLATE_NUM=$((i+1)) ;;
        "read_id") COL_READ_ID_NUM=$((i+1)) ;; # Assuming the header name for read IDs is 'read_id'
    esac
done
unset IFS # Reset IFS to default

# Verify that all required columns were found
MISSING_HEADERS=0
if [ "$COL_CHANNEL_NUM" -eq -1 ]; then echo "ERROR: Header 'channel' not found in ${BASECALLED_SEQ_SUMMARY}"; MISSING_HEADERS=1; fi
if [ "$COL_PASSES_FILTERING_NUM" -eq -1 ]; then echo "ERROR: Header 'passes_filtering' not found in ${BASECALLED_SEQ_SUMMARY}"; MISSING_HEADERS=1; fi
if [ "$COL_SEQUENCE_LENGTH_TEMPLATE_NUM" -eq -1 ]; then echo "ERROR: Header 'sequence_length_template' not found in ${BASECALLED_SEQ_SUMMARY}"; MISSING_HEADERS=1; fi
if [ "$COL_READ_ID_NUM" -eq -1 ]; then echo "ERROR: Header 'read_id' not found in ${BASECALLED_SEQ_SUMMARY}"; MISSING_HEADERS=1; fi

if [ "$MISSING_HEADERS" -eq 1 ]; then
    echo "Please check the header names in your sequencing summary file."
    exit 1
fi
echo "Column mapping successful:"
echo "  Channel column: ${COL_CHANNEL_NUM}"
echo "  Passes_filtering column: ${COL_PASSES_FILTERING_NUM}"
echo "  Sequence_length_template column: ${COL_SEQUENCE_LENGTH_TEMPLATE_NUM}"
echo "  Read_id column: ${COL_READ_ID_NUM}"
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
    echo "ERROR: Basecalled sequencing summary was not successfully identified or is not a file: ${BASECALLED_SEQ_SUMMARY}"
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
    TRMT_FASTQ_OUT="${FASTQ_TRMT_DIR}/${TRMT}/${SAMPLEID}.reads.${TRMT}.fastq.gz" 
    TRMT_LST_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.reads.${TRMT}.lst"

    mkdir -p "${FASTQ_TRMT_DIR}/${TRMT}"

    echo "  Output sequencing summary: ${TRMT_SEQSUM_OUT}"
    echo "  Output FASTQ (gzipped): ${TRMT_FASTQ_OUT}"
    echo "  Output read list: ${TRMT_LST_OUT}"

    head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${TRMT_SEQSUM_OUT}"

    # Define AWK conditional logic using dynamically found column numbers
    AWK_CONDITIONAL_LOGIC=""
    # Note: Shell expands variables like ${COL_CHANNEL_NUM} and ${MAXCHAN} before awk sees the script.
    # Awk will see literal numbers and strings in their place.
    if [ "${TRMT}" == "AS" ]; then
        AWK_CONDITIONAL_LOGIC='$'${COL_CHANNEL_NUM}' <= '"${MAXCHAN}"' && $'${COL_PASSES_FILTERING_NUM}' == "TRUE" && $'${COL_SEQUENCE_LENGTH_TEMPLATE_NUM}' >= 1000'
    elif [ "${TRMT}" == "Control" ]; then
        AWK_CONDITIONAL_LOGIC='$'${COL_CHANNEL_NUM}' > '"${MAXCHAN}"' && $'${COL_PASSES_FILTERING_NUM}' == "TRUE" && $'${COL_SEQUENCE_LENGTH_TEMPLATE_NUM}' >= 1000'
    else
        echo "ERROR: Unknown TRMT value: ${TRMT}"
        exit 1
    fi

    AWK_SCRIPT_FOR_SUM="${AWK_CONDITIONAL_LOGIC} {print}"
    echo "  Filtering sequencing summary and appending to ${TRMT_SEQSUM_OUT}"
    awk -F'\t' -v OFS='\t' "${AWK_SCRIPT_FOR_SUM}" "${BASECALLED_SEQ_SUMMARY}" >> "${TRMT_SEQSUM_OUT}"

    # Construct the full awk script for generating the read list, using dynamic column for read_id
    AWK_SCRIPT_FOR_LST="${AWK_CONDITIONAL_LOGIC} {print \$${COL_READ_ID_NUM}}"
    echo "  Generating read ID list: ${TRMT_LST_OUT}"
    awk -F'\t' "${AWK_SCRIPT_FOR_LST}" "${BASECALLED_SEQ_SUMMARY}" > "${TRMT_LST_OUT}"

    echo "  Generating FASTQ and gzipping: ${TRMT_FASTQ_OUT}"
    if [ -s "${CONCAT_READS_PASS}" ] && [ -s "${TRMT_LST_OUT}" ]; then
        seqtk subseq "${CONCAT_READS_PASS}" "${TRMT_LST_OUT}" | gzip -c > "${TRMT_FASTQ_OUT}"
        echo "  FASTQ for ${TRMT} generated and gzipped."
    else
        echo "  WARNING: Either ${CONCAT_READS_PASS} does not exist/is empty or ${TRMT_LST_OUT} is empty. Skipping FASTQ generation for ${TRMT}."
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
    ANALYSIS_ALL_SEQSUM_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.mapped.${TRMT}.ALL.seqsum.txt" 

    mkdir -p "${FASTQ_TRMT_MAPPED_DIR}/${TRMT}"

    echo "  Input FASTQ for mapping: ${TRMT_FASTQ_INPUT}"
    echo "  Output Mapped FASTQ (ALL, gzipped): ${ANALYSIS_ALL_FASTQ_OUT}"
    echo "  Output BAM (ALL mapped): ${ANALYSIS_ALL_BAM_OUT}"

    SAMTOOLS_SORT_TEMP="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.reads.tmp"
    echo "  Mapping reads to entire community..."
    minimap2 -ax map-ont -t "${THREADS}" "${REFERENCE}" "${TRMT_FASTQ_INPUT}" | \
        samtools sort -@ "${THREADS}" -T "${SAMTOOLS_SORT_TEMP}" -o - | \
        samtools view -F 2308 -b -o "${ANALYSIS_ALL_BAM_OUT}" -
    echo "  Reads mapped to ${ANALYSIS_ALL_BAM_OUT}."

    echo "  Indexing BAM file: ${ANALYSIS_ALL_BAM_OUT}"
    samtools index "${ANALYSIS_ALL_BAM_OUT}"

    echo "  Extracting list of all mapped read IDs: ${ANALYSIS_ALL_LIST_OUT}"
    samtools view "${ANALYSIS_ALL_BAM_OUT}" | cut -f1 | sort -u > "${ANALYSIS_ALL_LIST_OUT}" 

    echo "  Generating FASTQ of all mapped reads: ${ANALYSIS_ALL_FASTQ_OUT}"
    if [ -s "${ANALYSIS_ALL_LIST_OUT}" ]; then
        seqtk subseq "${TRMT_FASTQ_INPUT}" "${ANALYSIS_ALL_LIST_OUT}" | gzip -c > "${ANALYSIS_ALL_FASTQ_OUT}"
        echo "  FASTQ of all mapped reads generated."
    else
        echo "  WARNING: No reads mapped for ${TRMT} (list is empty). Skipping FASTQ generation for ALL mapped."
        gzip -c < /dev/null > "${ANALYSIS_ALL_FASTQ_OUT}"
        echo "  Created empty ${ANALYSIS_ALL_FASTQ_OUT}."
    fi

    echo "  Generating sequencing summary for all mapped reads: ${ANALYSIS_ALL_SEQSUM_OUT}"
    head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${ANALYSIS_ALL_SEQSUM_OUT}"
    if [ -s "${ANALYSIS_ALL_LIST_OUT}" ]; then
        grep -F -f "${ANALYSIS_ALL_LIST_OUT}" "${BASECALLED_SEQ_SUMMARY}" >> "${ANALYSIS_ALL_SEQSUM_OUT}"
        echo "  Sequencing summary for all mapped reads generated."
    else
        echo "  WARNING: No reads mapped for ${TRMT}. Mapped sequencing summary will only contain header."
    fi
    echo ""

    echo "  --- Processing individual organisms for Treatment: ${TRMT} ---"
    for ISO in "${ISO_LIST[@]}"; do
        CURRENT_ISO_REF_NAME="${!ISO}" 
        echo "    Processing Isolate: ${ISO} (Reference target(s): ${CURRENT_ISO_REF_NAME})"

        ANALYSIS_ISO_BAM_OUT="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.${ISO}.bam"
        ANALYSIS_ISO_LIST_OUT="${MAPPED_DIR}/${SAMPLEID}.${TRMT}.${ISO}.lst"
        ANALYSIS_ISO_FASTQ_OUT="${FASTQ_TRMT_ISO_DIR}/${TRMT}_${ISO}/${SAMPLEID}.mapped.${TRMT}.${ISO}.fastq.gz"
        ANALYSIS_ISO_SEQSUM_OUT="${FASTA_OUT_DIR}/${SAMPLEID}.mapped.${TRMT}.${ISO}.seqsum.txt" 

        mkdir -p "${FASTQ_TRMT_ISO_DIR}/${TRMT}_${ISO}"

        echo "      Output BAM (${ISO}): ${ANALYSIS_ISO_BAM_OUT}"
        echo "      Output FASTQ (${ISO}, gzipped): ${ANALYSIS_ISO_FASTQ_OUT}"

        if [ ! -s "${ANALYSIS_ALL_BAM_OUT}" ]; then
            echo "      WARNING: Parent BAM ${ANALYSIS_ALL_BAM_OUT} is empty or missing. Skipping isolate ${ISO}."
            samtools view -b -o "${ANALYSIS_ISO_BAM_OUT}" /dev/null 
            touch "${ANALYSIS_ISO_LIST_OUT}"
            gzip -c < /dev/null > "${ANALYSIS_ISO_FASTQ_OUT}"
            head -n1 "${BASECALLED_SEQ_SUMMARY}" > "${ANALYSIS_ISO_SEQSUM_OUT}"
            continue
        fi
        
        echo "      Filtering BAM for ${ISO}..."
        samtools view -@ "${THREADS}" -b "${ANALYSIS_ALL_BAM_OUT}" ${CURRENT_ISO_REF_NAME} -o "${ANALYSIS_ISO_BAM_OUT}"
        echo "      BAM filtered for ${ISO}."

        echo "      Indexing BAM file for ${ISO}: ${ANALYSIS_ISO_BAM_OUT}"
        samtools index "${ANALYSIS_ISO_BAM_OUT}"

        echo "      Extracting read ID list for ${ISO}: ${ANALYSIS_ISO_LIST_OUT}"
        samtools view -F 0x04 "${ANALYSIS_ISO_BAM_OUT}" | cut -f1 | sort -u > "${ANALYSIS_ISO_LIST_OUT}"

        echo "      Generating FASTQ for ${ISO}: ${ANALYSIS_ISO_FASTQ_OUT}"
        if [ -s "${ANALYSIS_ALL_FASTQ_OUT}" ] && [ -s "${ANALYSIS_ISO_LIST_OUT}" ]; then
            seqtk subseq "${ANALYSIS_ALL_FASTQ_OUT}" "${ANALYSIS_ISO_LIST_OUT}" | gzip -c > "${ANALYSIS_ISO_FASTQ_OUT}"
            echo "      FASTQ for ${ISO} generated."
        else
            echo "      WARNING: Input FASTQ (${ANALYSIS_ALL_FASTQ_OUT}) or list (${ANALYSIS_ISO_LIST_OUT}) is empty/missing. Skipping FASTQ generation for ${ISO}."
            gzip -c < /dev/null > "${ANALYSIS_ISO_FASTQ_OUT}"
            echo "      Created empty ${ANALYSIS_ISO_FASTQ_OUT}."
        fi
        
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
        INPUT_ISO_SEQSUM_FILE="${FASTA_OUT_DIR}/${SAMPLEID}.mapped.${TRMT}.${ISO}.seqsum.txt" 
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
