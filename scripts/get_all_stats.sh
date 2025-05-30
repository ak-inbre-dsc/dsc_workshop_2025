#!/bin/bash
# Description: Script to collect SeqKit stats for various FASTQ outputs
# from the adaptive sampling analysis pipeline.

# --- Script Behavior ---
set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.
set -o pipefail # Cause a pipeline to return the exit status of the last command in the pipe that failed.

# --- 0. Configuration & User-Defined Variables ---
echo "--- Configuration ---"

if [ "$#" -eq 0 ]; then # Check if no arguments were provided
    echo "Usage: $0 <SAMPLEID>"
    echo "Error: SAMPLEID must be provided as the first command-line argument."
    exit 1
fi
SAMPLEID="$1"
echo "SAMPLEID set from command line: ${SAMPLEID}"

THREADS=4 # Adjust based on your system

# Main Output Directory Structure (mirroring the main processing script)
MAIN_OUTPUT_ROOT="${HOME}/Output_AS" # Root for all sample outputs
OUTPUT_DIR="${MAIN_OUTPUT_ROOT}/${SAMPLEID}"

# Derived Paths for input FASTQ files and output summaries
SUMMARY_DIR="${OUTPUT_DIR}/summary_data"
FASTQ_TRMT_DIR="${OUTPUT_DIR}/fastq_treatment"           # Location of treatment-specific FASTQs
FASTQ_TRMT_ISO_DIR="${OUTPUT_DIR}/fastq_treatment_iso"   # Location of mapped isolate-specific FASTQs

echo "SAMPLEID: ${SAMPLEID}"
echo "THREADS: ${THREADS}"
echo "MAIN_OUTPUT_ROOT: ${MAIN_OUTPUT_ROOT}"
echo "OUTPUT_DIR (Sample Specific): ${OUTPUT_DIR}"
echo "SUMMARY_DIR (Stats Output): ${SUMMARY_DIR}"
echo "FASTQ_TRMT_DIR (Treatment FASTQs): ${FASTQ_TRMT_DIR}"
echo "FASTQ_TRMT_ISO_DIR (Isolate FASTQs): ${FASTQ_TRMT_ISO_DIR}"
echo "--- End Configuration ---"
echo ""

# --- Create Summary Directory ---
echo "Ensuring summary directory exists: ${SUMMARY_DIR}"
mkdir -p "${SUMMARY_DIR}"
echo "Summary directory ensured."
echo ""

# --- Check for SeqKit ---
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit could not be found. Please install it first."
    echo "Installation instructions: https://bioinf.shenwei.me/seqkit/download/"
    exit 1
fi
echo "SeqKit found."
echo ""

# --- 1. Generate Stats for Treatment-Specific FASTQ Files ---
echo "--- Generating Stats for Treatment-Specific FASTQ Files ---"
# Files like: Output_AS/SAMPLEID/fastq_treatment/AS/SAMPLEID.reads.AS.fastq.gz
#             Output_AS/SAMPLEID/fastq_treatment/Control/SAMPLEID.reads.Control.fastq.gz

TRMT_FASTQ_INPUT_PATTERN="${FASTQ_TRMT_DIR}/*/${SAMPLEID}.reads.*.fastq.gz"
TRMT_STATS_OUTPUT_FILE="${SUMMARY_DIR}/${SAMPLEID}.reads.all_treatments.seqkit_stats.tsv"

echo "Searching for treatment FASTQ files matching: ${TRMT_FASTQ_INPUT_PATTERN}"
echo "Output will be saved to: ${TRMT_STATS_OUTPUT_FILE}"

# Check if files exist before running seqkit
if ! ls ${TRMT_FASTQ_INPUT_PATTERN} 1> /dev/null 2>&1; then
    echo "WARNING: No files found matching treatment FASTQ pattern: ${TRMT_FASTQ_INPUT_PATTERN}"
    echo "Skipping stats generation for treatment FASTQs."
else
    seqkit stats \
        -a \
        -b \
        -j "${THREADS}" \
        -T \
        ${TRMT_FASTQ_INPUT_PATTERN} \
        -o "${TRMT_STATS_OUTPUT_FILE}"

    if [ $? -eq 0 ] && [ -s "${TRMT_STATS_OUTPUT_FILE}" ]; then
        echo "Successfully generated stats for treatment FASTQs!"
        echo "Summary saved to: ${TRMT_STATS_OUTPUT_FILE}"
        echo "--- Treatment FASTQ Stats Summary ---"
        cat "${TRMT_STATS_OUTPUT_FILE}"
        echo "------------------------------------"
    else
        echo "Error generating stats for treatment FASTQs or output file is empty."
    fi
fi
echo ""


# --- 2. Generate Stats for Mapped Isolate-Specific FASTQ Files ---
echo "--- Generating Stats for Mapped Isolate-Specific FASTQ Files ---"
# Files like: Output_AS/SAMPLEID/fastq_treatment_iso/AS_Bacillus/SAMPLEID.mapped.AS.Bacillus.fastq.gz
#             Output_AS/SAMPLEID/fastq_treatment_iso/Control_Listeria/SAMPLEID.mapped.Control.Listeria.fastq.gz

ISO_FASTQ_INPUT_PATTERN="${FASTQ_TRMT_ISO_DIR}/*_*/${SAMPLEID}.mapped.*.*.fastq.gz"
ISO_STATS_OUTPUT_FILE="${SUMMARY_DIR}/${SAMPLEID}.mapped.all_isolates.seqkit_stats.tsv"

echo "Searching for mapped isolate FASTQ files matching: ${ISO_FASTQ_INPUT_PATTERN}"
echo "Output will be saved to: ${ISO_STATS_OUTPUT_FILE}"

# Check if files exist before running seqkit
if ! ls ${ISO_FASTQ_INPUT_PATTERN} 1> /dev/null 2>&1; then
    echo "WARNING: No files found matching mapped isolate FASTQ pattern: ${ISO_FASTQ_INPUT_PATTERN}"
    echo "Skipping stats generation for mapped isolate FASTQs."
else
    seqkit stats \
        -a \
        -b \
        -j "${THREADS}" \
        -T \
        ${ISO_FASTQ_INPUT_PATTERN} \
        -o "${ISO_STATS_OUTPUT_FILE}"

    if [ $? -eq 0 ] && [ -s "${ISO_STATS_OUTPUT_FILE}" ]; then
        echo "Successfully generated stats for mapped isolate FASTQs!"
        echo "Summary saved to: ${ISO_STATS_OUTPUT_FILE}"
        echo "--- Mapped Isolate FASTQ Stats Summary ---"
        cat "${ISO_STATS_OUTPUT_FILE}"
        echo "-----------------------------------------"
    else
        echo "Error generating stats for mapped isolate FASTQs or output file is empty."
    fi
fi
echo ""

echo "--- All SeqKit Stats Generation Finished ---"