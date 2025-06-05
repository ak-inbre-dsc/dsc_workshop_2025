# R Script for Generating Violin Plots from Combined Sequencing Summary Data

# Install packages if you haven't already
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("ggsci") # For scientific journal color palettes
# install.packages("patchwork") # For combining plots (optional, if needed for Part 3)
# install.packages("here") # Optional, for easier path management

library(tidyverse)
library(ggplot2)
library(ggsci) # Load ggsci
# library(patchwork) # For arranging multiple plots
# library(here) # Optional: here::here("path/to/your/file.txt")

# --- 1. Load and Prepare Data ---
# Define the path to your combined sequencing summary file from Part 1b of the bash script
# This file should have a "TRMT" column and the metrics of interest for overall treatment reads.
combined_seqsum_file <- "DSC2025_RAD_Enrich.reads.all_trmts.combined.seqsum.txt"

# Check if the file exists
if (!file.exists(combined_seqsum_file)) {
  stop(paste("Combined treatment sequencing summary file not found:", combined_seqsum_file, 
             "\nPlease ensure the file path is correct and the file exists in your working directory."))
}

# Read the TSV file 
seqsum_data <- read_tsv(combined_seqsum_file, show_col_types = FALSE)

# Inspect the data
print("--- Raw Combined Treatment Sequencing Summary Data ---")
print(head(seqsum_data))
print(colnames(seqsum_data)) # Verify column names

# Data Wrangling for Treatment-Level Data:
seqsum_data_wrangled <- seqsum_data %>%
  mutate(
    TRMT = factor(TRMT, levels = c("Control", "AS")),
    sequence_length_template = as.numeric(sequence_length_template),
    mean_qscore_template = as.numeric(mean_qscore_template)
  )

# Check for NAs that might have been introduced if columns didn't exist or conversion failed
if(any(is.na(seqsum_data_wrangled$sequence_length_template))) {
  warning("NAs introduced in 'sequence_length_template' for treatment-level data. Check original column name and data type.")
}
if(any(is.na(seqsum_data_wrangled$mean_qscore_template))) {
  warning("NAs introduced in 'mean_qscore_template' for treatment-level data. Check original column name and data type.")
}
if(any(is.na(seqsum_data_wrangled$TRMT))) {
  warning("NAs introduced in 'TRMT' or levels don't match for treatment-level data. Check 'TRMT' column in your data.")
}


print("--- Wrangled Combined Treatment Sequencing Summary Data ---")
print(head(seqsum_data_wrangled))
print(summary(seqsum_data_wrangled$TRMT)) 

# --- 2. Visualize Treatment-Level Data with Violin Plots ---

# Plot 1: Violin plot of Sequence Length Template by Treatment
plot_seq_length_violin_trmt <- ggplot(seqsum_data_wrangled, aes(x = TRMT, y = sequence_length_template, fill = TRMT)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") + 
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.shape = NA) + 
  scale_fill_npg() + 
  scale_y_log10(labels = scales::comma) + 
  labs(title = "Distribution of Read Lengths (Template) by Treatment",
       x = "Treatment",
       y = "Sequence Length (bp, log10 scale)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_seq_length_violin_trmt)
ggsave("plot_trmt_seq_length_violin.png", plot = plot_seq_length_violin_trmt, width = 7, height = 6)


# Plot 2: Violin plot of Mean Q-score Template by Treatment
plot_mean_qscore_violin_trmt <- ggplot(seqsum_data_wrangled, aes(x = TRMT, y = mean_qscore_template, fill = TRMT)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.shape = NA) +
  scale_fill_npg() + 
  labs(title = "Distribution of Mean Read Quality (Template) by Treatment",
       x = "Treatment",
       y = "Mean Q-score (Template)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_mean_qscore_violin_trmt)
ggsave("plot_trmt_mean_qscore_violin.png", plot = plot_mean_qscore_violin_trmt, width = 7, height = 6)

print("--- Treatment-Level Violin Plots Finished ---")
print("---                                       ---")


# --- Part 3: Visualize Mapped Isolate-Specific Data with Violin Plots ---
# This part loads the combined summary of reads mapped to specific isolates
# and plots distributions of length and quality for each isolate, split by treatment.

# --- 3.1: Load and Prepare Mapped Isolate Sequencing Summary Data ---
# Define the path to your combined *mapped isolate* sequencing summary file.
# This file should have "TRMT", "ISO", and the metrics of interest.
# Filename confirmed by user.
mapped_iso_seqsum_file <- "DSC2025_RAD_Enrich.mapped.all_trmts.all_isos.combined.seqsum.txt" 

# Check if the file exists
if (!file.exists(mapped_iso_seqsum_file)) {
  stop(paste("Combined mapped isolate sequencing summary file not found:", mapped_iso_seqsum_file, 
             "\nPlease ensure the file path is correct and the file exists in your working directory.",
             "\nThis file is typically generated in Part 3 of your bash pipeline."))
}

# Read the TSV file
mapped_iso_data <- read_tsv(mapped_iso_seqsum_file, show_col_types = FALSE)

print("--- Raw Combined Mapped Isolate Sequencing Summary Data ---")
print(head(mapped_iso_data))
print(colnames(mapped_iso_data))

# Data Wrangling for Mapped Isolate Data:
# Ensure TRMT and ISO are factors, and target columns are numeric.
mapped_iso_data_wrangled <- mapped_iso_data %>%
  mutate(
    TRMT = factor(TRMT, levels = c("Control", "AS")),
    ISO = factor(ISO), # Convert ISO to factor for faceting
    sequence_length_template = as.numeric(sequence_length_template),
    mean_qscore_template = as.numeric(mean_qscore_template)
  )

# Check for NAs
if(any(is.na(mapped_iso_data_wrangled$sequence_length_template))) {
  warning("NAs introduced in 'sequence_length_template' for mapped isolate data. Check original column name and data type.")
}
if(any(is.na(mapped_iso_data_wrangled$mean_qscore_template))) {
  warning("NAs introduced in 'mean_qscore_template' for mapped isolate data. Check original column name and data type.")
}
if(any(is.na(mapped_iso_data_wrangled$TRMT))) {
  warning("NAs introduced in 'TRMT' or levels don't match for mapped isolate data.")
}
if(any(is.na(mapped_iso_data_wrangled$ISO))) {
  warning("NAs introduced in 'ISO' for mapped isolate data.")
}

print("--- Wrangled Combined Mapped Isolate Sequencing Summary Data ---")
print(head(mapped_iso_data_wrangled))
print(summary(mapped_iso_data_wrangled$TRMT))
print(summary(mapped_iso_data_wrangled$ISO))


# --- 3.2: Visualize Mapped Isolate Data with Violin Plots (Faceted by Isolate) ---

# Plot 3.1: Violin plot of Sequence Length Template by Treatment, Faceted by Isolate
plot_iso_seq_length_violin <- ggplot(mapped_iso_data_wrangled, 
                                     aes(x = TRMT, y = sequence_length_template, fill = TRMT)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.shape = NA) +
  scale_fill_npg() +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~ ISO, ncol = 3, scales = "free_y") + # Facet by Isolate, allow y-axis to vary
  labs(title = "Distribution of Mapped Read Lengths by Isolate and Treatment",
       x = "Treatment",
       y = "Sequence Length (bp, log10 scale)") +
  theme_minimal(base_size = 12) + # Adjust base size for facets
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels if crowded

print(plot_iso_seq_length_violin)
ggsave("plot_iso_seq_length_violin_faceted.png", plot = plot_iso_seq_length_violin, width = 10, height = 8)


# Plot 3.2: Violin plot of Mean Q-score Template by Treatment, Faceted by Isolate
plot_iso_mean_qscore_violin <- ggplot(mapped_iso_data_wrangled, 
                                       aes(x = TRMT, y = mean_qscore_template, fill = TRMT)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.shape = NA) +
  scale_fill_npg() +
  facet_wrap(~ ISO, ncol = 3, scales = "free_y") + # Facet by Isolate
  labs(title = "Distribution of Mean Mapped Read Quality by Isolate and Treatment",
       x = "Treatment",
       y = "Mean Q-score (Template)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_iso_mean_qscore_violin)
ggsave("plot_iso_mean_qscore_violin_faceted.png", plot = plot_iso_mean_qscore_violin, width = 10, height = 8)


print("--- R Script for Violin Plots (including Isolate-Specific) Finished ---")
print(paste("Plots saved to:", getwd())) 
