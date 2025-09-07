# R Script for Visualizing SeqKit Stats from Adaptive Sampling Pipeline

# Install packages if you haven't already
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("ggsci") # For scientific journal color palettes
# install.packages("here") # Optional, for easier path management

library(tidyverse)
library(ggplot2)
library(ggsci) # Load ggsci
# library(here) # Optional: here::here("path/to/your/file.tsv")

# --- Part 1: Analysis of Treatment-Specific FASTQ Stats ---
# This part focuses on the overall statistics for reads filtered by treatment (AS vs. Control)

# --- 1.1: Load and Prepare Treatment Stats Data ---
# Define the path to your treatment stats file.
# Replace with your actual file name/path if different.
# Assumes the file is in the same directory as this R script.
treatment_stats_file <- "DSC2025_RAD_Enrich.reads.all_treatments.seqkit_stats.tsv"

# Check if the file exists
if (!file.exists(treatment_stats_file)) {
  stop(paste("Treatment stats file not found:", treatment_stats_file, 
             "\nPlease ensure the file path is correct and the file exists in your working directory."))
}

# Read the TSV file
treatment_data <- read_tsv(treatment_stats_file, show_col_types = FALSE)

# Inspect the data
print("--- Raw Treatment Stats Data ---")
print(head(treatment_data))

# Data Wrangling: Extract SAMPLEID and Treatment (TRMT) from the 'file' column
# Example filename format: DSC2025_RAD_Enrich.reads.AS.fastq.gz
treatment_data_wrangled <- treatment_data %>%
  mutate(
    # Extract SAMPLEID (everything before ".reads.")
    SAMPLEID = gsub("\\.reads\\..*", "", file),
    # Extract Treatment (AS or Control)
    TRMT = str_extract(file, "(?<=reads\\.)(AS|Control)(?=\\.fastq\\.gz)"),
    # Convert TRMT to factor to control plotting order
    TRMT = factor(TRMT, levels = c("Control", "AS")) 
  ) %>%
  # Convert relevant columns to numeric if they are not already (SeqKit -T should handle this)
  mutate(
    num_seqs = as.numeric(num_seqs),
    sum_len = as.numeric(sum_len),
    N50 = as.numeric(N50),
    avg_qual = as.numeric(AvgQual) # Corrected: seqkit column is avg_qual
  )

# Rename avg_qual to AveQual for consistency with the user's request
if ("avg_qual" %in% colnames(treatment_data_wrangled)) {
  treatment_data_wrangled <- treatment_data_wrangled %>%
    rename(AveQual = avg_qual)
}


print("--- Wrangled Treatment Stats Data ---")
print(head(treatment_data_wrangled))
print(summary(treatment_data_wrangled$TRMT)) # Check factor levels

# --- 1.2: Visualize Treatment-Specific Differences ---

# Plot 1: Number of Reads (num_seqs) by Treatment
plot_num_seqs_trmt <- ggplot(treatment_data_wrangled, aes(x = TRMT, y = num_seqs, fill = TRMT)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_npg() + # Apply ggsci NPG palette
  labs(title = "Number of Reads by Treatment",
       x = "Treatment",
       y = "Number of Reads (num_seqs)") +
  scale_y_continuous(labels = scales::comma) + # Format y-axis labels
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_num_seqs_trmt)
ggsave("plot_num_seqs_by_treatment.png", plot = plot_num_seqs_trmt, width = 6, height = 5)


# Plot 2: Total Yield (sum_len) by Treatment
plot_sum_len_trmt <- ggplot(treatment_data_wrangled, aes(x = TRMT, y = sum_len, fill = TRMT)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_npg() + # Apply ggsci NPG palette
  labs(title = "Total Yield (Sum of Lengths) by Treatment",
       x = "Treatment",
       y = "Total Yield (bp)") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_sum_len_trmt)
ggsave("plot_sum_len_by_treatment.png", plot = plot_sum_len_trmt, width = 6, height = 5)


# Plot 3: N50 by Treatment
plot_n50_trmt <- ggplot(treatment_data_wrangled, aes(x = TRMT, y = N50, fill = TRMT)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  scale_fill_npg() + # Apply ggsci NPG palette
  labs(title = "N50 by Treatment",
       x = "Treatment",
       y = "N50 (bp)") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_n50_trmt)
ggsave("plot_n50_by_treatment.png", plot = plot_n50_trmt, width = 6, height = 5)


# Plot 4: Average Quality (AveQual) by Treatment
plot_avg_qual_trmt <- ggplot(treatment_data_wrangled, aes(x = TRMT, y = AveQual, fill = TRMT)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  scale_fill_npg() + # Apply ggsci NPG palette
  labs(title = "Average Read Quality by Treatment",
       x = "Treatment",
       y = "Average Quality (AveQual)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_avg_qual_trmt)
ggsave("plot_avg_qual_by_treatment.png", plot = plot_avg_qual_trmt, width = 6, height = 5)


# --- Part 2: Analysis of Mapped Isolate-Specific FASTQ Stats ---
# This part focuses on the proportion of yield (sum_len) mapping to each isolate, compared between treatments.

# --- 2.1: Load and Prepare Mapped Isolate Stats Data ---
# Define the path to your mapped isolate stats file.
isolate_stats_file <- "DSC2025_RAD_Enrich.mapped.all_isolates.seqkit_stats.tsv"

# Check if the file exists
if (!file.exists(isolate_stats_file)) {
  stop(paste("Isolate stats file not found:", isolate_stats_file,
             "\nPlease ensure the file path is correct and the file exists in your working directory."))
}

# Read the TSV file
isolate_data <- read_tsv(isolate_stats_file, show_col_types = FALSE)

print("--- Raw Mapped Isolate Stats Data ---")
print(head(isolate_data))

# Data Wrangling: Extract SAMPLEID, Treatment (TRMT), and Isolate (ISO) from the 'file' column
# Example filename format: DSC2025_RAD_Enrich.mapped.AS.Bacillus.fastq.gz
isolate_data_wrangled <- isolate_data %>%
  extract(
    file, 
    into = c("SAMPLEID_temp", "TRMT", "ISO", "ext"), 
    regex = "^(.*?)\\.mapped\\.(AS|Control)\\.([^.]+?)\\.(fastq\\.gz)$",
    remove = FALSE 
  ) %>%
  select(-ext, -SAMPLEID_temp) %>% 
  mutate(
    SAMPLEID = gsub("\\.mapped\\..*", "", file),
    TRMT = factor(TRMT, levels = c("Control", "AS")), # Convert TRMT to factor for plotting order
    num_seqs = as.numeric(num_seqs), # Keep num_seqs if needed for other plots
    sum_len = as.numeric(sum_len)   # Ensure sum_len is numeric for yield calculations
  )

print("--- Wrangled Mapped Isolate Stats Data (for Yield Analysis) ---")
print(head(isolate_data_wrangled))
print(summary(isolate_data_wrangled$TRMT)) # Check factor levels


# Calculate total yield (sum_len) per treatment (for calculating proportions based on yield)
total_yield_per_treatment <- isolate_data_wrangled %>%
  group_by(TRMT) %>%
  summarise(total_trmt_yield = sum(sum_len, na.rm = TRUE), .groups = "drop")

print("--- Total Yield per Treatment (Mapped Isolates) ---")
print(total_yield_per_treatment)

# Join total yield back and calculate proportions based on yield
isolate_yield_proportions <- isolate_data_wrangled %>%
  left_join(total_yield_per_treatment, by = "TRMT") %>%
  mutate(proportion_yield = if_else(total_trmt_yield > 0, sum_len / total_trmt_yield, 0))
  # Using if_else to handle cases where total_trmt_yield might be 0 to avoid NaN, setting proportion to 0.

print("--- Mapped Isolate Yield Proportions ---")
print(head(isolate_yield_proportions))

# --- 2.2: Visualize Mapped Isolate-Specific Differences (Proportions by Yield) ---

# Plot 5: Proportion of Mapped Yield per Isolate by Treatment
plot_iso_yield_proportions <- ggplot(isolate_yield_proportions, aes(x = TRMT, y = proportion_yield, fill = ISO)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_npg() + # Apply ggsci NPG palette (handles 7+ colors)
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of Mapped Yield by Isolate and Treatment",
       x = "Treatment",
       y = "Proportion of Total Mapped Yield for Treatment",
       fill = "Isolate") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_iso_yield_proportions)
ggsave("plot_isolate_yield_proportions_by_treatment.png", plot = plot_iso_yield_proportions, width = 8, height = 6)

# Optional: Faceted bar plot for easier comparison of individual isolates' yield proportions
plot_iso_yield_proportions_faceted <- ggplot(isolate_yield_proportions, aes(x = TRMT, y = proportion_yield, fill = TRMT)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_npg() + # Apply ggsci NPG palette
  facet_wrap(~ ISO, scales = "free_y", ncol = 3) + # Free_y allows y-axis to adapt per isolate
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of Mapped Yield by Isolate (Faceted by Isolate)",
       x = "Treatment",
       y = "Proportion of Total Mapped Yield for Treatment",
       fill = "Treatment") +
  theme_minimal(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none", 
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_iso_yield_proportions_faceted)
ggsave("plot_isolate_yield_proportions_faceted.png", plot = plot_iso_yield_proportions_faceted, width = 10, height = 7)


# Optional: Compare absolute total yield per isolate (sum_len)
plot_iso_yield_absolute <- ggplot(isolate_data_wrangled, aes(x = ISO, y = sum_len, fill = TRMT)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_npg() + # Apply ggsci NPG palette
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Total Mapped Yield by Isolate and Treatment",
       x = "Isolate",
       y = "Total Yield (bp)",
       fill = "Treatment") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_iso_yield_absolute)
ggsave("plot_isolate_yield_absolute_by_treatment.png", plot = plot_iso_yield_absolute, width = 10, height = 6)


print("--- R Script for Visualization Finished ---")
print(paste("Plots saved to:", getwd())) # Print current working directory where plots are saved

