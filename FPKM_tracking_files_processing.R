
### Script to process FPKM tracking files to make counts matrix ---------------------------------------------------------

# FPKM tracking files are generated from Cufflinks.
# Cufflinks is a bioinformatics tool designed to analyze transcriptome sequencing (RNA-Seq) data to estimate transcript abundances, identify novel isoforms, and perform differential expression analysis.

# Dataset link -
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-71493

setwd("C:/Ankur_IIITD/Bioinfo_Basics/microarray_raw_data_processing/E_GEOD_71493_tumor")

# Load required packages
library(dplyr)
library(readr)

# List all FPKM tracking files
file_list <- list.files(pattern = "*.fpkm_tracking", full.names = TRUE)

# Extract sample IDs 
sample_ids <- gsub("/.fpkm_tracking", "", basename(file_list))
sample_ids <- gsub("_.+", "", sample_ids)   # Keep only GSM ID part

# Initialize an empty list to store data
expression_list <- list()

# Loop through files and extract FPKM values
for (i in seq_along(file_list)) {
  file <- file_list[i]
  sample_id <- sample_ids[i]
  
  # Read the data
  df <- read_tsv(file, col_types = cols()) %>%
    select(gene_short_name, FPKM) %>%
    distinct(gene_short_name, .keep_all = TRUE) %>%   # Ensure unique genes
    rename(!!sample_id := FPKM)
  
  expression_list[[sample_id]] <- df
}

# Merge all data frames by gene_short_name
expression_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_short_name"), expression_list)
expression_matrix <- as.data.frame(expression_matrix)

# Set gene_short_name as row names
rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix <- expression_matrix[ , -1] 

# Save the expression matrix
write.csv(expression_matrix, "tumor_expression_matrix.csv", row.names = TRUE)

# Print first few rows
head(expression_matrix)
