
### Script to process raw text files to make expression matrix -----------------------------------------------------

# This data is generated from transcription profiling by array.

# Dataset link -
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MAXD-9


# Load necessary library
library(dplyr)
library(biomaRt)

setwd("C:/Ankur_IIITD/Bioinfo_Basics/microarray_raw_data_processing/E_MAXD-9_Paca2_Cells")

# Function to read and process each file
process_file <- function(file_path) {
  df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Print first few values to check formatting
  print("Original Reporter.identifier values:")
  print(head(df$`Reporter.identifier`))
  
  # Extract the gene identifier
  df$`Reporter.identifier` <- gsub("^.*\\[Reporter\\]([^/]+)/.*$", "\\1", df$`Reporter.identifier`)
  
  # Check if extraction is working correctly
  print("Extracted Reporter.identifier values:")
  print(head(df$`Reporter.identifier`))
  
  # Ensure the column exists before selection
  if (!"Reporter.identifier" %in% colnames(df)) {
    stop("Error: 'Reporter.identifier' column not found after regex extraction.")
  }
  
  # Ensure expression column exists
  expr_col <- grep("GenePix", colnames(df), value = TRUE)
  if (length(expr_col) == 0) {
    stop("Error: Expression column ('GenePix.F635.Median...B635') not found in file.")
  }
  
  # Select relevant columns
  df_selected <- df %>% dplyr::select(`Reporter.identifier`, all_of(expr_col))
  
  # Rename expression column based on the sample name
  sample_name <- gsub(".txt", "", basename(file_path))
  colnames(df_selected)[2] <- sample_name
  
  # Aggregate duplicate genes by taking the mean
  df_selected <- df_selected %>%
    group_by(`Reporter.identifier`) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
  
  return(df_selected)
}

# List all .txt files in the directory
file_list <- list.files(pattern = "*.txt")
file_list

# Read and process all files
expr_list <- lapply(file_list, process_file)

# Merge all data frames by "Reporter.identifier"
expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "Reporter.identifier"), expr_list)    # it will create tibble

# Use Ensembl's biomaRt to fetch gene symbols for RefSeq Protein IDs
convert_refseq_to_symbol <- function(refseq_ids) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Connect to Ensembl database
  
  # Query the database to get Gene Symbols
  gene_map <- getBM(
    attributes = c("refseq_mrna", "hgnc_symbol"),  # RefSeq Protein ID -> Gene Symbol
    filters = "refseq_mrna",
    values = refseq_ids,
    mart = mart
  )
  
  return(gene_map)
}

# Extract the first column from expr_matrix (RefSeq Protein IDs)
refseq_ids <- expr_matrix[,1]

# Convert RefSeq IDs to Gene Symbols
gene_mapping <- convert_refseq_to_symbol(refseq_ids)

# Merge gene symbols with expr_matrix
expr_matrix <- expr_matrix %>%
  left_join(gene_mapping, by = c("Reporter.identifier" = "refseq_mrna"))

# Convert Reporter.identifier to gene symbols where available
expr_matrix <- expr_matrix %>%
  mutate(Reporter.identifier = ifelse(!is.na(hgnc_symbol), hgnc_symbol, Reporter.identifier)) %>%
  dplyr::select(-hgnc_symbol)  # Explicitly use dplyr::select()

# Save the expression matrix to a CSV file
write.csv(expr_matrix, "expression_matrix.csv", row.names = TRUE)


