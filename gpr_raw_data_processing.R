
### Script to process .gpr files to make expression matrix -------------------------------------------------------------------

# .gpr file contains data from a miRNA microarray experiment designed to measure the expression levels of specific human miRNAs in a biological sample. 

# Dataset link -
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3628

# miRNA dataset link -
# https://www.mirbase.org/download/


###################################### Brief information about the data ##########################################

# 1. Spot Identification & Position
# Block – The block number where the spot is located on the array.
# Column – The column number within the block.
# Row – The row number within the block.
# Name – The gene or probe name.
# ID – The sequence ID associated with the probe.
# X – The X-coordinate of the spot in the scanned image.
# Y – The Y-coordinate of the spot in the scanned image.
# Dia. – The diameter of the spot in microns.

# 2. Fluorescence Intensity at 635 nm (Red Channel)
# F635 Median – Median fluorescence intensity of the spot.
# F635 Mean – Mean fluorescence intensity of the spot.
# F635 SD – Standard deviation of fluorescence intensity within the spot.
# F635 CV – Coefficient of variation (SD / Mean) of the fluorescence intensity.

# 3. Background Fluorescence (635 nm)
# B635 – Background intensity of the spot.
# B635 Median – Median background intensity of the local region.
# B635 Mean – Mean background intensity of the local region.
# B635 SD – Standard deviation of background intensity.
# B635 CV – Coefficient of variation for the background.

# 4. Quality Control Metrics
# % > B635+1SD – Percentage of pixels in the spot with intensity greater than background mean + 1 standard deviation.
# % > B635+2SD – Percentage of pixels in the spot with intensity greater than background mean + 2 standard deviations.
# F635 % Sat. – Percentage of pixels in the spot that are saturated (i.e., at maximum detector value).

# 5. Ratio-Based Measurements (635 nm / 532 nm)
# Ratio of Medians (635/2) – Ratio of median intensities at 635 nm vs. 532 nm.
# Ratio of Means (635/2) – Ratio of mean intensities at 635 nm vs. 532 nm.
# Median of Ratios (635/2) – Median of individual pixel ratios.
# Mean of Ratios (635/2) – Mean of individual pixel ratios.
# Ratios SD (635/2) – Standard deviation of the ratio values.
# Rgn Ratio (635/2) – Region-based intensity ratio.
# Rgn R2 (635/2) – Goodness-of-fit (R²) for the ratio calculation.

# 6. Spot Area & Shape
# F Pixels – Number of pixels used for the foreground (spot).
# B Pixels – Number of pixels used for the background calculation.
# Circularity – Shape factor measuring how circular the spot is (1 = perfect circle).

# 7. Summed Intensities
# Sum of Medians (635/2) – Sum of median intensities at 635 nm and 532 nm.
# Sum of Means (635/2) – Sum of mean intensities at 635 nm and 532 nm.
# Log Ratio (635/2) – Log-transformed ratio of median intensities.

# 8. Processed Expression Values
# F635 Median - B635 – Background-subtracted median fluorescence intensity.
# F635 Mean - B635 – Background-subtracted mean fluorescence intensity.
# F635 Total Intensity – Total fluorescence intensity for the spot.

# 9. Signal-to-Noise & Flags
# SNR 635 – Signal-to-noise ratio at 635 nm (F635 Mean / B635 SD).
# Flags – Flagged values indicating spot quality issues (e.g., poor-quality spots).
# Normalize – Placeholder for normalization status (typically unused).
# Autoflag – Automatically flagged spots based on quality control criteria.


#######################################################################################################################


# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)

setwd("C:/Ankur_IIITD/Bioinfo_Basics/microarray_raw_data_processing/E_MTAB_3628_miRNA_Profile")

# List all .gpr files in the working directory
gpr_files <- list.files(pattern = "*.gpr")
gpr_files

# Initialize an empty list to store data from each GPR file
data_list <- list()

for (file in gpr_files) {
  raw_data <- readLines(file)   # read full content
  header_line <- grep("Block", raw_data)   # search for the header line
  df <- read.delim(file, skip = header_line - 1, header = TRUE, check.names = FALSE)   # read data from correct point
  
  # Ensure required columns exist
  if (!("Name" %in% colnames(df)) || !("F635 Median - B635" %in% colnames(df))) {
    stop(paste("Required columns not found in file:", file))
  }
  
  # Select relevant columns
  df_selected <- df[, c("Name", "F635 Median - B635")]
  colnames(df_selected) <- c("miRNA", gsub(".gpr", "", basename(file)))
  
  # Aggregate values by miRNA (taking the mean in case of duplicates)
  df_selected <- df_selected %>%
    group_by(miRNA) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Store in list
  data_list[[file]] <- df_selected
}

# Merge all data frames by miRNA column
expression_matrix <- Reduce(function(x, y) full_join(x, y, by = "miRNA"), data_list)

# Replace NA values with 0 (optional)
expression_matrix[is.na(expression_matrix)] <- 0



################## Probe IDs to Gene IDs Mapping (Additional) #####################

# Read mirna data file
mirna_file <- readLines("miRNA.dat")

# Extract miRNA IDs
miRNA_IDs <- grep("^ID", mirna_file, value = TRUE) %>%
  str_extract("ID\\s+\\S+") %>%
  gsub("ID\\s+", "", .)

# Extract Accession IDs
access_num <- grep("^AC", mirna_file, value = TRUE) %>%
  str_extract("AC\\s+\\S+") %>%
  gsub("AC\\s+", "", .)%>%
  gsub(";", "", .)

# Create a data frame with mappings
mirna_mapping <- data.frame(Probe_ID = miRNA_IDs, Accession = access_num, stringsAsFactors = FALSE)

# Function to replace miRNA probe IDs with their corresponding accession numbers
replace_probes_with_accession <- function(expression_matrix, mirna_mapping) {
  expression_matrix <- expression_matrix %>%
    mutate(miRNA = sapply(miRNA, function(x) {
      # Check if any Probe_ID in mirna_mapping is a prefix of x (ignoring suffixes)
      matched_probe <- mirna_mapping$Probe_ID[str_detect(x, paste0("^", mirna_mapping$Probe_ID, "\\b"))]
      
      # Get the corresponding Accession ID
      matched_accession <- mirna_mapping$Accession[mirna_mapping$Probe_ID == matched_probe]
      
      # If a match is found, return the accession ID, otherwise keep the original value
      ifelse(length(matched_accession) > 0, matched_accession, x)
    }))
  return(expression_matrix)
}

# Apply the function to update expression_matrix
expression_matrix <- replace_probes_with_accession(expression_matrix, mirna_mapping)









