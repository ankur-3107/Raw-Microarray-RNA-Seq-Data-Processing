
### Script to process CEL data files to make counts data (sample * genes) -------------------------------------------------------

# .CEL raw data files are generated from affymetrix microarray
# To process .CEL files, oligo and affy packages can be use
# For annotations hgu133plus2.db package can be use

## Dataset link -
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-25066

setwd("C:/Ankur_IIITD/Bioinfo_Basics/microarray_raw_data_processing/E_GEOD_25066_breast_cancer")

# Install the required packages
# BiocManager::install(c("affy", "oligo", "Biobase", "annotate", "hgu133plus2.db"))

library(affy)  # For reading CEL files
library(oligo) # For processing Affymetrix data
library(Biobase)  # to implement base functions
library(annotate)  # to summarize genomic annotations  
library(hgu133plus2.db) # Annotation package for the microarray platform

# List all CEL files
cel.files <- list.files(pattern = "*.CEL")

# Read all files at a time
raw.data <- read.celfiles(cel.files)

# RMA Normalization (Robust Multi-array Analysis)
rma.data <- rma(raw.data)

# Extract expression matrix
expr.data <- exprs(rma.data)

# Modify column names to keep only sample IDs
colnames(expr.data) <- sub("_.*", "", colnames(expr.data))

# Map probe IDs to gene symbols
probe.ids <- rownames(expr.data)
gene.symbols <- mapIds(hgu133plus2.db,
                       keys = probe.ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first"
                       )

# Set gene symbols as row names
rownames(expr.data) <- gene.symbols

# Remove probes without gene symbols
expr.data <- na.omit(expr.data)
head(expr.data)

write.csv(expr.data, "hESC_expression_data.csv", row.names = TRUE)
