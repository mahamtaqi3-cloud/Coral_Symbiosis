# Install BiocManager if you haven't already
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install the core bioinformatics packages needed for your project
BiocManager::install(c("DESeq2", "ggplot2", "dplyr", "pheatmap", "clusterProfiler"))
options(download.file.method = "libcurl")
BiocManager::install("DESeq2")
library(DESeq2)
# Install the GEOquery package to fetch data from NCBI
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
# Load the library
library(GEOquery)
# Download the study metadata
gse <- getGEO("GSE97613", GSEMatrix = TRUE)
show(gse)
# Access the first dataset in the list
eset <- gse[[1]]
# Get the expression data (this is the table of gene values)
expression_data <- exprs(eset)
# View the first few rows to see what the data looks like
head(expression_data)
# Get the sample information (e.g., conditions/treatments)
sample_info <- pData(eset)
# View the metadata to identify which columns define your "Control" vs "Heat Stress" groups
View(sample_info)
# Example: Defining your experimental design
# Create a data frame where 'condition' is your variable (e.g., "treated" vs "not treated")
sample_info <- data.frame(
  row.names = colnames(your_count_matrix),
  condition = factor(c("not_treated", "not_treated", "treated", "treated"))
)
# Use the data we already extracted earlier
sample_info <- data.frame(
  row.names = colnames(expression_data),
  condition = factor(c("not_treated", "not_treated", "treated", "treated", "treated", "treated"))
)
# Load ggplot2 for plotting
library(ggplot2)
# Create a simple scatter plot of your samples
# We'll use the column names (sample IDs) and their condition
pca_data <- data.frame(sample = rownames(sample_info), condition = sample_info$condition)
ggplot(pca_data, aes(x = sample, y = condition, color = condition)) +
  geom_point(size = 5) +
  theme_minimal() +
  labs(title = "Experimental Design Overview")
# Convert expression_data to a matrix
count_matrix <- round(as.matrix(expression_data))
# Load DESeq2
library(DESeq2)
# Create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)
# Install and load the necessary package
install.packages("readxl")
library(readxl)
# Read your file (ensure the file is in your project folder)
# Replace 'GSE141133_TMM_Mg.xlsx' with your exact file name
data <- read_excel("GSE141133_TMM_Mg.xlsx")
getwd()
library(readxl)
data <- read_excel("GSE141133_TMM_Mg.xlsx")
# Convert data to a matrix, using the first column as row names
# Note: Check if the first column contains your Gene IDs
mat <- as.matrix(data[, -1])
rownames(mat) <- data[[1]]
# Ensure pheatmap is installed
if(!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)
# Plot the heatmap
pheatmap(mat, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
# Assuming the first column is your Gene IDs, select only columns 2 through 7
mat_numeric <- as.matrix(data[, 2:7])
# Assign the Gene IDs from the first column to the row names of your new matrix
rownames(mat_numeric) <- data[[1]]
# Now try the heatmap again with the strictly numeric matrix
library(pheatmap)
pheatmap(mat_numeric, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
# Create a clean numeric matrix
mat_numeric <- as.matrix(data[, 2:7])
# Force all values to be treated as numbers
storage.mode(mat_numeric) <- "numeric"
# Re-apply the row names
rownames(mat_numeric) <- data[[1]]
library(pheatmap)
pheatmap(mat_numeric, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
dim(mat_clean)
# Now use the cleaned matrix for the heatmap
pheatmap(mat_clean, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
# Now use the cleaned matrix for the heatmap
pheatmap(mat_clean, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
mat_clean <- mat_numeric[complete.cases(mat_numeric), ]
library(readxl)
data <- read_excel("GSE141133_TMM_Mg.xlsx")
mat_numeric <- as.matrix(data[, 2:7])
library(readxl)
data <- read_excel("GSE141133_TMM_Mg.xlsx")
mat_numeric <- as.matrix(data[, 2:7])
# 1. Load the data
library(readxl)
data <- read_excel("GSE141133_TMM_Mg.xlsx")
# 2. Select columns 2 through 7 and convert to matrix in one step
mat_numeric <- data.matrix(data[, 2:7])
# 2. Select columns 2 through 7 and convert to matrix in one step
mat_numeric <- data.matrix(data[, 2:7])
# 3. Add row names
rownames(mat_numeric) <- data[[1]]
# 4. Clean NAs and plot
mat_clean <- mat_numeric[complete.cases(mat_numeric), ]
library(pheatmap)
pheatmap(mat_clean, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
# Replace NA, NaN, and Inf with 0
mat_clean[is.na(mat_clean)] <- 0
mat_clean[is.nan(mat_clean)] <- 0
mat_clean[is.infinite(mat_clean)] <- 0
library(pheatmap)
pheatmap(mat_clean, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
rm(data, mat, mat_numeric)
gc()
mat_subset <- mat_clean[1:500, ]
mat_subset
memory.limit(size = 8000)
library(pheatmap)
pheatmap(mat_subset, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
mat_subset[is.na(mat_subset)] <- 0
mat_subset[is.nan(mat_subset)] <- 0
mat_subset[is.infinite(mat_subset)] <- 0
library(pheatmap)
pheatmap(mat_subset, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
# Keep only rows where all values are not the same (sd > 0)
mat_subset <- mat_subset[apply(mat_subset, 1, sd) > 0, ]
# Ensure there are no NAs left
mat_subset <- mat_subset[complete.cases(mat_subset), ]
library(pheatmap)
pheatmap(mat_subset, scale = "row", show_rownames = FALSE, main = "Expression Heatmap")
