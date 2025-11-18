
# 1. Load Libraries
library(GEOquery)
library(limma)
library(dplyr)      # Used for metadata manipulation
library(ggplot2)

# 2. Load Dataset
gse <- getGEO("GSE6740", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])


# 3. Check for Missing Values
cat("Missing values in expression data:", sum(is.na(expression_data)), "\n")


# 4. Boxplot Before Normalization
png(file = file.path("GSE6740_boxplot_before_normalization.png"))
boxplot(expression_data, main = "GSE6740 - Before Normalization", las = 2, outline = FALSE)
dev.off()


# 5. Apply Normalization
expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")


# 6. Boxplot After Normalization
png(file = file.path("GSE6740_boxplot_after_normalization.png"))
boxplot(expression_data, main = "GSE6740 - After Normalization", las = 2, outline = FALSE)
dev.off()

# 7. Filter Low-Expressed Genes
gene_means <- rowMeans(expression_data)
threshold <- quantile(gene_means, probs = 0.25)
expression_data <- expression_data[gene_means > threshold, ]
cat("Dimensions after filtering:", dim(expression_data), "\n")



# 8. Extract Subtype and Cell Type (Modified)
# Extract Cell Type (CD4 or CD8)
cell_type <- ifelse(grepl("CD4", metadata$title), "CD4+", "CD8+")
cell_type <- as.factor(cell_type)
# Extract Patient Group (A, C, L, N)
patient_group <- dplyr::case_when(
  grepl("Acute", metadata$title) ~ "Acute (A)",
  grepl("Chronic", metadata$title) ~ "Chronic (C)",
  grepl("non-progressive", metadata$title) ~ "Non-Progressor (L)",
  grepl("Uninfected", metadata$title) ~ "Uninfected (N)",
  TRUE ~ "Unknown" # Fallback
)
patient_group <- as.factor(patient_group)
# Print summaries of the new factors
cat("Cell Type distribution:\n")
table(cell_type)
cat("\nPatient Group (Subtype) distribution:\n")
table(patient_group)

# 9. Generate PCA Plot (Modified)
pca_result <- prcomp(t(expression_data), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[, 1], 
  PC2 = pca_result$x[, 2], 
  Patient_Group = patient_group,
  Cell_Type = cell_type
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Patient_Group, shape = Cell_Type)) +
  geom_point(size = 3) + # Made points slightly larger
  ggtitle("GSE6740 PCA: Clustering by Patient Group and Cell Type") +
  theme_minimal() +
  labs(color = "Patient Group", shape = "Cell Type")

ggsave(file.path("GSE6740_pca_plot.png"), plot = pca_plot)


# 10. Save Processed Data (Modified)
save(expression_data, metadata, patient_group, cell_type, 
     file = file.path("GSE6740_processed_data_step1.RData"))
