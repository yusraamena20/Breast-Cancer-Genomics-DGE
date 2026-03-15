# ------------------------------------------------------------------------------
# Project: Differential Gene Expression (DGE) Analysis of Breast Cancer
# Dataset: GSE15852 (NCBI GEO)
# Objective: Identify top DEGs between Tumor and Normal breast tissue.
# ------------------------------------------------------------------------------

# 1. Load Necessary Libraries
# Note: If missing, install via BiocManager::install(c("GEOquery", "limma", "ggplot2", "pheatmap", "org.Hs.eg.db"))
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)

# 2. Load and Preprocess Data
# Using comment.char="!" to skip GEO-specific metadata headers
data <- read.delim("GSE15852_series_matrix.txt", header = TRUE, row.names = 1, comment.char = "!")

# Log2 Transformation (Crucial for normalizing skewed intensity values)
data_log <- log2(data + 1)

# 3. Define Experimental Groups
# 43 Tumor samples followed by 43 Normal samples
group <- factor(c(rep("Tumor", 43), rep("Normal", 43)))

# 4. Statistical Analysis (limma Pipeline)
design <- model.matrix(~0 + group)
colnames(design) <- c("Normal", "Tumor")

fit <- lmFit(data_log, design)
cont_matrix <- makeContrasts(Tumor-Normal, levels=design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

# Extract full results table
results <- topTable(fit2, sort.by="P", n=Inf)

# 5. Visualizations
# A. PCA Plot to verify group clustering
pca_result <- prcomp(t(data_log))
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], Group = group)

ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="PCA: Tumor vs Normal Breast Tissue Cluster Analysis")

# B. Heatmap of Top 50 Genes
top_50_genes <- rownames(results)[1:50]
heatmap_data <- data_log[top_50_genes, ]
pheatmap(heatmap_data,
         annotation_col = data.frame(Group = group, row.names = colnames(data_log)),
         show_colnames = FALSE,
         scale = "row",
         main = "Heatmap: Top 50 Differentially Expressed Genes",
         color = colorRampPalette(c("blue", "white", "red"))(100))

