library(DESeq2) 
library(limma)
library(sva)
library(ggrepel)
library(ggplot2) 
library(plotly)
library(plyr)
library(dplyr)
library(GGally)
library(tidyverse)
library(scales)
library(pheatmap)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(edgeR)
library(preprocessCore)
library(RColorBrewer) 
library(M3C)
library(gplots)
library(VennDiagram)
library(GeneOverlap)
library(EnhancedVolcano)
library(ggpubr)
library("AnnotationDbi")
library("org.Hs.eg.db")

# here i am setting my directory
setwd("/Users/ayeshawad/Desktop/UNC/Orientation")

#here i am reading my files
coldata <- read.csv("RNA-seq_coldata_human_adult_2023-08-04.csv", check.names = "FALSE")
cts <- read.csv("CD_vs_NIBD_colon_uninflamed_raw_counts.csv",row.names = 1, check.names = "FALSE")

#these are my filters for whatever i am looking for (always change). Also used filter and Base R to show how to do both
coldata <- coldata %>%
  filter(tissue == "colon" & inflammation == "uninflamed")
coldata <- coldata[!is.na(coldata$inflammation) & coldata$inflammation != "", ]

#this is just cleaning the data that we read wrongly. no need to do this all the time
coldata[, 1] <- make.unique(as.character(coldata[, 1]))
rownames(coldata) <- coldata$sample
coldata$sample <- NULL  # Remove the "gene_id" column from the dataframe

# here i will filter the samples that are in both the coldata and counts files
coldata <- coldata[rownames(coldata) %in% colnames(cts), ]

# here i will sort the rows in coldata and columns in counts to be in same order
coldata = coldata[ with(coldata, order(rownames(coldata))),]
cts = cts[, with(cts, order(colnames(cts)))]



#these check if they match and are in the same order
all(colnames(cts) %in% rownames(coldata))
all(colnames(cts) == rownames(coldata))

# now we need to construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                       colData = coldata,
                       design = ~ disease)

#we can filter and do a lot of things here, but we can also jump straight into the analysis

dds = DESeq(dds)

resultsCell1 <- results(dds, contrast = c("disease","CD","NIBD"))
resultsCell1ord <- resultsCell1[order(resultsCell1$padj),]
resultsCell1 <- as.data.frame(resultsCell1ord)[, ]

# Get the gene symbol for the EnsEMBL id
# Remove version numbers correctly and then ensure no leading or trailing spaces
ens.str <- gsub("\\..*", "", rownames(resultsCell1))  # Remove the dot and everything after it

# Print to check if the transformation worked
print("Cleaned Ensembl IDs:")
print(head(ens.str))

# Check valid Ensembl keys from the org.Hs.eg.db database
valid_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")

# Print a few valid Ensembl IDs for comparison
print("Valid Ensembl IDs from the database:")
print(head(valid_keys))

# Map the symbols using the cleaned Ensembl IDs
resultsCell1$symbol <- mapIds(org.Hs.eg.db,
                              keys = ens.str,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")
resultsCell1$entrez <- mapIds(org.Hs.eg.db,
                              keys=ens.str,
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")

# Sort by the adjusted p-value
resultsCell1ord <- resultsCell1[order(resultsCell1$padj),]
resultsCell1ordDF <- as.data.frame(resultsCell1ord)[, ]

# Write out the results to a file - may want to include full path
write.table(resultsCell1ordDF, file = "DESEQ_Cell1_vs_control_20230713.txt", quote = FALSE, eol = "\n", sep = "\t", na = "NA")

# Create a volcano plot
# Create the volcano plot as before
tab = data.frame(logFC = resultsCell1$log2FoldChange, negLogPval = -log10(resultsCell1$pvalue), symbol = resultsCell1$symbol)
head(tab)
par(mar = c(5, 4, 4, 4))
plot(tab$logFC, tab$negLogPval, pch = 16, cex = 0.6, 
     xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))

# Draw significance lines
lfc = 2
pval = 0.01
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab$logFC[signGenes], tab$negLogPval[signGenes], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)

# Identify the top 10 most differentially expressed genes (top 5 up and 5 down based on log2FC and padj)
top_upregulated <- head(tab[order(-tab$logFC, tab$negLogPval), ], 5)
top_downregulated <- head(tab[order(tab$logFC, tab$negLogPval), ], 5)

# Combine the top upregulated and downregulated genes
top_genes <- rbind(top_upregulated, top_downregulated)

# Add labels to the most differentially expressed genes using `text()` with better positioning and larger size
if (nrow(top_genes) > 0) {
  # Use `pos = 4` for labels to the right, and add a slight offset using `adj`
  text(top_genes$logFC, top_genes$negLogPval, labels = top_genes$symbol, 
       pos = 4, cex = 1, col = "blue", font = 2, adj = 0.3)
} else {
  print("No top genes identified for labeling.")
}
# Next comparison can go here....

# Save the R object that be directly loaded in the future, if you want (I don't usually do this)
# saveRDS(dds, file = "All_samples_20230713.rds")
# Read back in the object
# dds <- readRDS("All_samples_20230713.rds")


# For PCA plots, want to perform a batch correction first using limmaBatchCorrect
# First do a variance stabilizing transformation
vsd = vst(dds)

# Extract the matrix of normalized expression values
mat= assay(vsd)

#checking something 
# Perform PCA using prcomp on the transposed matrix
pca_res <- prcomp(t(mat))  # Transpose because prcomp expects samples as rows

# Extract the PCA loadings (eigenvectors)
pca_loadings <- pca_res$rotation  # Loadings matrix

# Display the top contributing genes for PC1 and PC2
top_pc1_genes <- head(sort(abs(pca_loadings[, "PC1"]), decreasing = TRUE), 10)
top_pc2_genes <- head(sort(abs(pca_loadings[, "PC2"]), decreasing = TRUE), 10)

top_genes_df <- data.frame(
  Gene = rownames(pca_loadings),
  PC1_Loading = pca_loadings[, "PC1"],
  PC2_Loading = pca_loadings[, "PC2"]
)

# Filter for only the top 10 genes in each category
top_genes_df <- top_genes_df[rownames(top_genes_df) %in% c(names(top_pc1_genes), names(top_pc2_genes)), ]

# Convert Ensembl IDs to gene symbols using mapIds()
top_genes_df$Symbol <- mapIds(org.Hs.eg.db,
                              keys = gsub("\\..*", "", top_genes_df$Gene),  # Remove version numbers
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")

# Replace NA symbols with the original Ensembl IDs (fallback)
top_genes_df$Symbol[is.na(top_genes_df$Symbol)] <- top_genes_df$Gene[is.na(top_genes_df$Symbol)]

# Print the dataframe with symbols for verification
print("Top Contributing Genes for PC1 and PC2 with Symbols:")
print(top_genes_df)

# Plot top loadings for PC1 using gene symbols
top_pc1_genes_df <- top_genes_df[order(abs(top_genes_df$PC1_Loading), decreasing = TRUE), ]
barplot(top_pc1_genes_df$PC1_Loading, names.arg = top_pc1_genes_df$Symbol, las = 2,
        main = "Top Genes Contributing to PC1", col = "skyblue", cex.names = 0.8)

# Plot top loadings for PC2 using gene symbols
top_pc2_genes_df <- top_genes_df[order(abs(top_genes_df$PC2_Loading), decreasing = TRUE), ]
barplot(top_pc2_genes_df$PC2_Loading, names.arg = top_pc2_genes_df$Symbol, las = 2,
        main = "Top Genes Contributing to PC2", col = "lightgreen", cex.names = 0.8)
# Perform the batch/covariate removal, if necessary
#mat <- removeBatchEffect(mat, batch = batch, covariate = cbind())

# Store the transformed and corrected matrix
assay(vsd) = mat

# Write out the matrix, if desired
write.table(mat, file = "All_samples_corrected_20230713.txt", quote = FALSE, eol = "\n", sep = "\t", na = "NA")

# Can directly create the PCA plot and save to a file
# This PCA format is very basic
pdf("PCA_all_samples.pdf");
# Color samples based on "cell_type" in coldata file, use the top 1000 most variable genes
plotPCA(vsd, intgroup = c("disease"), ntop = 1000)
dev.off()

# Alternatively, can use ggplot to create PCA plot to give more flexibility in formatting
# Need to extract sample information for the plot
disease = coldata$disease;
sample = coldata$sample;

# Capture PCA data to increase how data is displayed, including putting on labels of points
# ntop - specifies how many genes used in the PCA analysis
pcaData <- plotPCA(vsd, intgroup=c("disease"), ntop = 1000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
head(pcaData)
print(pcaData)
print(percentVar)
head(assay(vsd))

# Write the output from the PCA analysis - the principle components - if desired 
write.table(pcaData, file = "PCA_data_all_samples_20230713.txt", quote = FALSE, eol = "\n", sep = "\t", na = "NA")

# Create the PCA plot using ggplot that will be plotted to the screen for immediate review
pdf("PCA_all_samples_ggplot.pdf");
ggplot(pcaData, aes(PC1, PC2, color=disease, label=disease)) +
  geom_point(size=3) +
  geom_text(aes(label=sample),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

#INDEPENDENT ANALYSIS 
# Merge PCA data with coldata (ensure `coldata` has a `sample` column)
pcaData$sample <- rownames(pcaData)
coldata$sample <- rownames(coldata)  # Create a sample column in coldata
merged_data <- merge(pcaData, coldata, by = "sample")

# List of columns to focus on
columns_of_interest <- c("group", "sex", "subtype", "IECAve", "CD3Ave", "BCellAve", "NeutrophilAve", 
                         "NaturalKillerAve", "FibroAve", "ageSampleCollection")

# Select only the columns of interest from merged_data
subset_data <- merged_data[, c("PC1", "PC2", columns_of_interest), drop = FALSE]

# Perform correlation analysis for specified columns
# Ensure numeric conversion where necessary (e.g., converting categorical variables to numeric factors)
correlation_results_pc1 <- sapply(subset_data[, columns_of_interest], function(x) cor(subset_data$PC1, as.numeric(factor(x)), use = "complete.obs"))
correlation_results_pc2 <- sapply(subset_data[, columns_of_interest], function(x) cor(subset_data$PC2, as.numeric(factor(x)), use = "complete.obs"))

# Combine the correlation results into a single dataframe for easy plotting
correlation_df <- data.frame(
  Factor = columns_of_interest,
  PC1_Correlation = correlation_results_pc1,
  PC2_Correlation = correlation_results_pc2
)

# Create a correlation matrix for heatmap
cor_matrix <- as.matrix(correlation_df[, -1])  # Remove the 'Factor' column

# Set row names to metadata factors
rownames(cor_matrix) <- correlation_df$Factor

# Load the pheatmap library if not already loaded
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# Create the heatmap for PC1 and PC2 correlations, with clear distinction
pheatmap(cor_matrix,
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Correlation Heatmap of Metadata Factors with PC1 and PC2",
         legend_breaks = c(-1, -0.5, 0, 0.5, 1),
         fontsize_row = 12)

