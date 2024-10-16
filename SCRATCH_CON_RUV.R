# The following represents code for standard downstream analyses of
# RNA-seq and ATAC-seq data, including determining factors of unwanted
# variation (RUV-seq) and various visualization techniques
#Load required libraries
library(DESeq2)
library(RUVSeq)
library(limma)
library(edgeR)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db) # For mouse, change to "org.Mm.eg.db"
library(data.table)
library(vsn)
library(edgeR)
library(dendextend)
library(matrixStats)
library(corrplot)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(plotly)
library(forcats)
library(cowplot)
library(dplyr)
library(GenomicFeatures)

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

# Remove low expressed genes
cts <- cts[apply(cts, 1, function(x) length(which(x>=5))/length(x) >0.25),] # at least 5 reads in at least 25% of samples
# or 
keep <- filterByExpr(cts, design=NULL)
cts <- cts[keep,]


#these check if they match and are in the same order
all(colnames(cts) %in% rownames(coldata))
all(colnames(cts) == rownames(coldata))

# Perform initial DESeq analysis (no correction for unwanted variation)
# Create a DESeq2 object. Set the design formula to include known covariates and comparison variable
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ disease)
dds$disease <- relevel(dds$disease, ref = "NIBD")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomLRT(dds, maxit = 1000, reduced = ~1)
uncorrected_results <- results(dds)


############## RUV control genes method 2 - no known contrast ################################
# Create model matrix to run limma removeBatchEffect
vsd <- vst(dds) 
mm <- model.matrix(~1, colData(vsd))
sex <- as.numeric(as.factor(coldata$sex))
batch <- as.numeric(as.factor(coldata$batch))

# Remove batch effect
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, covariates=cbind(sex), design=mm)
# Pull out the normalized matrix
nc <- as.data.frame(assay(vsd))
#row.names(nc) <- genes$gene

# Calculate the average gene expression of each genes and take the top 5000 highest expressed
nc$avg_exp <- rowMeans(nc,na.rm=TRUE)
nc <- nc %>% arrange(desc(avg_exp))
nc <- nc[1:5000,]
nc <- dplyr::select(nc, -avg_exp)

# Calculate the variance of each genes, and choose the lowest 1000 genes as the negative control gene
nc$cv <- rowSds(as.matrix(nc))/rowMeans(nc,na.rm=TRUE)
nc <- nc %>% arrange(cv)
nc <- nc[1:1000,]
nc <- nc %>% dplyr::select(-cv)

# Create a new Expression Set and perform an upper quartile normalization
nc <- round(nc*1000) #NEED TO FIND BETTER SCALING PARAMETER
ruv_matrix <- newSeqExpressionSet(as.matrix(nc), phenoData = data.frame(coldata, row.names = rownames(coldata)))
ruv_matrix <- betweenLaneNormalization(ruv_matrix, which="upper")
ctrlgene_names <- rownames(nc)
#####################################################################################


############################# Run RUVg #############################################
# re-run for various k as necessary or run in a loop

for (k in 1:20) {
  ruv<- RUVg(ruv_matrix, ctrlgene_names, k=k) # k is the number of factors
  plotRLE(ruv, outline=FALSE, main = paste("RLE plot for k=", k, sep=""), col=as.numeric(as.factor(coldata$disease)))
  #plotPCA(ruv, main = paste("PCA plot for k=", k, sep=""), col = as.numeric(as.factor(coldata$disease)))
}

#Add factors of variation to coldata
ruv_factors <- pData(ruv)[,-(1:ncol(coldata)-1)]
colData(dds) <- cbind(colData(dds), ruv_factors)
design(dds) <- ~disease + W_2 + W_3 + W_14 + W_16
dds$disease <- relevel(dds$disease, ref = "NIBD")
####################################################################################

############################## Run DESeq2 ###########################################
dds <- DESeq(dds)
results <- results(dds)

# Add gene symbol
ENSEMBL_org <- keys(org.Hs.eg.db, keytype = "ENSEMBL") # ENSEMBL keys don't include decimal places
results_ens.str <- substr(rownames(results), 1, 18)
results_ens.str <- gsub("\\..*","",results_ens.str) 
results$symbol <- mapIds(org.Hs.eg.db,
                         keys=results_ens.str,
                         column=c("SYMBOL"),
                         keytype="ENSEMBL",
                         multiVals="first")
resultsOrdered <- results[order(results$padj),]
write.table(resultsOrdered, file="DESeq_results.txt", sep="\t", quote=FALSE) 
####################################################################################

####################### Transform and Remove Batch Effects ##########################
# Transform data to use for visualizations 
vsd <- vst(dds) 
mm <- model.matrix(~ 1 + coldata$disease)

# Remove non-numeric columns from ruv_factors to fit in the assay thing
ruv_factors <- ruv_factors[, !colnames(ruv_factors) %in% c("montreal_perianal_disease")]
assay(vsd) <- removeBatchEffect(assay(vsd), batch=coldata$batch, 
                                covariates=cbind(as.numeric(as.factor(coldata$sex)), ruv_factors), 
                                design=mm)

######################################################################################


########################### PCA Plot #################################################
# Run PCA and write the output from the PCA analysis - the principle components
pcaData <- plotPCA(vsd, intgroup=c("disease"), ntop = 1000, returnData=TRUE)
write.table(pcaData, file = "PCA_data.txt", quote = FALSE, eol = "\n", sep = "\t", na = "NA")

# Create the PCA plot using ggplot that will be plotted to the screen for immediate review
pdf("PCA_ggplot.pdf");
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=disease, label="disease")) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()
#####################################################################################

############### RUV Factor Correlation Matrix #######################################
# Create correlation figure vs other phenotype data

# Make matrix of RUV Factors
fuv <- ruv_factors[, c("W_2", "W_3","W_14", "W_16"), drop = FALSE]
colnames(fuv) <- c("W_2", "W_3","W_14", "W_16")
fuv = data.frame(fuv)

##Make another matrix with covariates of interest (cbind), and correlate it with RUV factors (cor)
RUVMatrix = cor(cbind(as.numeric(as.factor(coldata$disease)),coldata$TIN_med_new),
                fuv, use = "complete.obs", method = "pearson")

rownames(RUVMatrix) = #add names to rows based on your covariates of interest
  c("Disease", "TIN_med")

#Make correlation plot (corrplot R package)
corrplot(RUVMatrix)
###########################################################################################

################### Volcano Plots #########################################################
tab = data.frame(logFC = results$log2FoldChange, negLogPval = -log10(results$padj))
head(tab)
par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~pvalue))
lfc = 2 # Draw horizontal line at logFC = 2
pval = 0.01 # Draw vertical line at p-value = 0.01
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
points(tab[signGenes, ], pch = 16, cex = 0.8, col = "red") # Color significant genes red
abline(h = -log10(pval), col = "green3", lty = 2) 
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
#############################################################################################



