### DESeq2 @ CSHL
### Using Travor's data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158493
### Ref paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10263444/

# Install the latest version of DEseq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

# Me experimenting with some extra packages you can ignore 
#BiocManager::install("apeglm")
#BiocManager::install("mvtnorm")
#library("apeglm")
#library("mvtnorm")

# load the DESeq2 library
library("DESeq2")

#set the working directory
setwd("/Users/irisvanterve/Desktop/CSHL_Computational_Genomics_course/Project")

# Read in the raw read counts, I remove the first row as it contains overhead
rawCounts <- read.delim("GSE158491_bulkCD4_counts.txt",sep="\t", header=TRUE, skip=1, row.names = "Geneid")
View(rawCounts)

# Here I'm keeping all columns, but remove column 2-6 as it contains data we don't need in the count matrix 
rawCounts <- rawCounts[,-(1:5)]
View(rawCounts)

# Read in the sample information
sampleData <- read.delim("bulkrnacol.txt",sep="\t", header=TRUE)
View(sampleData)

# sampleData$condition needs to be in factor format
str(sampleData) # check format
sampleData$condition <- factor(sampleData$condition) #change format
str(sampleData) # check format again

# Generate the datamatrix
dds <- DESeqDataSetFromMatrix(countData = rawCounts, colData = sampleData, design= ~ condition)
dds <- DESeq(dds)

# output the names of the groups that we compare
resultsNames(dds) #lists the coefficients
res <- results(dds, name="condition_Fetal_Spleen_naiveCD4_vs_Adult_PBMC_naiveCD4") #needs to be one group that comes out of the line above

# Or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="condition_Fetal_Spleen_naiveCD4_vs_Adult_PBMC_naiveCD4", type="apeglm")

# Make variance stabilizing transformation before PCA&heatmap
vsd <- vst(dds, blind=FALSE) # FALSE means transformation will consider treatment groups

## Principle Component Analysis
library("ggplot2")
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

## Heatmap
library(pheatmap)
n <- 500  # Select the top n most variable genes
select <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:n]
data <- assay(vsd)[select,]
data <- t(scale(t(data))) #scaling the rows to make the gene expression comparable across genes

# Create annotation for heatmap
annotations <- data.frame(Condition = dds$condition)
rownames(annotations) <- colnames(dds)

# Change colors to make it similar to paper
condition_colors <- c("Adult" = "yellow3", "Fetal" = "darkgreen", "Newborn"="lightblue")

# Rename conditions to make it shorter and similar to paper
annotations$Condition <- factor(annotations$Condition, 
                                levels = c("Fetal_Spleen_naiveCD4", "Newborn_UmbCordBld_naiveCD4","Adult_PBMC_naiveCD4"), 
                                labels = c("Fetal","Newborn", "Adult"))

# Make the actual heatmap
pheatmap(data, 
        cluster_cols = TRUE, 
         cluster_rows = TRUE,
         scale = 'none', # I scaled the data already above
         show_rownames = FALSE, 
         show_colnames = FALSE, 
        main = "Heatmap of Differentially Expressed Genes", #title of the plot
        annotation_col = annotations,
        color = colorRampPalette(c("red", "white", "blue"))(50), #nr specifies the number of colors to be used in the color palette
        annotation_colors = list(Condition = condition_colors))
         

## Volcano plot
# chose two groups to compare
res <- results(dds, contrast=c("condition", "Fetal_Spleen_naiveCD4", "Adult_PBMC_naiveCD4"))

# data adjustment to prevent NAs
res$padj[is.na(res$padj)] <- 1  # Replace NA with 1 (as -log10(1) = 0)
res$padj[is.infinite(res$padj)] <- 1  # Replace Inf with 1

# make sure res is a dataframe
res_df <- as.data.frame(res)

# Set thresholds
threshold_pvalue <- 0.05
threshold_log2FoldChange <- 5 

# Create a new column to identify genes to label > we need the actual gene names here which we now don't have
res_df$label <- ifelse(res_df$pvalue < threshold_pvalue & abs(res_df$log2FoldChange) > threshold_log2FoldChange, rownames(res_df), "")

# Create the volcano plot (+example on how to save the plot)
volcano <- volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point() +
  #geom_text(aes(label = label), vjust = 1.5, hjust = 0.5, check_overlap = TRUE, angle = 45) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", title = "Volcano Plot: Fetal vs Adult") +
  geom_hline(yintercept = -log10(threshold_pvalue), color = "grey", linetype = "dashed") +
  geom_vline(xintercept = c(-threshold_log2FoldChange, threshold_log2FoldChange), color = "grey", linetype = "dashed")
print(volcano_plot)

ggsave("Volcano Plot: Fetal vs Adult.png", plot = volcano, bg="white")
