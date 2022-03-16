setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")
library("DESeq2")
# Input : cts <- count data, coldata <- coldata
counts <- read.table("1.Data/count_table.txt", sep = "\t", header = T)
coldata <- read.table("1.Data/group_table.txt", sep = "\t", header = T)

# PCA Plot
# SummarizedExperiment(assays = counts, colData = coldata)

# DEG analysis
## Divide by day

## Make DESeq Dataset

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~)

## Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
sig_gene <- res[which(abs(res$log2FoldChange)>=1) & res$padj<0.01),]