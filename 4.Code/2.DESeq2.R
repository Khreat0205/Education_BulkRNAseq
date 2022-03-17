# Set working directory
setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")

# Which library to use
library("DESeq2")

# Input : 1.count data
cts <- read.table("1.Data/count_table.txt", sep = "\t", header = T)
counts <- cts[,-c(1,2)]
featuredata <- cts[,c(1,2)]

# Input : 2.coldata
coldata <- read.table("1.Data/group_table.txt", sep = "\t", header = T)
## Match names
colnames(counts) <- coldata$sample
## Set as factor
coldata$celltype <- factor(coldata$celltype)
coldata$mouse <- factor(coldata$mouse)

# PCA Plot

# DEG analysis
## Divide data by day
counts_list <- list(day8 = subset(counts,select=grepl("^08", colnames(counts))),
                    day12 = subset(counts,select=grepl("^12", colnames(counts))),
                    day16 = subset(counts,select=grepl("^16", colnames(counts))),
                    day24 = subset(counts,select=grepl("^24", colnames(counts))))
coldata_list <- list(day8 = coldata[grepl("^08", coldata$sample),],
                     day12 = coldata[grepl("^12", coldata$sample),],
                     day16 = coldata[grepl("^16", coldata$sample),],
                     day24 = coldata[grepl("^24", coldata$sample),])

dds <- list()
res <- list()
sig_gene <- list()
for (i in 1:4){
  ## Make DESeq2 dataset
  dds[[i]] <- DESeqDataSetFromMatrix(countData = counts_list[[i]],
                                     colData = coldata_list[[i]],
                                     design = ~ mouse + celltype)
  mcols(dds[[i]]) <- data.frame(mcols(dds[[i]]), featuredata)
  ## Note on factor levels
  dds[[i]]$celltype <- relevel(dds[[i]]$celltype, ref="GCTFH")
  ## Run DESeq2
  dds[[i]] <- DESeq(dds[[i]])
  ## See results
  res[[i]] <- results(dds[[i]], saveCols=1:2)
  ## Significant gene (cutoff : FC>=2 % adjusted p value<0.01)
  sig_gene[[i]] <- res[[i]][which(abs(res[[i]]$log2FoldChange)>=1 & res[[i]]$padj<0.01),]
  }
