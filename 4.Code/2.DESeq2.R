library("DESeq2")
# Input : cts <- count data, coldata <- coldata

# PCA Plot
SummarizedExperiment(assays = cts, colData = coldata)


# DEG analysis
## Make DESeq Dataset
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~)

## Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
sig_gene <- res[which(abs(res$log2FoldChange)>=1) & res$padj<0.01),]