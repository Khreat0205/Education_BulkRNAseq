# Set working directory
# If you use R project file, you would not need this step.
setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")

# Which library to use
library("DESeq2")
library("ggplot2")
library("ComplexHeatmap")

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
coldata$day <- factor(coldata$day)


# PCA Plot
## normalize data & remove batch effects
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ day+celltype+day:celltype)
dds <- DESeq(dds)
## draw PCA plot
pcaData <- plotPCA(vst(dds, blind = F), intgroup=c('day','celltype'), returnData=TRUE)
pcaPercentVar <- round(100 * attr(pcaData, "percentVar"))
pcaPlot <- ggplot(pcaData, aes(PC1, PC2, color=day, shape=celltype)) +
            geom_point(size=3) +
            xlab(paste0("PC1: ",pcaPercentVar[1],"% variance")) +
            ylab(paste0("PC2: ",pcaPercentVar[2],"% variance")) +
            coord_fixed()+
            theme_classic()
pcaPlot

# Correlation plot
## correlation
corData <- cor(counts, method = 'pearson')
## draw plot
corPlot <- Heatmap(matrix = corData,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if(!is.na(corData[i,j]))
                       grid.text(round(corData[i, j],digits = 2), x, y, gp = gpar(fontsize = 10))
                   },
                   cluster_rows = T, cluster_columns = T,
                   heatmap_legend_param = list(title="Correlation"),
                   rect_gp = gpar(col = 'white', lwd = 3))
draw(corPlot)

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
  dds[[i]]$celltype <- relevel(dds[[i]]$celltype, ref="TFHlike")
  ## Run DESeq2
  dds[[i]] <- DESeq(dds[[i]])
  ## See results
  res[[i]] <- results(dds[[i]], saveCols=1:2)
  ## Significant gene (cutoff : FC>=2 % adjusted p value<0.01)
  sig_gene[[i]] <- res[[i]][which(abs(res[[i]]$log2FoldChange)>=1 & res[[i]]$padj<0.01),]
  ## Order by adjusted p value
  sig_gene[[i]] <- sig_gene[[i]][order(sig_gene[[i]]$padj),]
}

