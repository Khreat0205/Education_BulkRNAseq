# Set working directory
# If you use R project file, you would not need this step.
# setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")

# Which library to use
library("DESeq2")
library("ggplot2")
library("ComplexHeatmap")

# Input : 1.count data
cts <- read.table("1.Data/count_table.txt", sep = "\t", header = T)
## Raw count data per sample
counts <- cts[,-c(1,2)]
## Ensembl & Gene symbol information
featuredata <- cts[,c(1,2)]

# Input : 2.coldata
coldata <- read.table("1.Data/group_table.txt", sep = "\t", header = T)
## Match names
colnames(counts) <- coldata$sample
## Set as factor
coldata$celltype <- factor(coldata$celltype)
coldata$mouse <- factor(coldata$mouse)
coldata$day <- factor(coldata$day)

# View data
head(counts)
head(coldata)

# PCA Plot
## normalize data & remove batch effects
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design = ~ day+celltype+day:celltype)
dds <- DESeq(dds)

saveRDS(dds, file="3.Results/deseqdatset.rds")


## draw PCA plot
pcaData <- plotPCA(vst(dds, blind = F), intgroup=c('day','celltype'), returnData=TRUE)
pcaPercentVar <- round(100 * attr(pcaData, "percentVar"))
pcaPlot <- ggplot(pcaData) +
            geom_point(mapping = aes(PC1, PC2, color=day, shape=celltype), size=3) +
            xlab(paste0("PC1: ",pcaPercentVar[1],"% variance")) +
            ylab(paste0("PC2: ",pcaPercentVar[2],"% variance")) +
            coord_fixed()+
            theme_classic()
pcaPlot

ggsave(pcaPlot, filename = "./3.Results/Visualization/PCA.pdf", dpi = 600)

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
png(filename = './3.Results/Visualization/Correlation.png', width = 900, height = 600)
draw(corPlot)
dev.off()


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

## Before for loop, assign empty list
dds <- list()
res <- list()
sig_gene <- list()
sig_gene_bypadj <- list()
sig_gene_byFC <- list()

## First, Let's see Day 8
### Make DESeq2 dataset
dds[[1]] <- DESeqDataSetFromMatrix(countData = counts_list[[1]],
                                   colData = coldata_list[[1]],
                                   ## correcting for effect of mouse + DEG between celltype
                                   design = ~ mouse + celltype)
### Add metadata column (in this case, gene name)
mcols(dds[[1]]) <- data.frame(mcols(dds[[1]]), featuredata)
### Note on factor levels (which condition to use as reference) (in this case, FC=GCTFH/TFHlike)
dds[[1]]$celltype <- relevel(dds[[1]]$celltype, ref="TFHlike")
### Run DESeq2
dds[[1]] <- DESeq(dds[[1]])
### See results
res[[1]] <- results(dds[[1]], saveCols=1:2)
### Significant gene (cutoff : FC>=2 % adjusted p value<0.01)
sig_gene[[1]] <- res[[1]][which(abs(res[[1]]$log2FoldChange)>=1 & res[[1]]$padj<0.01),]
### Order by absolute value of log2FC (high to low)
sig_gene_byFC[[1]] <- sig_gene[[1]][order(abs(sig_gene[[1]]$log2FoldChange),decreasing=T),]
### Order by adjusted p value (low to high)
sig_gene_bypadj[[1]] <- sig_gene[[1]][order(sig_gene[[1]]$padj),]

## Run for loop
for (i in 2:4){
  dds[[i]] <- DESeqDataSetFromMatrix(countData = counts_list[[i]],
                                     colData = coldata_list[[i]],
                                     design = ~ mouse + celltype)
  mcols(dds[[i]]) <- data.frame(mcols(dds[[i]]), featuredata)
  dds[[i]]$celltype <- relevel(dds[[i]]$celltype, ref="TFHlike")
  dds[[i]] <- DESeq(dds[[i]])
  res[[i]] <- results(dds[[i]], saveCols=1:2)
  sig_gene[[i]] <- res[[i]][which(abs(res[[i]]$log2FoldChange)>=1 & res[[i]]$padj<0.01),]
  sig_gene_byFC[[i]] <- sig_gene[[i]][order(abs(sig_gene[[i]]$log2FoldChange),decreasing=T),]
  sig_gene_bypadj[[i]] <- sig_gene[[i]][order(sig_gene[[i]]$padj),]
}


saveRDS(res, file="3.Results/result.rds")
saveRDS(sig_gene, file="3.Results/significant_gene.rds")
