# Set working directory
# If you use R project file, you would not need this step.
# setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")


# Which library to use
library("DESeq2")
library("ggplot2")
library("ComplexHeatmap")


# Input : 1.count data
cts <- read.table("1.Data/count_table.txt", sep = "\t", header = T)

head(cts)

## Raw count data per sample
# Remove 1st, 2nd column of cts and assign it as counts
counts <- cts[,-c(1,2)]

head(counts)

## Ensembl & Gene symbol information
# Extract 1st, 2nd column of cts and assign it as featuredata
featuredata <- cts[,c(1,2)]

head(featuredata)

# Input : 2.coldata
coldata <- read.table("1.Data/group_table.txt", sep = "\t", header = T)

head(coldata)

## Match names
colnames(counts) <- coldata$sample

## Set as factor
# factor() converts any value into categorical data
coldata$celltype <- factor(coldata$celltype)
coldata$mouse <- factor(coldata$mouse)
coldata$day <- factor(coldata$day)

coldata$celltype


# PCA Plot
## Normalization & Batch correction
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = coldata,
                              # correcting for effect of day (first value) + DEG between celltype (last value)
                              design = ~ day+celltype)
dds <- DESeq(dds)

saveRDS(dds, file="3.Results/deseqdatset.rds")


## Draw PCA plot
pcaData <- plotPCA(vst(dds, blind = F), # vst() returns tranformed, normalized counts from DESeqDataset
                   intgroup=c('day','celltype'), returnData=TRUE)

pcaData

pcaPercentVar <- round(100 * attr(pcaData, "percentVar")) # from pcaData, extract "percentVar" (percent of variance)
pcaPlot <- ggplot(pcaData) +
            # draw dot plot (as axis x:PC1, y:PC2)
            geom_point(mapping = aes(PC1, PC2, color=day, shape=celltype), size=3) + 
            # write axis label
            xlab(paste0("PC1: ",pcaPercentVar[1],"% variance")) +
            ylab(paste0("PC2: ",pcaPercentVar[2],"% variance")) +
            # draw plot with a ratio of 1:1
            coord_fixed() +
            # draw plot with classic theme
            theme_classic()

pcaPlot

ggsave(pcaPlot, filename = "./3.Results/Visualization/PCA.pdf", dpi = 600)


# Correlation plot
## Calculate correlation
corData <- cor(counts, method = 'pearson') # we used pearson's correlation, but you can use kendall rank correlation and spearman's correlation as well

## Draw correlation plot
corPlot <- Heatmap(matrix = corData,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     # if value in corData is not <NA>
                     if(!is.na(corData[i,j]))
                       # write correlation value in each heatmap body
                       grid.text(round(corData[i, j],digits = 2), x, y, gp = gpar(fontsize = 10))
                   },
                   # do clustering and show it in row and column
                   cluster_rows = T, cluster_columns = T,
                   # legend title as "Correlation"
                   heatmap_legend_param = list(title="Correlation"),
                   # draw white rectangular line around heatmap body
                   rect_gp = gpar(col = 'white', lwd = 3))

## Save as png file
png(filename = './3.Results/Visualization/Correlation.png', width = 900, height = 600) 
draw(corPlot)
dev.off()


# DEG analysis
## Divide data by day
counts_list <- list(# subset() selects columns of which matched pattern ("^08" - by grepl()) is in column names of counts
                    day8 = subset(counts,select=grepl("^08", colnames(counts))), 
                    day12 = subset(counts,select=grepl("^12", colnames(counts))),
                    day16 = subset(counts,select=grepl("^16", colnames(counts))),
                    day24 = subset(counts,select=grepl("^24", colnames(counts))))
coldata_list <- list(# find indices of matched pattern in "sample" column of coldata & leave only selected rows
                     day8 = coldata[grepl("^08", coldata$sample),],
                     day12 = coldata[grepl("^12", coldata$sample),],
                     day16 = coldata[grepl("^16", coldata$sample),],
                     day24 = coldata[grepl("^24", coldata$sample),])

head(counts_list)
head(coldata_list)

## Before for loop, assign empty list
dds <- list()
res <- list()
sig_gene <- list()
sig_gene_bypadj <- list()
sig_gene_byFC <- list()

## First, Let's see Day 8
### Make DESeq2 dataset
dds[[1]] <- DESeqDataSetFromMatrix(# To assess the first element in the list, use list[[1]]
                                   countData = counts_list[[1]],
                                   colData = coldata_list[[1]],
                                   # correcting for effect of mouse (first value) + DEG between celltype (last value)
                                   design = ~ mouse + celltype)

head(counts_list[[1]])
head(coldata_list[[1]])

### Add metadata column (in this case, gene name)
# mcols() return metadata columns, so add featuredata (Gene information) beside them, and renew them
mcols(dds[[1]]) <- data.frame(mcols(dds[[1]]), featuredata)

### Note on factor levels (which condition to use as reference) (in this case, FC=GCTFH/TFHlike)
dds[[1]]$celltype <- relevel(dds[[1]]$celltype, ref="TFHlike")

### Run DESeq2 (Normalization & Batch correction)
dds[[1]] <- DESeq(dds[[1]])

### See results
# results() returns log2FoldChange, pvalue, etc. from DESeq2Dataset
# saveCols=1:2 (parameter) means to save 1st, 2nd column of metadata beside the extracted results
res[[1]] <- results(dds[[1]], saveCols=1:2)

head(res[[1]])

### Significant gene (cutoff : FC>=2 & adjusted p value<0.01)
# which() returns indices which fit the written condition
# (in this case, values of log2FoldChange column has to be more than 1 and values of padj column has to be less than 0.01)
# abs() returns absolute value
sig_gene[[1]] <- res[[1]][which(abs(res[[1]]$log2FoldChange)>=1 & res[[1]]$padj<0.01),]

head(sig_gene[[1]])

### Order by absolute value of log2FC (high to low)
# order() returns row indices in a certain sorted manner
# decreasing=T (parameter) means to sort in a descending manner
sig_gene_byFC[[1]] <- sig_gene[[1]][order(abs(sig_gene[[1]]$log2FoldChange),decreasing=T),]

head(sig_gene_byFC[[1]])

### Order by adjusted p value (low to high)
# since the default paramter is to sort in a ascending manner, you don't have to write any other parameter
sig_gene_bypadj[[1]] <- sig_gene[[1]][order(sig_gene[[1]]$padj),]

head(sig_gene_bypadj[[1]])


## Run for loop
# i is assigned as 2, 3, 4 in each turn
# so write i in the place where the value changes in each turn
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


## Save RDS file
saveRDS(res, file="3.Results/result.rds") # saveRDS(object, file="filename.rds")
saveRDS(sig_gene, file="3.Results/significant_gene.rds")
