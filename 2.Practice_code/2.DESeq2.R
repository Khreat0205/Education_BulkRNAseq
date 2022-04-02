# Set working directory
# If you use R project file, you would not need this step.
setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")

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


# PCA Plot
#<<<<<<< HEAD
## Sunah
#=======
## normalize data & remove batch effects
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design = ~ day+celltype+day:celltype)
dds <- DESeq(dds)

saveRDS(dds, file="3.Results/deseqdatset.rds")


## draw PCA plot
pcaData <- plotPCA(          , intgroup=c('day','celltype'), returnData=TRUE)
pcaPercentVar <- round(100 * attr(pcaData, "percentVar"))
pcaPlot <- ggplot(        ) +
  geom_point(mapping = aes(PC1, PC2, color=day, shape=celltype), size=3) +
  xlab(paste0("PC1: ",pcaPercentVar[1],"% variance")) +
  ylab(paste0("PC2: ",pcaPercentVar[2],"% variance")) +
  coord_fixed()+
  theme_classic()
pcaPlot

ggsave(        , filename = "./3.Results/Visualization/PCA.pdf", dpi = 600)

# Correlation plot
## correlation
corData <- cor(       , method = 'pearson')
## draw plot
corPlot <- Heatmap(matrix =       ,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if(!is.na(corData[i,j]))
                       grid.text(round(corData[i, j],digits = 2), x, y, gp = gpar(fontsize = 10))
                   },
                   cluster_rows = T, cluster_columns = T,
                   heatmap_legend_param = list(title="Correlation"),
                   rect_gp = gpar(col = 'white', lwd = 3))
png(filename = './3.Results/Visualization/Correlation.png', width = 900, height = 600)
draw(       )
dev.off()

#>>>>>>> ec698ff07c187fa7edbfc9ecd0d4ff263211d65a

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

## Run for loop
for (i in 1:4){
  ################  TO-DO ###################
  # Make DESeq2 dataset
  # Fill in which input to use as count data, and column data (metadata)
  
  dds[[i]] <- DESeqDataSetFromMatrix(countData =           [[i]],
                                     colData =             [[i]],
                                     ## correcting for effect of mouse + calculating DEG between celltype
                                     design = ~ mouse + celltype)
  
  ################  TO-DO ###################
  
  ## Add metadata column (in this case, gene name)
  mcols(dds[[i]]) <- data.frame(mcols(dds[[i]]), featuredata)
  ## Note on factor levels (which condition to use as reference) (in this case, FC=GCTFH/TFHlike)
  dds[[i]]$celltype <- relevel(dds[[i]]$celltype, ref="TFHlike")
  ## Run DESeq2
  dds[[i]] <- DESeq(dds[[i]])
  ## See results
  res[[i]] <- results(dds[[i]], saveCols=1:2)
  ## Significant gene (cutoff : FC (Fold Change) >=2 % adjusted p value<0.01)
  sig_gene[[i]] <- res[[i]][which(abs(res[[i]]$log2FoldChange)>=1 & res[[i]]$padj<0.01),]
  ## Order by absolute value of log2FC (high to low)
  sig_gene_byFC[[i]] <- sig_gene[[i]][order(abs(sig_gene[[i]]$log2FoldChange),decreasing=T),]
  
  ################  TO-DO ###################
  # As we have learned how to use "order" function by sorting by absolute value of log2FC,
  # Let's sort significant genes by adjusted p value and assign as "sig_gene_bypadj"
  # Tips
  #   - order() returns row indices in a certain sorted manner
  #   - adjusted p value is saved in "padj" column
  #   - since the default paramter of order() is to sort in ascending order, you don't have to use any paramter in this situation
  
  sig_gene_bypadj[[i]] <- sig_gene[[i]][order(                  ),]
  
  ###############  TO-DO ###################
}


saveRDS(res, file="3.Results/result.rds")
saveRDS(sig_gene, file="3.Results/significant_gene.rds")
