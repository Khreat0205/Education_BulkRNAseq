# Set working directory
# If you use R project file, you would not need this step.
# setwd("C:/Users/user/Documents/GitHub/Education_BulkRNAseq/")

# Which library to use
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)
library(ComplexHeatmap)

# Load data from previous script
dds <- readRDS(file="3.Results/deseqdatset.rds")
res <- readRDS(file="3.Results/result.rds")
sig_gene <- readRDS(file="3.Results/significant_gene.rds")

# Check if there is anything wrong with the data
head(dds)
head(res)
head(sig_gene)

# Create directory for saving plots
# dir.create('3.Results/Visualization/')



### Draw Volcano plot

# Volcano plot depicts genes up- or down-regulated with fold-change 2 and adjusted p value < 0.01 in
# CD90neg/lo GCTfh cells relative to CD90hi GCTfh-like cells at day 8, 12, 16 or day 24

# Let's draw volcano plot of day 8 first
# We will use the function named EnhancedVolcano in EnhancedVolcano package

day8.res <- res[[1]] # GCTFH cell vs TFH-like cell at day 8
head(day8.res)

EnhancedVolcano(

                ################  TO-DO ###################
                # Please fill parameters of this function to draw volcano plot!
                # Erase ### and fill proper codes below
                #
                # Explanation
                #    - toptable : data-frame of test statistics
                #    - lab : A column in toptable containing variable names
                #    - x : A column name in toptable containing log2 fold change
                #    - y : A column name in toptable containing p-value
                #
                # Tip
                #    - The $ operator can be used to select a column
                #    - You can use "colnames(x)" to see column names in data frame x.

                toptable = ###,
                lab = ### ,
                x = ### ,
                y = ### ,

                ################  TO-DO ###################


                pCutoff = 0.01, # Cut-off for statistical significance. (in raw value)
                FCcutoff = 2,   # Cut-off for absolute log2 fold-change
                title = 'Day 8', # Title of the plot
                subtitle = "Fold Change > 2, padj < 0.01", # Subtitle of the plot
                ylim = c(0, 15), # Limits of the y-axis
                pCutoffCol = 'padj',
                gridlines.minor = F,
                legendPosition = 'none')




# If this picture looks ok, let's use a for-loop to plot all 4 timepoints.
# You can save plot with ggsave function in ggplot2 packages.

exp.title <- c('Day8', 'Day12', 'Day16', 'Day24') # Vector for set titles in loop

for(i in 1:4){
  EnhancedVolcano(toptable = res[[i]],
                  lab = res[[i]]$mgi_symbol,
                  x ='log2FoldChange',
                  y = 'pvalue',
                  pCutoff = 0.01,
                  FCcutoff = 2,
                  title = exp.title[i],
                  subtitle = "Fold Change > 2, padj < 0.01",
                  ylim = c(0, 15),
                  pCutoffCol = 'padj',
                  gridlines.minor = F,
                  legendPosition = 'none')


  ggsave(filename = paste0('3.Results/Visualization/', exp.title[[i]], '_volcano.pdf'),
         device = 'pdf',
         width = 5,
         height= 8,
         units = 'in')
}




### Heatmap plot
# Heatmap graph showing the fold-change of Tfh-related gene expression in CD90neg/lo GCTfh
# over CD90hi GCTfh-like cells in each individual pLN.
# Numbers indicate the log2 fold-change in gene expression (n = 3 at each time point)

# Lets define the target gene in interest
Tfh.related.gene <- c('Padi4', 'Ascl2', 'S1pr2', 'Sostdc1',
                      'Asb2', 'Pou2af1', 'Ctla4', 'Fam43a', 'Pde3b',
                      'Pdcd1', 'Cebpb', 'Bcl6', 'P2rx4', 'Cd69', 'Icos',
                      'Id2','Cd27', 'Slc26a11', 'Cxcr5',
                      'Sh2d1a', 'Batf', 'Il4', 'Il6ra', 'Btla',
                      'Maf', 'Il21')

# Find index of target gene from the result file
idx <- match(Tfh.related.gene, res[[1]]$mgi_symbol)

# Count data can be extracted with the function counts() from dds object
countData <- data.frame(counts(dds, normalized=TRUE))

# Leave only the rows corresponding to the gene of interest.
Tfh.related.countData <- countData[idx, ]

# Set row names.
rownames(Tfh.related.countData) <- Tfh.related.gene

head(Tfh.related.countData)
dim(Tfh.related.countData)
colnames(Tfh.related.countData)


################  TO-DO ###################
#   Calculate log2 Fold change of each compare and assign as "log2FC.data"
#   Erase ### in square bracket and fill proper code below
#
#   Tips
#     - Calculation should be log2(TFH / TFHLike) for each compare.
#     - Check column names with "colnames(x)" function.
#     - Use "df[a:b, c:d]" to select slice of data frame.
#          ex1) df[1, ]  : Return row 1
#          ex2) df[, 5] : Return column 5
#          ex3) df[1:5, 2] : Return Rows 1:5 and column 2
#     - log2(x) computes the logarithm of the given column in base 2.
#           ex ) log2(10)

log2FC.data <- log2(Tfh.related.countData[ ### ] / Tfh.related.countData[ ### ])

###############  TO-DO ###################

head(log2FC.data)


# Set column names
exp.title <- c('Day8', 'Day12', 'Day16', 'Day24')
new.column.name <- paste0(rep(exp.title, each=3), '_', rep(c(1,2,3),4))
colnames(log2FC.data) <- new.column.name

# Set column split information. This will be used in next function
column_split_info <- factor(rep(exp.title, each=3), levels = exp.title)

# Draw heatmap plot with Heatmap function in ComplexHeatmap packages
h <- Heatmap(matrix = log2FC.data, # Target matrix for plotting
             name = 'log2FC', # Name of legend
             cluster_columns = F, # Whether to cluster column
             column_split = column_split_info, # How to split column
             cell_fun = function(j, i, x, y, width, height, fill) { # Function to plot values of each cell
               grid.text(sprintf("%.2f", log2FC.data[i, j]), x, y, gp = gpar(fontsize = 8))
             })
draw(h)


# You can save plot with these function
pdf('3.Results/Visualization/Heatmap.pdf', width = 6, height = 10)
draw(h)
dev.off()
