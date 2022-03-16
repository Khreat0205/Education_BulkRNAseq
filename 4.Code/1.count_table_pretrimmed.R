library(data.table)
library(biomaRt)

gctfh_files_path <- list.files("1.Data/raw_count_txt/GCTFH/",full.names = T)
gctfh_files <- list.files("1.Data/raw_count_txt/GCTFH/",full.names = F)
gct_sample <- gsub("^.*_D|.count.txt","",gctfh_files)
gctfh_info <- function(x){
  sample_celltype <- "GCTFH"
  sample_number <- substring(x,9,10)
  sample_day <- substring(x,1,2)
  result <- c(sample_celltype,sample_number,sample_day)

  return(result)

}

gct_group_table <- data.table(t(sapply(gct_sample,gctfh_info)),keep.rownames = T)
colnames(gct_group_table) <- c("sample","celltype","mouse","day")


tfh_files_path <- list.files("1.Data/raw_count_txt/TFHlike//",full.names=T)
tfh_files <- list.files("1.Data/raw_count_txt/TFHlike//")
tfh_sample <- gsub("^.*_D|.count.txt","",tfh_files)

tfh_info <- function(x){
  sample_celltype <- "TFHlike"
  sample_number <- substring(x,11,12)
  sample_day <- substring(x,1,2)
  result <- c(sample_celltype,sample_number,sample_day)
  return(result)
}

tfhlike_group_table <- data.table(t(sapply(tfh_sample,tfh_info)),keep.rownames = T)
colnames(tfhlike_group_table) <- c("sample","celltype","mouse","day")

group_table <- rbind(gct_group_table,tfhlike_group_table)
file_path <- c(gctfh_files_path,
                tfh_files_path)
group_table <- cbind(group_table,file_path)


for(i in 1:nrow(group_table)){

  counts <- fread(group_table$file_path[i],header=F)
  counts <- counts[grepl("ENS",counts$V1),]
  colnames(counts) <- c("ensembl_gene_id",group_table$sample[i])

  if(i ==1 ){
    count_table <- counts
  }else {
    count_table <- merge(count_table, counts, by="ensembl_gene_id",all=T)
  }
}

filtered_row <- rowSums(count_table[,-1] >= 10) > 1
filtered_count_table <- count_table[filtered_row,]



datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

G_list <- getBM(filters ="ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=filtered_count_table$ensembl_gene_id,mart= ensembl)

final_count_table <- merge(G_list,filtered_count_table,by="ensembl_gene_id",all=T)


fwrite(group_table,file = "1.Data/group_table.txt",sep="\t")
fwrite(final_count_table,file = "1.Data/count_table.txt",sep="\t")

