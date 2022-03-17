library(data.table)
library(biomaRt)
#### CD90neg/lo GCTfh cell 샘플들의 파일 경로 불러오기 및 파일명 기반 샘플 이름 지정
gctfh_files_path <- list.files("1.Data/raw_count_txt/GCTFH/",full.names = T)
gctfh_files <- list.files("1.Data/raw_count_txt/GCTFH/",full.names = F)
gct_sample <- gsub("^.*_D|.count.txt","",gctfh_files)
#### 샘플 이름에 따른 그룹 정보 추출
gctfh_info <- function(x){
  sample_celltype <- "GCTFH"
  sample_number <- substring(x,9,10)
  sample_day <- substring(x,1,2)
  result <- c(sample_celltype,sample_number,sample_day)

  return(result)

}
#### 추출된 그룹 정보를 테이블로 형성
sample_to_info <- sapply(gct_sample,gctfh_info)
head(sample_to_info)

gct_group_table <- data.table(t(sample_to_info),keep.rownames = T)
head(gct_group_table)

colnames(gct_group_table) <- c("sample","celltype","mouse","day")

#### CD90high GCTfh-like cell 샘플들의 파일 경로 불러오기 및 파일명 기반 샘플 이름 지정
tfh_files_path <- list.files("1.Data/raw_count_txt/TFHlike//",full.names=T)
tfh_files <- list.files("1.Data/raw_count_txt/TFHlike//")
tfh_sample <- gsub("^.*_D|.count.txt","",tfh_files)
#### 샘플 이름에 따른 그룹 정보 추출
tfh_info <- function(x){
  sample_celltype <- "TFHlike"
  sample_number <- substring(x,11,12)
  sample_day <- substring(x,1,2)
  result <- c(sample_celltype,sample_number,sample_day)
  return(result)
}
#### 추출된 그룹 정보를 테이블로 형성
tfhlike_group_table <- data.table(t(sapply(tfh_sample,tfh_info)),keep.rownames = T)
colnames(tfhlike_group_table) <- c("sample","celltype","mouse","day")

#### GCTfh cell, GCTfh-like cell 그룹 테이블 병합
group_table <- rbind(gct_group_table,tfhlike_group_table)
file_path <- c(gctfh_files_path,
               tfh_files_path)
#### GCTfh cell, GCTfh-like cell 그룹 정보 테이블에 파일 경로 삽입
group_table <- cbind(group_table,file_path)

#### group table에 삽입된 파일 경로를 하나씩 읽어 합치는 반복문
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


#### 최소 10 이상 count 된 유전자만 남김
filtered_row <- rowSums(count_table[,-1] >= 10) > 0
filtered_count_table <- count_table[filtered_row,]


#### 유전자 이름이 ensembl gene id로 되어 있어 MGI symbol로 변환
## ensembl
ensembl <- useEnsembl(biomart = "ensembl")
datasets <- listDatasets("ensembl")
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

G_list <- getBM(filters ="ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=filtered_count_table$ensembl_gene_id,mart= ensembl)
## 매치된 gene symbol 정보와 count table 병합
final_count_table <- merge(G_list,filtered_count_table,by="ensembl_gene_id",all=T)

## 파일 저장: 그룹 테이블
fwrite(group_table,file = "1.Data/group_table.txt",sep="\t")
## 파일 저장: count 테이블
fwrite(final_count_table,file = "1.Data/count_table.txt",sep="\t")

getwd()
