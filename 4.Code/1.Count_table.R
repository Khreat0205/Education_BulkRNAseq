library(data.table)
library(biomaRt)
getwd()
#####################################################################################
#### CD90neg/lo GCTfh cell 샘플들의 파일 경로 불러오기 및 파일명 기반 샘플 이름 지정
gctfh_files_path <- list.files("1.Data/raw_count_txt/GCTFH/",full.names = T)
head(gctfh_files_path)

gctfh_files <- list.files("1.Data/raw_count_txt/GCTFH/",full.names = F)
head(gctfh_files)

# ex) GSM4413011_D08GCTFH-01.count.txt
# ^ = 개행 문자, 문자열의 시작점
# . = 모든 문자
# * = 앞의 패턴 반복
# | = OR, 또는을 뜻함
gct_sample <- gsub("^.*_D|.count.txt","",gctfh_files)
## 1) ^.*_D = GSM4413011_D
## 2) .count.txt = .count.txt
head(gct_sample)

#### 샘플 이름에 따른 그룹 정보 추출
gctfh_info <- function(x){
  # x = "08GCTFH-01"
  sample_celltype <- "GCTFH"
  sample_number <- substring(x,9,10)
  # substring("08GCTFH-01",9,10) # "01"
  sample_day <- substring(x,1,2)
  # substring("08GCTFH-01",1,2) # "08"
  result <- c(sample_celltype,sample_number,sample_day)
  # result # c("GCTFH", "01", "08")
  return(result)

}
#### 추출된 그룹 정보를 테이블로 형성
head(gct_sample)
sample_to_info <- sapply(gct_sample,gctfh_info)
head(sample_to_info)

gct_group_table <- data.table(t(sample_to_info),keep.rownames = T)
head(gct_group_table)

colnames(gct_group_table) <- c("sample","celltype","mouse","day")
head(gct_group_table)


#####################################################################################
#### CD90high GCTfh-like cell 샘플들의 파일 경로 불러오기 및 파일명 기반 샘플 이름 지정
tfh_files_path <- list.files("1.Data/raw_count_txt/TFHlike//",full.names=T) #
head(tfh_files_path)

tfh_files <- list.files("1.Data/raw_count_txt/TFHlike//") #
head(tfh_files)

tfh_sample <- gsub("^.*_D|.count.txt","",tfh_files) #
head(tfh_sample)

#### 샘플 이름에 따른 그룹 정보 추출
tfh_info <- function(x){
  sample_celltype <- "TFHlike"
  # x
  # 08TFHlike-01
  # 1234567890ab
  sample_number <- substring(x,11,12) #
  sample_day <- substring(x,1,2) #
  result <- c(sample_celltype,sample_number,sample_day)
  return(result)
}

#### 추출된 그룹 정보를 테이블로 형성
tfhlike_group_table <- data.table(t(sapply(tfh_sample,tfh_info)),keep.rownames = T)
colnames(tfhlike_group_table) <- c("sample","celltype","mouse","day")

#### GCTfh cell, GCTfh-like cell 그룹 테이블 병합
group_table <- rbind(gct_group_table,tfhlike_group_table)
head(group_table)
file_path <- c(gctfh_files_path,
               tfh_files_path)
head(file_path)
#### GCTfh cell, GCTfh-like cell 그룹 정보 테이블에 파일 경로 삽입
group_table <- cbind(group_table,file_path)
head(group_table)

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

head(count_table)
#### 최소 10 이상 count 된 유전자만 남김
filtered_row <- rowSums(count_table[,-1] >= 10) > 0
head(filtered_row)

filtered_count_table <- count_table[filtered_row,]
head(filtered_count_table)

#### 유전자 이름이 ensembl gene id로 되어 있어 MGI symbol로 변환

# ensembl <- useEnsembl(biomart = "ensembl")
# ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
# saveRDS(ensembl,"3.Results/ensembl.rds")

ensembl <- readRDS("3.Results/ensembl.rds")
head(ensembl@attributes$name)
ensembl@attributes$name[grep("symbol",ensembl@attributes$name)]
# Mouse Genome Informatics
# http://www.informatics.jax.org/

# G_list <- getBM(filters ="ensembl_gene_id",
#                 attributes= c("ensembl_gene_id","mgi_symbol"),
#                 values=filtered_count_table$ensembl_gene_id,
#                 mart= ensembl)
# saveRDS(G_list,"3.Results/g_list.rds")
G_list <- readRDS("3.Results/g_list.rds")
head(G_list)

## 매치된 gene symbol 정보와 count table 병합
final_count_table <- merge(G_list,filtered_count_table,by="ensembl_gene_id",all=T)
head(final_count_table)


## 파일 저장: 그룹 테이블
fwrite(group_table,file = "1.Data/group_table.txt",sep="\t")
## 파일 저장: count 테이블
fwrite(final_count_table,file = "1.Data/count_table.txt",sep="\t")

