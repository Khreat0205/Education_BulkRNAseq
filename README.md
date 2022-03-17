# Education_BulkRNAseq
목표:

## I. 분석 예제 데이터

GEO series dataset: [GSE147035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147035 ) 

논문: Primary germinal center-resident T follicular helper cells are a physiologically distinct subset of CXCR5hiPD-1hi T follicular helper cells. Yeh et al. *Immunity*. 2022. ([링크](https://www.sciencedirect.com/science/article/pii/S1074761321005513))

다운 받아야할 파일: R 설치 파일, Rstudio 설치 파일, Count matrix, tar/gz 압축 풀수 있는 반디집, 패키지 설치 스크립트 (?)



## II. 분석 내용	

1. **GEO Data processing**
2. **Data Description**: PCA, Correlation analysis. (+ Quality check)
3. **Differentially Expressed Gene analysis**: DESeq2
4. **Visualization**: Volcano plot, Heatmap plot
5. **Gene Ontology analysis**: DAVID, Toppgene



## III. 실습 전 준비사항

1. 사전 설치 필요: 0.Pre-install 폴더 내 파일 활용하여 설치
   1. R
   2. Rstudio
   3. R packages (DESeq2, BiocManager,data.table, biomaRt)
2. 사전 다운로드
   1. 이 git을 압축파일로 다운로드 하거나, 링크(드랍박스)에서 다운로드



## IV. 진행 순서

1. 실습생 기본 환경 세팅 확인 
2. GEO 소개 및 다운로드, 처리 방법 (in R)
3. PCA 분석, 상관관계 그래프를 통한 데이터 확인 (in R)
4. DESeq2를 활용한 Differentially expressed genes 분석 (in R)
5. DEG 분석 시각화 (Volcano plot, Heatmap plot) (in R)
6. Gene ontology 및 관련 도구 소개, Gene ontology 분석 (in Web)
7. Q & A



## V. Contacts

scientist0205@snu.ac.kr

