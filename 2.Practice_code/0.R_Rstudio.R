x <- c("a","b","c")
x
z <- c(1,2,3,4,5,6,7,8,9,10)
z
i <- c("1","2","3","4","5")
i
j <- c("a",1,2)
j

coldata <- read.table("1.Data/group_table.txt", sep = "\t", header = T)
countdata <- read.table("1.Data/count_table.txt", sep = "\t", header = T)
countdata
?read.table
head(coldata)
?head
head(coldata)
coldata[1,1]
coldata[1,2]
coldata[2,2]
