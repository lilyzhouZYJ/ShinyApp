library(rhdf5)
library(tidyverse)
library(tibble)
library(dplyr)

#load the hdf5 file
files <- list.files(pattern="*.hdf5", path = 'hdf5/', full.names = TRUE)
intervals <- h5read(files[1], name = "intervals/transposed_index_start_end")

df <- tibble(chr = (intervals[,1] + 1), start = intervals[,2], stop = intervals[,3])
for (file in files){
  a <- h5read(file, name = "counts/values")
  pref <- str_extract(file[1], 'M_[^.]+')
  df[[pref]] <- a
}

#load into df2 + calculate rpkm
df2 <- tibble(chr = (intervals[,1] + 1), start = intervals[,2], stop = intervals[,3])
for (file in files){
  a <- h5read(file, name = "counts/values")
  total <- sum(a)
  pref <- str_extract(file[1], 'M_[^.]+')
  df2[[pref]] <- a*(10^9)/total
}
df2 <- df2 %>% mutate(targetsize = stop - start)
df2[,4:174] <- df2[,4:174]/df2$targetsize

#filter based on median
df2[["median"]] <- apply(df2[,4:174],1,median)
df2 <- filter(df2, df2$median>=1)

#calculate z-score
df2[["sd"]] <- apply(df2[,4:174], 1, sd)
df2[["mean"]] <- apply(df2[,4:174], 1, mean)
df2[,4:174] <- (df2[,4:174]-df2$mean)/df2$sd

#SVD
mat <- df2[-c(1:3,175:178)]  #leave only relevant data values
mat <- data.frame(t(mat))  #transpose into exon by sample
s <- svd(mat)  #svd calculation

#scree plot of singular values to determine the number of components to delete
scree <- data.frame(d=s$d, num=1:171)
ggplot(scree, aes(x=num, y=d)) + geom_line() + geom_point() + geom_smooth(formula="y~log(x)")

#delete no components
D <- diag(s$d)
df3 <- s$u %*% D %*% t(s$v)
df3 <- t(df3)
cnames <- colnames(df2)[4:174]
colnames(df3) <- cnames
head <- df2[1:3]
df3 <- cbind(head, df3)
df3$index <- c(1:nrow(df3))
dataset <- df3 %>% unite(exon, chr, start, stop, sep="_") %>% gather(sample, zrpkm, -exon, -index)
dataset <- cbind(head,dataset)

#svd: delete 10 components
s$d10 <- s$d
s$d10[c(1:10)] <- 0
D10 <- diag(s$d10)
df3_10 <- s$u %*% D10 %*% t(s$v)
df3_10 <- t(df3_10)
cnames <- colnames(df2)[4:174]
colnames(df3_10) <- cnames
head <- df2[1:3]
df3_10 <- cbind(head, df3_10)
df3_10$index <- c(1:nrow(df3_10))
dataset_10 <- df3_10 %>% unite(exon, chr, start, stop, sep="_") %>% gather(sample, zrpkm, -exon, -index)
dataset_10 <- cbind(head,dataset_10)

#split dataset
dataset <- as_tibble(dataset)
dataset_level <- unique(dataset$exon)
dataset$exon <- factor(dataset$exon, levels=dataset_level)
dataset_10 <- as_tibble(dataset_10)
dataset_10_level <- unique(dataset_10$exon)
dataset_10$exon <- factor(dataset_10$exon, levels=dataset_10_level)
for (i in 1:24){
  name <- paste('dataset_chr', i, sep = '')
  assign(name, filter(dataset, chr == i))
}
for (i in 1:24){
  name <- paste('dataset_10_chr', i, sep = '')
  assign(name, filter(dataset_10, chr == i))
}


#genes to regions
install.packages("BiocManager")
BiocManager::install('biomaRt')
library('biomaRt')

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

chromosomes=c(1:22,"X","Y")
gene_list <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                  'start_position', 'end_position'),
                   filters = c('chromosome_name','with_hgnc'), 
                   values = list(chromosomes,TRUE),
                   mart = mart)
gene_list <- group_by(as_tibble(gene_list), chromosome_name)
gene_list <- arrange(gene_list, start_position, .by_group = TRUE)


       