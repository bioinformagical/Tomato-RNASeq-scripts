### Installing Required Packages (only done once)
#source("https://bioconductor.org/biocLite.R")
BiocManager::install("EBSeq")
BiocManager::install("EBSeqHMM")
install.packages("blockmodeling")
install.packages("tidyverse")

### Load Required Pacakges
#library(data.table)
#library(DESeq2)
library(RColorBrewer)
library(EBSeq)
library(EBSeqHMM)
library(tidyverse)
#install.packages("datawizard")
library(datawizard)

## The Purge
rm(list = ls())

### Set WorkDir and load my own R functions
setwd("/Users/robreid/Dropbox (UNC Charlotte)/R/EB-seq_Muday_timecourse/")
source("/Users/robreid/Dropbox (UNC Charlotte)/R/EB-seq_Muday_timecourse/RNAseq_fnxs.R")



##Load data first
## Based on vignette:   https://www.biostat.wisc.edu/~ningleng/EBSeqHMM/EBSeqHMM_vignette_v1.pdf
rawCounts = read_delim('/Users/robreid/bitbucket/flavonoid-rnaseq/72_F3H_PollenTube/results/muday-144-SL4_counts-salmon.txt',delim="\t", col_names = TRUE)
str(rawCounts)

## save description
desc = rawCounts$description
geneid=rawCounts$gene_name


#### Need to make a "Conditions" Vector that describes what time point each 
###
###.   THAT THIS MIGHT BREAK ON OTHER INPUT DATA, WE END UP SORTING ALPHABETICAL
###
###
CondVector=c("T15","T15","T15","T30","T30","T30","T45","T45","T45","T75","T75","T75")
print(CondVector)
#making time points for each of the samples
Conditions=factor(CondVector, levels=c("T15","T30","T45","T75"))
str(Conditions)
#Downstream analysis by EBSeq-HMM requires the conditions to be specified as a factor. In particular,
#levels of the factor need to be sorted along the time/spatial course
levels(Conditions)


###### Choosing a subset via Tidyverse!!!
x = as_tibble(rawCounts)
a28<-select(x,matches("A.28"))
a34<-select(x,matches("A.34"))
v28<-select(x,matches("V.28"))
v34<-select(x,matches("V.34"))
f28<-select(x,matches("F.28"))
f34<-select(x,matches("F.34"))
#a28_tbl <- select(x,matches("A.28"))

#y = as_tibble(rownames_to_column(justCounts)) %>%
#  select(matches("A.28"))
#select(matches("row"),matches("A.28")) 

fetch_UPDownGeneCalls <- function(tibby) {
  dim(tibby)
  tibby <- tibby[,order(colnames(tibby))]  ##Reordering the columns to match the conditions
  mat <- as.data.frame(tibby)
  rownames(mat) <- geneid
  head(mat)
  Sizes=MedianNorm(mat)
  qSizes=QuantileNorm(mat,.75)
  print(qSizes)
  GeneNormData=GetNormalizedMat(mat, Sizes)
  EBSeqHMMGeneOut=EBSeqHMMTest(Data=data.matrix(mat), sizeFactors=Sizes, Conditions=Conditions,UpdateRd=50)
  dim(EBSeqHMMGeneOut)
  GeneDECalls=GetDECalls(EBSeqHMMGeneOut, FDR=.05)
  
  #print the gene COnf calls summary
  GeneConfCalls=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=T)
  print(GeneConfCalls$NumEach) ##We print these, but results are not returned. Should make a separate function.
  list_data <- list(GeneDECalls,GeneConfCalls$NumEach,GeneNormData)
  return(list_data)
}


a28gc <- fetch_UPDownGeneCalls(a28)
a34gc <- fetch_UPDownGeneCalls(a34)
v28gc <- fetch_UPDownGeneCalls(v28)
v34gc <- fetch_UPDownGeneCalls(v34)
f28gc <- fetch_UPDownGeneCalls(f28)
f34gc <- fetch_UPDownGeneCalls(f34)

write.table(a28gc[1], file=sprintf("a28gc_ebseqGeneCalls-SL4.txt"))
write.table(a34gc[1], file=sprintf("a34gc_ebseqGeneCalls-SL4.txt"))
write.table(v28gc[1], file=sprintf("v28gc_ebseqGeneCalls-SL4.txt"))
write.table(v34gc[1], file=sprintf("v34gc_ebseqGeneCalls-SL4.txt"))
write.table(f28gc[1], file=sprintf("f28gc_ebseqGeneCalls-SL4.txt"))
write.table(f34gc[1], file=sprintf("f34gc_ebseqGeneCalls-SL4.txt"))

## Write out Number of genes per path
write.table(a28gc[2], file=sprintf("a28gc_numIneachPath-SL4.txt"))
write.table(a34gc[2], file=sprintf("a34gc_numIneachPath-SL4.txt"))
write.table(v28gc[2], file=sprintf("v28gc_numIneachPath-SL4.txt"))
write.table(v34gc[2], file=sprintf("v34gc_numIneachPath-SL4.txt"))
write.table(f28gc[2], file=sprintf("f28gc_numIneachPath-SL4.txt"))
write.table(f34gc[2], file=sprintf("f34gc_numIneachPath-SL4.txt"))



fetch_EEgeneCalls <- function(tibby) {
  dim(tibby)
  tibby <- tibby[,order(colnames(tibby))]  ##Reordering the columns to match the conditions
  mat <- as.data.frame(tibby)
  rownames(mat) <- geneid
  head(mat)
  Sizes=MedianNorm(mat)
  qSizes=QuantileNorm(mat,.75)
  print(qSizes)
  GeneNormData=GetNormalizedMat(mat, Sizes)
  EBSeqHMMGeneOut=EBSeqHMMTest(Data=data.matrix(mat), sizeFactors=Sizes, Conditions=Conditions,UpdateRd=50)
  dim(EBSeqHMMGeneOut)
  GeneDECalls=GetDECalls(EBSeqHMMGeneOut, FDR=.05)
  
  #Even expression
  AllPaths=GetAllPaths(EBSeqHMMGeneOut)
  print(AllPaths)
  
  #print the gene Conf calls summary
  GeneConfCalls=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=F)
  print(GeneConfCalls$NumEach) ##We print these, but results are not returned. Should make a separate function.
  list_data <- list(GeneDECalls,GeneConfCalls$NumEach,GeneNormData)
  return(list_data)
}


testy <- fetch_EEgeneCalls(a34)
a34m <- as.data.frame(a34)
rownames(a34m) <- geneid
Sizes=MedianNorm(a34m)
EBSeqHMMGeneOut=EBSeqHMMTest(Data=data.matrix(a34m), sizeFactors=Sizes, Conditions=Conditions,UpdateRd=2)
GeneDECalls=GetDECalls(EBSeqHMMGeneOut, FDR=.05,)


a28mat <- as.data.frame(a28)
rownames(a28mat) <- geneid
head(a28mat)
a28mat <- a28mat[,order(colnames(a28mat))]
EBSeqHMMGeneOut=EBSeqHMMTest(Data=data.matrix(a34m), sizeFactors=Sizes, Conditions=Conditions,UpdateRd=2)

#To convert tibble back to data frame
####.   x <- column_to_rownames(mtcars_tbl)
#a28mat <- column_to_rownames(a28_tbl)

  
#EBSeq-HMM requires library size factors to adjust for sequencing depth differences among different samples.
Sizes=MedianNorm(a28mat)
#Could also do quartiles:
#QuantileNorm(GeneMat,.75)
qSizes=QuantileNorm(a28mat,.75)

## Collect the size factors as text files. Can see which normalizing has most effect.
write.table(Sizes, "Normalized_Matrix.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
write.table(qSizes, "QuantileMatrix.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

GeneNormData=GetNormalizedMat(a28mat, Sizes)
#Suppose we are particularly interested in Gene_23, we may apply:
PlotExp(GeneNormData, Conditions, Name="PRAM_106382")

#Running EBSeq-HMM on gene expression estimates
# may be used to estimate parameters and the posterior probability (PP) of being
#in each expression path. For example, here we run five iterations of the Baum-Welch algorithm by setting UpdateRd=5
EBSeqHMMGeneOut=EBSeqHMMTest(Data=data.matrix(a28mat), sizeFactors=Sizes, Conditions=Conditions,UpdateRd=50)

### Detection of DE genes and inference of gene’s most likely path
#DE genes are defined as those showing significant change in at least one condition. Under a target FDR = 0.05, we call genes with
#  posterior probability PP(remain constant)< 0.05 as DE genes.
GeneDECalls=GetDECalls(EBSeqHMMGeneOut, FDR=.05)
head(GeneDECalls)
#DE under 5% target FDR.
str(GeneDECalls)

#######   Clustering DE genes into expression paths ###########
#To cluster DE genes into expression paths, we consider DE genes with confident assignments. By default,
#a gene will be called as a confident assignment to its most likely path if its maximum PP is greater than
#0.5. A user may change this threshold by specifying cutoff.

GeneConfCalls=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=T)
#str(GeneConfCalls$EachPath)
print(GeneConfCalls$EachPath[1:4])

Path4=GeneConfCalls$EachPath[["Down-Down-Up"]]
print(Path4)

#See how many genes in each cluster
print(GeneConfCalls$NumEach)










#To check if a gene is in the list of FDR surviving DE
"Gene_23"%in%rownames(GeneDECalls)
GeneDECalls["Gene_23",]

#######   Clustering DE genes into expression paths ###########
#To cluster DE genes into expression paths, we consider DE genes with confident assignments. By default,
#a gene will be called as a confident assignment to its most likely path if its maximum PP is greater than
#0.5. A user may change this threshold by specifying cutoff.

GeneConfCalls=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=T)
#str(GeneConfCalls$EachPath)
print(GeneConfCalls$EachPath[1:4])

Path4=GeneConfCalls$EachPath[["Down-Down-Up"]]
print(Path4)

#See how many genes in each cluster
print(GeneConfCalls$NumEach)
str(GeneConfCalls$EachPathNames)

#Assume we are only interested in the monotone increasing path and monotone decreasing path (1st
# and 16th in the list of dynamic paths), we may define Paths=AllPaths[c(1,16)] in GetConfidentCalls
# function to obtain the simplified list:
AllPaths=GetAllPaths(EBSeqHMMGeneOut)
print(AllPaths)
GeneConfCallsTwoPaths=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, Paths=AllPaths[c(1,2)])
print(GeneConfCallsTwoPaths)

## Data Visualization
# /*
# EBSeq-HMM also provides functions to visualize genes/isoforms of interest.
# Prior to generating the plots, we first need to obtain a normalized expression matrix that adjusts for library
# sizes. The normalized matrix may be obtained by the GetNormalizedMat function:
# */
AllPathsWithEE=GetAllPaths(EBSeqHMMGeneOut, OnlyDynamic=F)
print(AllPathsWithEE)

GeneConfCallsTwoPaths=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05, Paths=AllPaths)
print(GeneConfCallsTwoPaths)

GeneNormData=GetNormalizedMat(a28mat, Sizes)
print(GeneConfCallsTwoPaths$EachPath[["Up-Up-Down"]])
#GeneOfInterest=GeneConfCallsTwoPaths$EachPathGeneNames[["Down-Down-Down-Down"]]
#print(GeneOfInterest)

#diagnostic QQ plot
par(mfrow=c(3,2))
QQP(EBSeqHMMGeneOut,GeneLevel=T)

par(mfrow=c(4,4))
DenNHist(EBSeqHMMGeneOut, GeneLevel=T)

