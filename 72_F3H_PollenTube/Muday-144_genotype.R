# Clear & Prep Workspace 
rm(list = ls())

library("DESeq2")
library(ggplot2)
library(ggrepel)
library(htmltools)
library(DESeq2)
library(stringr)
library(dplyr)
library("gridExtra")
library(ggpubr)
library(ggforce)

counts = read.csv('results/muday-144-SL5_counts-salmon.txt',sep = "\t",stringsAsFactors = F,
                  header=T, check.names = FALSE, row.names = "gene_name")
counts = counts[,-ncol(counts)]



## MetaData funtion
library(stringr)
makeExperimentMetaDataFrame = function(fname="results/muday-144-SL5_counts-salmon.txt") {
  experiment_data = read.csv(fname,
                             sep = "\t",
                             stringsAsFactors = F, 
                             header=T, 
                             check.names = FALSE, 
                             row.names = "gene_name")
  sample_names = colnames(experiment_data)
  boolean_vector = str_detect(sample_names,'[VFA]\\.\\d\\d\\.\\d\\d\\.\\d')
  sample_names = sample_names[boolean_vector]
  genotype = sapply(strsplit(sample_names,"\\."),function(x){x[[1]]})
  time = sapply(strsplit(sample_names,"\\."),function(x){x[[3]]})
  temperature = sapply(strsplit(sample_names,"\\."),function(x){x[[2]]})
  to_return = data.frame(genotype,time,temperature)
  rownames(to_return) = sample_names
  return(to_return)
}

## DESeq function 
getDEgenes <- function(counts, sampleData) {
  cts <- as.matrix(counts) 
  coldata <- sampleData
  coldata <- coldata[,c("genotype", "time", "temperature")]
  coldata$genotype <- factor(coldata$genotype, levels = c("A", "F", "V"))
  coldata$time <- factor(coldata$time, levels = c("15", "30", "45", "75"))
  coldata$temperature <- factor(coldata$temperature, levels = c("28", "34"))
  
  cts <- cts[, rownames(coldata)]
  dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                colData = coldata,
                                design = ~ genotype + temperature + genotype:temperature)
  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  resultsNames(dds)
  res05 <- results(dds, alpha=0.05)
  numSignGenes<- sum(res05$pvalue < 0.05, na.rm=TRUE) # Not adjusted
  numSignGenes<- sum(res05$padj < 0.05, na.rm=TRUE) # Adjusted
  res05 <- res05[order(res05$pvalue),]
  res05top25 <- res05[1:25,]
  res05top25 
  vsd <- vst(dds, blind=FALSE)
}


## PCA Plot function 

PCA_plots <- function(vsd, sampleData, title) {
  options(ggrepel.max.overlaps = Inf)
  pcaData <- plotPCA(vsd, intgroup=c("temperature", "genotype", "time"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=time, shape=genotype)) +
    #geom_label_repel(aes(label = rownames(sampleData))) +
    scale_color_manual(values =  c("15" ="#56B4E9","30" = "#009E73",
                                   "45" = "#E69F00", "75" = "#CC79A7")) +
    scale_fill_manual(values =  c("15" ="#56B4E9","30" = "#009E73",
                                  "45" = "#E69F00", "75" = "#CC79A7")) +
    geom_point(size =4, aes(fill =time, alpha=temperature)) +
    geom_point(size=4) +
    scale_shape_manual(values= c("A"= 23, "F"= 22, "V"= 25)) + 
    scale_alpha_manual(values=c("28"=0.1, "34"=1)) +
    #geom_mark_circle(aes(color = as.factor(time)), expand = unit(0.5,"mm"))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(title) +
    theme_bw() +
    theme(aspect.ratio = 1)+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 15))
}

## Run metadata function
sampleData = makeExperimentMetaDataFrame()


##################################################################################################


# View DESeq2 Design line by line

cts <- as.matrix(counts) 
coldata <- sampleData
coldata <- coldata[,c("genotype", "time", "temperature")]

## Should I group the factors????

coldata$genotype <- factor(coldata$genotype, levels = c("A", "F", "V"))
coldata$time <- factor(coldata$time, levels = c("15", "30", "45", "75"))
coldata$temperature <- factor(coldata$temperature, levels = c("28", "34"))
  
cts <- cts[, rownames(coldata)]

## Check Design and Interactions
model.matrix(~genotype + temperature + genotype:temperature, coldata)


## Run analysis
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                colData = coldata,
                                design = ~genotype + temperature + genotype:temperature)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
dds <- DESeq(dds, minReplicatesForReplace=Inf)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
results(dds) # View last output Design with results
resultsNames(dds) # View all of the Design and interactions


## Organizing the data to view specific significant genes
res05 <- results(dds, alpha=0.05)
numSignGenes<- sum(res05$pvalue < 0.05, na.rm=TRUE) # Not adjusted
numSignGenes<- sum(res05$padj < 0.05, na.rm=TRUE) # Adjusted
res05 <- res05[order(-res05$log2FoldChange),] # Order by fold change
res05top25 <- res05[1:10,] # Take top 10 of upregulated genes
res05top25  # Analysis includes all times, genotypes, and temperatures
res05last25 <- tail(res05, n = 10) # Take top 10 of downregulated genes
res05last25

## Heat map of top 10 most upregulated and top 10 most downregulated significant genes?



####################################################################################################

# Run DESeq function with each time point

## 15 minutes
min15 <- select(counts,contains("15"))
sampleData$type <- rownames(sampleData)
meta15 <- dplyr::filter(sampleData, grepl('15', type))
vsd = getDEgenes(min15, meta15)
PCA_15<-PCA_plots(vsd, meta15, "15 minutes All Genotypes")
PCA_15

## 30 Minutes 

min30 <- select(counts,contains("30"))
sampleData$type <- rownames(sampleData)
meta30 <- dplyr::filter(sampleData, grepl('30', type))
vsd = getDEgenes(min30, meta30)
PCA_30<-PCA_plots(vsd, meta30, "30 minutes All Genotypes")
PCA_30


## 45 Minutes 

min45 <- select(counts,contains("45"))
sampleData$type <- rownames(sampleData)
meta45 <- dplyr::filter(sampleData, grepl('45', type))
vsd = getDEgenes(min45, meta45)
PCA_45<-PCA_plots(vsd, meta45, "45 minutes All Genotypes")
PCA_45


## 75 Minutes 

min75 <- select(counts,contains("75"))
sampleData$type <- rownames(sampleData)
meta75 <- dplyr::filter(sampleData, grepl('75', type))
vsd = getDEgenes(min75, meta75)
PCA_75<-PCA_plots(vsd, meta75, "75 minutes All Genotypes")
PCA_75


## All PCA plots combined

combined<- ggpubr::ggarrange(PCA_15, PCA_30, PCA_45, PCA_75, # list of plots
                             labels = "AUTO",
                             font.label = list(size = 30), 
                             common.legend = T, # COMMON LEGEND
                             legend = "right", # legend position
                             align = "hv", # Align them both, horizontal and vertical
                             nrow = 2, 
                             ncol = 2)
combined

## Note: Will need to fix issue with all of colored times not showing up on combined legend 


## Next make sure the DESeq function collects result data to make heatmaps 

#########################################################################################

# Next create heatmaps for each timepoint top 10 results: we want to use the top 10 up-regulated and 10 down-regulated logfoldchange from each time on the x axis for each genotype and on the y-axis the genes?


# Online resources for heatmap: 

## https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-heatmap2/tutorial.html#create-heatmap-of-top-genes

## https://sebastianraschka.com/Articles/heatmaps_in_r.html
