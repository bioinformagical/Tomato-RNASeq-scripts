# Functions that perform statistical analyses

library(stringr)
library(readr)
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)

source("Common.R")


getDeGenes = function(counts,group1_name,group2_name,
                      method="DESeq2") {
  to_return = NULL
  if (method=="DESeq2") {
    to_return = compareTwoGroupsDESeq2(counts,group1_name,
                                        group2_name) 
  }
  if (method=="edgeR") {
     to_return = compareTwoGroupsEdgeR(counts,group1_name,
                                       group2_name)
  }
  return(to_return)
}

compareTwoGroupsDESeq2 = function(counts,group1_name,group2_name) {
  group1_indexes=grep(group1_name,names(counts)) 
  group2_indexes=grep(group2_name,names(counts))
  indexes=c(group1_indexes,group2_indexes)
  
  condition = factor(c(rep(group1_name,length(group1_indexes)),
                       rep(group2_name,length(group2_indexes))))
  m = round(as.matrix(counts[,indexes]))
  dds = DESeqDataSetFromMatrix(m,
                               DataFrame(condition),
                               ~ condition)
  featureData = data.frame(gene=rownames(m))
  mcols(dds) = DataFrame(mcols(dds),featureData)
  dds = DESeq(dds, minReplicatesForReplace=Inf)
  results_table = results(dds,
                          alpha=0.05)
  results_table = results_table[!is.na(results_table$padj),]
  results_table$gene=row.names(results_table)
  results_table$group1=group1_name
  results_table$group2=group2_name
  results_table$description = counts[results_table$gene,
                                     "description"]
  o = order(results_table$pvalue)
  results_table = results_table[o,]
  d = data.frame(results_table)
  return(d)
}


compareTwoGroupsEdgeR  = function(counts,group1_name,group2_name) {
  group1_indexes=grep(group1_name,names(counts)) 
  group2_indexes=grep(group2_name,names(counts))
  indexes=c(group1_indexes,group2_indexes)
  group = factor(c(rep(group1_name,length(group1_indexes)),
                   rep(group2_name,length(group2_indexes))))
  little_DGEList = DGEList(counts=counts[,indexes],group=group)
  keep = filterByExpr(little_DGEList)
  little_DGEList = little_DGEList[keep,,keep.lib.sizes=FALSE]
  little_DGEList = calcNormFactors(little_DGEList)
  design=model.matrix(~group)
  little_DGEList = estimateDisp(little_DGEList,design)
  fit = glmFit(little_DGEList, design)
  lrt=glmLRT(fit,coef=2)
  results=lrt$table
  results$Q=p.adjust(results$PValue,method="fdr")
  results$gene=row.names(results)
  results$group1=group1_name
  results$group2=group2_name
  results$description = counts[results$gene,"description"]
  o = order(results$Q)
  results = results[o,]
}

volcano_plot<- function(d, title, method="DESeq2",Q=0.05) {
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if (method=="DESeq2") {
    lfc_colname = "log2FoldChange"
    padj_colname = "padj"
    p_colname = "pvalue"
  }
  if (method=="edgeR") {
    lfc_colname = "logFC"
    padj_colname = "Q"
    p_colname = "PValue"
  }
  pCutoffCol = padj_colname
  if (!any(d[,pCutoffCol]<=Q)) {
    Q = round(min(d$padj),2)+0.01
  }
  lfcmin = min(d[,lfc_colname])
  lfcmax = max(d[,lfc_colname])
  #x_max = -log10(min(d[,p_colname]))
  keyvals <- ifelse(
    d[,lfc_colname] < -1, '#56B4E9',
    ifelse(d[,lfc_colname] > 1, '#CC79A7','black'))
  names(keyvals)[keyvals == '#CC79A7'] = 'Up-regulated'
  names(keyvals)[keyvals == 'black'] = 'Mid-regulated'
  names(keyvals)[keyvals == '#56B4E9'] = 'Down-regulated'
  d$p_colname = ifelse(d$p_colname < -log10(20),d$p_colname,1e-20)
  
  v = EnhancedVolcano(d,
                      lab = d$gene,
                      x = lfc_colname,
                      y = p_colname,
                      selectLab = d[which(names(keyvals) %in% c('Up-regulated', 'Down-regulated')),
                                    "gene"],
                      colCustom = keyvals,
                      colAlpha = 1,
                      shape = c(4, 1, 6, 3),
                      pCutoff = Q,  
                      FCcutoff = 1.0,
                      pCutoffCol = pCutoffCol, # 
                      gridlines.minor=FALSE, 
                      gridlines.major=FALSE,
                      col = colorBlindGrey8,
                      #xlim = c(-10, 10),
                      #ylim = c(0,pmax),
                      #ylab = "-log(p)",
                      title = title, 
                      subtitle = paste("Horizontal dashed line shows FDR",Q),
                      legendPosition = "right")
  return(v)
}


