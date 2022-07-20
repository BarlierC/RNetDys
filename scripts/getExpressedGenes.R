#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)

#Get genes conserved at least in PercExp % of the cells
#Return filtered matrix
dt_filter_gene_expression <- function(dt,PercExp){
  #Binarize data
  dtb <- dt
  dtb[dtb>0] <- 1
  #Count for each gene number of cell expressing it: percentage
  cg <- apply(dtb,1,function(x){sum(x)/length(dtb)*100})
  names(cg) <- rownames(dtb)
  #Keep only the ones stricly greater than the threshold
  cg_tokeep <- which(cg>PercExp)
  #Filter the data
  dt <- dt[cg_tokeep,]
  return(dt)
}

#Load scRNAseq matrice: genes in rows (gene name) and cells in columns
df <- fread(args[1])
df <- as.data.frame(df)
colnames(df)[1] <- "Genes"
df <- df[!duplicated(df[,c("Genes")]),]
rownames(df) <- df$Genes
df$Genes <- NULL
mtx <- as.matrix(df)
#Basic QC
#Remove cells with no expression
mtx <- mtx[rowSums(mtx)>0,]
#Remove genes with no expression
mtx <- mtx[,colSums(mtx)>0]
df_filtered <- dt_filter_gene_expression(df,50)
genes_expressed <- rownames(df_filtered)

#Load network
net <- fread(args[2])

#Format Network (Source,Target,[...])
chipatac_net <- data.table("Source"=net[,8],"Target"=net[,7])
colnames(chipatac_net) <- c("Source","Target")
chipatac_net <- chipatac_net[!duplicated(chipatac_net[,c("Source","Target")]),]

#Keep only genes found as expressed in the network
chipatac_net_expGenes <- chipatac_net[which(chipatac_net$Source %in% genes_expressed & chipatac_net$Target %in% genes_expressed),]

write.table(chipatac_net_expGenes,args[3],sep="\t",row.names=F,col.names=T,quote=F)
