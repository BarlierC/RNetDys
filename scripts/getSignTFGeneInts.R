#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stats)
library(reshape2)

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

#Load network: 2 columns Source, Target
net <- fread(args[2])

#Correlation between TFs (Source) and Genes (Targets)
mtx <- mtx[which(rownames(mtx) %in% unique(c(net$Source,net$Target))),]
corMtx <- cor(t(mtx),method="pearson")
dfCor <- reshape2::melt(corMtx)
colnames(dfCor) <- c("Source","Target","CorV")

#Merge with TF-Gene network
netF <- merge(net,dfCor,by=c("Source","Target"),all.x=T)

#Save
write.table(netF,args[3],sep="\t",row.names=F,col.names=T,quote=F)

