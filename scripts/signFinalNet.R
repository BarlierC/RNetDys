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

#Load final network: 5 columns Source, Target, IntType, CorV, Sign
net <- fread(args[2])

#Correlation between all TFs and Genes/promoters
mtx <- mtx[which(rownames(mtx) %in% unique(c(net$Source,net$Target))),]
corMtx <- cor(t(mtx),method="pearson")
dfCor <- reshape2::melt(corMtx, as.is = TRUE)
colnames(dfCor) <- c("TF","Promoter","CorV")

#Data.frame TF-Enhancer-Promoter
dfTFe <- net[which(net$IntType == "TF-Enhancer"),] #TF-Enhancer
dfTFe <- dfTFe[,1:2]
colnames(dfTFe) <- c("TF","Enhancer")
dfep <- net[which(net$IntType == "Enhancer-Promoter"),] #Enhancer-Promoter
dfep <- dfep[,1:2]
colnames(dfep) <- c("Enhancer","Promoter")
dfTFep <- merge(dfTFe,dfep,by="Enhancer",allow.cartesian=TRUE)
rm(dfTFe)
rm(dfep)

#Add correlation
dfTFepS <- merge(dfTFep,dfCor,by=c("TF","Promoter"),all.x=T)

#Define sum correlation & fill network
nettfg <- net[which(net$IntType == "TF-Promoter"),]
nettfg <- nettfg[,c("Source","Target","CorV","IntType")]
net <- net[which(net$IntType %in% c("Enhancer-Promoter","TF-Enhancer")),]
net <- net[,c("Source","Target","IntType")]

#Enhancer-Promoter
ep <- aggregate(dfTFepS[,"CorV"],by=list("Enhancer"=dfTFepS$Enhancer,"Promoter"=dfTFepS$Promoter),sum)
colnames(ep) <- c("Source","Target","CorV")
ep$IntType <- "Enhancer-Promoter"
fnet = rbind(nettfg,ep)
rm(ep)
#TF-Enhancer
tfe <- aggregate(dfTFepS[,"CorV"],by=list("TF"=dfTFepS$TF,"Enhancer"=dfTFepS$Enhancer),sum)
colnames(tfe) <- c("Source","Target","CorV")
tfe$IntType <- "TF-Enhancer"
fnet = rbind(fnet,tfe)
rm(tfe)
rm(dfTFepS)
rm(nettfg)

#Fill in Sign
fnet$Sign <- rep(NA,nrow(fnet))
fnet$Sign[which(fnet$CorV < 0)] <- "-"
fnet$Sign[which(fnet$CorV > 0)] <- "+"

#Save final net
write.table(fnet,args[3],sep="\t",row.names=F,col.names=T,quote=F)