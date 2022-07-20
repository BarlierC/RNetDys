#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)

#Intersected enhancer/atac/chip file
bed <- fread(args[1])
colnames(bed) <- c("EnhChr","EnhStart","EnhEnd","ChipChr","ChipStart","ChipEnd","TargetGene","TFB")

#TFGene Network
net = fread(args[2])

#Genes expressed
expg <- unique(c(net$Source,net$Target))

#Keep only if TF expressed
net <- bed[which(bed$TFB %in% expg),]

#Make unique
net <- net[,c(1,2,3,8)]
rm(bed)
gc()
net <- net[!duplicated(net[,c(1,2,3,4)]),]

write.table(net,args[3],sep="\t",row.names=F,col.names=T,quote=F)
