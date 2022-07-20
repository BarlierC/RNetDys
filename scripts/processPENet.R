#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)

pred <- fread(args[1])
pred <- pred[!duplicated(pred[,c("Source","Target")]),]
#For benchmarking purpose - removed
pred <- pred[which(pred$Score > 10),]
pred <- pred[which(pred$EnhScore >= 0.2),]

pred_bed <- pred[,c(1,2,3)]
pred_bed$Chr <- as.character(sapply(pred_bed$Source,function(x) strsplit(x,":")[[1]][1]))
pred_bed$Start <- as.numeric(sapply(pred$Source, function(x) strsplit(strsplit(x,":")[[1]][2],"-")[[1]][1]))
pred_bed$End <- as.numeric(sapply(pred$Source,function(x) strsplit(x,"-")[[1]][2]))
pred_bed <- pred_bed[,c(4,5,6,2,3)]

write.table(pred_bed,args[2],sep="\t",row.names=F,col.names=F,quote=F)
