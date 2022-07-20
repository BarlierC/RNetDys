#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(reshape2)

### Network
net = fread(args[1],header=T)

#Promoters
getfs <- unique(c(
          net$Target[which(net$IntType == "TF-Promoter")],
          net$Target[which(net$IntType == "Enhancer-Promoter")]
          ))

### ChIP-seq
chip = fread(args[2],header=F)

#Filter promoters of the network
chip <- chip[which(chip[,4] %in% getfs),]
chip <- chip[,c(1,2,3,4)]
chip <- chip[!duplicated(chip[,1:4]),]
colnames(chip) <- c("Chr","Start","Stop","Gene")

#Get enhancers of the network
enhNet <- net[which(net$IntType == "Enhancer-Promoter"),]
enh_df <- data.frame("Chr"=rep(NA,nrow(enhNet)),
                     "Start"=rep(NA,nrow(enhNet)),
                     "Stop"=rep(NA,nrow(enhNet)),
                     "Gene"=rep("ENHANCER",nrow(enhNet)))
enh_df$Chr <- sapply(enhNet$Source,function(x){strsplit(x,split=":")[[1]][1]})
enh_df$Start <- sapply(enhNet$Source,function(x){strsplit(strsplit(x,split=":")[[1]][2],split="-")[[1]][1]})
enh_df$Stop <- sapply(enhNet$Source,function(x){strsplit(strsplit(x,split=":")[[1]][2],split="-")[[1]][2]})
enh_df <- enh_df[!duplicated(enh_df[,c("Chr","Start","Stop","Gene")]),]

#Combine promoters & enhancers regions
coord <- rbind(chip,enh_df)

#Bed file of enhancers
write.table(coord,args[3],sep="\t",row.names=F,col.names=F,quote=F)
