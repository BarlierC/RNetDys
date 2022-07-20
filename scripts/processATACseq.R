#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(stringr)
library(data.table)

#Load matrix - Important to have it cleaned first !
atac_mtx <- read.table(args[1],header=T,row.names=1)

#Extract rownames(peaks coordinates to transform into bed file, coordinates needs to be such as chr1:1234-3451
if(colnames(atac_mtx)[1] == "X"){
	peaks <- atac_mtx$X
}else{
	peaks <- rownames(atac_mtx)
}

#Bed file
bed_df <- data.frame("Chr"=rep(NA,length(peaks)),"Start"=rep(NA,length(peaks)),"Stop"=rep(NA,length(peaks)))
bed_df$Chr <- sapply(peaks,function(x){
		strsplit(x,split=":")[[1]][1]
	})
bed_df$Start <- sapply(peaks,function(x){
		strsplit(strsplit(x,split=":")[[1]][2],split="-")[[1]][1]
	})
bed_df$Stop <- sapply(peaks,function(x){
		strsplit(strsplit(x,split=":")[[1]][2],split="-")[[1]][2]
	})

write.table(bed_df,args[2],sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)