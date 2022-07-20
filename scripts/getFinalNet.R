#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(igraph)

#TF-Gene (promoter region) net
TFGenenet <- fread(args[1])

#Enhancer-promoter net
EPNet <- fread(args[2])

#TF - Enhancer net
TFENet <- fread(args[3])
TFENet$Enh <- paste(TFENet$EnhChr,paste(TFENet$EnhStart,TFENet$EnhEnd,sep="-"),sep=":")
TFENet <- TFENet[,c("TFB","Enh")]
colnames(TFENet) <- c("Source","Target")

#Directed combined Network

#TF-Gene
TFGenenet$IntType <- "TF-Promoter" #TF-Promoter region interactions
TFGenenet$Sign <- NA
TFGenenet$Sign[which(TFGenenet$CorV < 0)] <- "-"
TFGenenet$Sign[which(TFGenenet$CorV > 0)] <- "+"
TFGenenet <- TFGenenet[,c("Source","Target","IntType","CorV","Sign")]

#TF-Enhancer
TFENet$IntType <- "TF-Enhancer" #TF-Enhancer region interactions
#Get list of enhancer for which at least one TF is binding = active enhancer
active_enhancers <- unique(TFENet$Target)
on_promot <- unique(c(TFGenenet$Target,TFGenenet$Source))
TFENet$CorV <- NA
TFENet$Sign <- NA
TFENet <- TFENet[,c("Source","Target","IntType","CorV","Sign")]

#Enhancer-Promoter
EPNet <- EPNet[,c(1,2,5)] #Enhancer-Promoter interactions
EPNet$IntType <- "Enhancer-Promoter"
EPNet$CorV <- NA
EPNet$Sign <- NA
EPNet <- EPNet[,c("Source","Target","IntType","CorV","Sign")]
#Remove enhancer-promoter interaction for which enhancer is not assumed to be active (no TF binding)
EPNet <- EPNet[which(EPNet$Source %in% active_enhancers),]
#Remove enhancer-promoter interaction for which promoter gene is not expressed
EPNet <- EPNet[which(EPNet$Target %in% on_promot),]
#Remove TF-Enhancers for which enhancer is not found
TFENet <- TFENet[which(TFENet$Target %in% EPNet$Source),]

#Combination of subnetworks
S2EPreg_net <- rbind(TFGenenet,TFENet)
S2EPreg_net <- rbind(S2EPreg_net,EPNet)

#Final Network
write.table(S2EPreg_net,args[4],sep="\t",row.names=F,col.names=T,quote=F)
