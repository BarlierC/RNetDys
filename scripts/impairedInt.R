#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(reshape2)

### Mapped SNPs
mapSnps = fread(args[1],header=F)
if(ncol(mapSnps) == 11){
  #Remove last columns added by intersection
  mapSnps <- mapSnps[,-11]
  colnames(mapSnps) <- c("Chr","Start","Stop","Gene","SNPChr","SNPStart","SNPend","Ref","Alt","MAF")
}else{
  #Remove last columns added by intersection
  mapSnps <- mapSnps[,-12]
  #RS ID
  colnames(mapSnps) <- c("Chr","Start","Stop","Gene","SNPChr","SNPStart","SNPend","Ref","Alt","MAF","RSID")
}
mapSnps$MAF <- NULL
mapSnps$SNP <- paste(mapSnps$SNPChr,paste(mapSnps$SNPStart,mapSnps$SNPend,sep="-"),sep=":")
mapSnps$SNPChr <- NULL
mapSnps$SNPStart <- NULL
mapSnps$SNPend <- NULL
mapSnps$Coord <- paste(mapSnps$Chr,paste(mapSnps$Start,mapSnps$Stop,sep="-"),sep=":")
mapSnps$Chr <- NULL
mapSnps$Start <- NULL
mapSnps$Stop <- NULL

if(ncol(mapSnps)==5){
  mapSnps <- mapSnps[,c("Gene","Coord","SNP","Ref","Alt")]
}else{
  mapSnps <- mapSnps[,c("Gene","Coord","SNP","Ref","Alt","RSID")]
}

#SNP on promoter region
promotSnp <- mapSnps[which(mapSnps$Gene != "ENHANCER"),]
#SNP on enhancer region
enhSnp <- mapSnps[which(mapSnps$Gene == "ENHANCER"),]
enhSnp$Gene <- NULL

### Network
net = fread(args[2],header=T)

# a/ TF-Promoter
netTFpro <- net[which(net$IntType == "TF-Promoter"),]
if(nrow(netTFpro)>0){
	# Add SNPs
  if(ncol(promotSnp) == 6){
    colnames(promotSnp) <- c("Target","Coordinates","SNP","Ref","Alt","RSID")
    netTFpro <- merge(netTFpro,promotSnp,by="Target",allow.cartesian=TRUE,all.x=F,all.y=F)
    colnames(netTFpro) <- c("Target","Source","CorV","IntType","Sign","PromoterCoord","SNP","Ref","Alt","RSID")
    netTFpro <- netTFpro[,c("Source","Target","CorV","IntType","Sign","SNP","Ref","Alt","RSID")]
  }else{
    colnames(promotSnp) <- c("Target","Coordinates","SNP","Ref","Alt")
    netTFpro <- merge(netTFpro,promotSnp,by="Target",allow.cartesian=TRUE,all.x=F,all.y=F)
    colnames(netTFpro) <- c("Target","Source","CorV","IntType","Sign","PromoterCoord","SNP","Ref","Alt")
    netTFpro <- netTFpro[,c("Source","Target","CorV","IntType","Sign","SNP","Ref","Alt")]
  }
}
if(ncol(netTFpro) == 9){
  netTFpro <- netTFpro[!duplicated(netTFpro[,1:9]),]
}else{
  netTFpro <- netTFpro[!duplicated(netTFpro[,1:8]),]
}

# b/ Enhancer-Promoter
netEnhPro <- net[which(net$IntType == "Enhancer-Promoter"),]
if(nrow(netEnhPro)>0){
	# Add SNPs
  if(ncol(enhSnp) == 5){
    colnames(enhSnp) <- c("Source","SNP","Ref","Alt","RSID")
    netEnhPro <- merge(netEnhPro,enhSnp,by="Source",allow.cartesian=TRUE,all.x=F,all.y=F)
    colnames(netEnhPro) <- c("Source","Target","CorV","IntType","Sign","SNP","Ref","Alt","RSID")
  }else{
    colnames(enhSnp) <- c("Source","SNP","Ref","Alt")
    netEnhPro <- merge(netEnhPro,enhSnp,by="Source",allow.cartesian=TRUE,all.x=F,all.y=F)
    colnames(netEnhPro) <- c("Source","Target","CorV","IntType","Sign","SNP","Ref","Alt")
  }
}
if(ncol(netEnhPro) == 9){
  netEnhPro <- netEnhPro[!duplicated(netEnhPro[,1:9]),]
}else{
  netEnhPro <- netEnhPro[!duplicated(netEnhPro[,1:8]),]
}

# Merge the two informations
if(nrow(netTFpro)>0 | nrow(netEnhPro)>0){
	#Order
	if(ncol(enhSnp) == 5){
	  netEnhPro <- netEnhPro[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt","RSID")]
	}else{
	  netEnhPro <- netEnhPro[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt")]
	}
  #TF-enhancer: add TF for which enhancers is impaired
  netTfEnh <- net[which(net$IntType == "TF-Enhancer"),]
  netTfEnh <- netTfEnh[,c("Source","Target","IntType","CorV","Sign")]
  tmp <- netEnhPro
  colnames(tmp)[1:2] <- c("Enh","Gene")
  tmp$IntType <- NULL
  tmp$CorV <- NULL
  tmp$Sign <- NULL
  tmp$Gene <- NULL
  colnames(netTfEnh)[1:2] <- c("TF","Enh")
  netTfEnh <- merge(netTfEnh,tmp,by="Enh")
  colnames(netTfEnh)[1:2] <- c("Target","Source")
  if(ncol(enhSnp) == 5){
    netTfEnh <- netTfEnh[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt","RSID")]
  }else{
    netTfEnh <- netTfEnh[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt")]
  }
}else{
	print("No impaired interactions found")
}

#Bind
if(nrow(netTFpro)>0 & nrow(netEnhPro)>0){
  netImpaired <- rbind(netTFpro,netEnhPro)
  netImpaired <- rbind(netImpaired,netTfEnh)
  #Remove duplicates
  if(ncol(enhSnp) == 5){
    netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt","RSID")]),]
  }else{
    netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt")]),]
  }
}else if(nrow(netTFpro)>0 & nrow(netEnhPro)==0){
  netImpaired <- netTFpro
  #Remove duplicates
  if(ncol(promotSnp) == 6){
    netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt","RSID")]),]
  }else{
    netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt")]),]
  }
}else if(nrow(netTFpro)==0 & nrow(netEnhPro)>0){
  netImpaired <- netEnhPro
  netImpaired <- rbind(netImpaired,netTfEnh)
  #Remove duplicates
  if(ncol(enhSnp) == 5){
    netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt","RSID")]),]
  }else{
    netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt")]),]
  }
}else{
  #Empty df
  netImpaired <- netEnhPro
}

#Impaired interactions to use for binding affinity
write.table(netImpaired,args[3],sep="\t",row.names=F,col.names=T,quote=F)

netImpaired <- netImpaired[!duplicated(netImpaired[,c("Source","Target","IntType","CorV","Sign","SNP","Ref","Alt")]),]
write.table(netImpaired,args[4],sep="\t",row.names=F,col.names=T,quote=F)

