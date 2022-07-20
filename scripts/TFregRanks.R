#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stringr)

#Impaired reg interactions
netImp <- fread(args[1],header=T)

#SNPs info
snpsInfo <- fread(args[2],header=F)
if(ncol(snpsInfo) == 6){
  colnames(snpsInfo) <- c("Chr","Start","Stop","Ref","Alt","MAF")
  snpsInfo$SNP <- paste(snpsInfo$Chr,paste(snpsInfo$Start,snpsInfo$Stop,sep="-"),sep=":")
  snpsInfo <- snpsInfo[,c("SNP","Ref","Alt","MAF")]
  #Merge info
  netImp <- merge(netImp,snpsInfo,by=c("SNP","Ref","Alt"))
}else{
  colnames(snpsInfo) <- c("Chr","Start","Stop","Ref","Alt","MAF","RSID")
  snpsInfo$SNP <- paste(snpsInfo$Chr,paste(snpsInfo$Start,snpsInfo$Stop,sep="-"),sep=":")
  snpsInfo <- snpsInfo[,c("SNP","Ref","Alt","MAF","RSID")]
  netImp <- merge(netImp,snpsInfo,by=c("SNP","Ref","Alt","RSID"))
}

#Unique TFs
utfs <- unique(netImp$Source[which(netImp$IntType != "Enhancer-Promoter")])
#Final df init
finalDf <- data.frame("TF"=rep(NA,length(utfs)),"Score"=rep(NA,length(utfs)),"Rank"=rep(NA,length(utfs)))
netImpf <- netImp[which(netImp$IntType != "Enhancer-Promoter"),]

#For each TF - Compute score rank
for(i in seq(1,length(utfs))){

  #Select impaired interactions TF
  tmp <- netImpf[which(netImpf$Source == utfs[i]),]
  tmp$log2FC <- log2(tmp$FC)

  #Elements regulatory elements for score computation

  #Unique targets
  ut <- unique(tmp$Target)
  #1/ number of regulatory elements for which the TF binding is impaired
  re <- length(ut)

  #2 number of downstream genes that could be impacted
  ng <- c()
  #TF - promoter (direct)
  ng <- c(ng,unique(tmp$Target[which(tmp$IntType == "TF-Promoter")]))
  #TF - enhancer - promoter(s) (indirect)
  ng <- c(ng,unique(netImp$Target[which(netImp$Source %in% tmp$Target)]))
  #Number of unique downstream genes
  ng <- length(unique(ng))

  #3/ Score weighted MAF
  smaf <- c()
  
  #For each unique regulator impaired
  for(j in seq(1,length(ut))){

    #Regulatory region impaired
    rri <- tmp[which(tmp$Target == ut[j]),]
    #MAF regulator
    mafr <- sum(rri$MAF)

    #For each SNP/binding impairment score
    for(x in seq(1,nrow(rri))){
      tmpsc <- abs(rri$FC[x]) * (rri$MAF[x] * mafr)
      #Add
      smaf <- c(smaf,tmpsc)
      rm(tmpsc)
    }
  }

  #Score
  sc = re * (ng/re) * sum(smaf)

  #Fill in
  finalDf[i,] <- c(utfs[i],sc,NA)
  rm(sc)
}

#Ranks
finalDf <- finalDf[order(finalDf$Score,decreasing=T),]
finalDf$Rank <- seq(1,nrow(finalDf))

#Write
write.table(finalDf,args[3],sep="\t",row.names = F,col.names = T,quote = F)
