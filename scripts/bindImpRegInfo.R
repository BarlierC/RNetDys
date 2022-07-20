#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stringr)

#General impairment
netImp <- fread(args[1],header=T)

#Impairment with TF binding information
tfImp <- fread(args[2],header=T)

#Format final data.frame

#TF-enhancers
tfEnhImp <- tfImp[str_detect(tfImp$Target,"chr"),]
if(nrow(tfEnhImp)>0){
  tfEnhImp$IntType <- "TF-Enhancer"
  tfEnhImp$RSID <- NULL
  tfEnhImp <- merge(tfEnhImp,netImp,by=c("Source","Target","IntType","SNP","Ref","Alt"))
  tfEnhImp <- tfEnhImp[,c("Source","Target","IntType","Sign","SNP","Ref","Alt","FC","RSID")]
  
  #Retrieve gene regulated by enhancer impaired with SNP & TF binding impairment
  enh <- unique(tfImp$Target[str_detect(tfImp$Target,"chr")])
  enhNetImp <- netImp[which(netImp$IntType == "Enhancer-Promoter" & netImp$Source %in% enh),]
  enhNetImp <- enhNetImp[!duplicated(enhNetImp[,c("Source","Target","IntType","Sign")]),]
  enhNetImp$FC <- NA
  enhNetImp$SNP <- NA
  enhNetImp$Ref <- NA
  enhNetImp$Alt <- NA
  enhNetImp$RSID <- NA
  enhNetImp <- enhNetImp[,c("Source","Target","IntType","Sign","SNP","Ref","Alt","FC","RSID")]
}
#TF-genes
tfGenesImp <- tfImp[str_detect(tfImp$Target,"^[A-Z]"),]
if(nrow(tfGenesImp)>0){
  tfGenesImp$IntType <- "TF-Promoter"
  tfGenesImp$RSID <- NULL
  tfGenesImp <- merge(tfGenesImp,netImp,by=c("Source","Target","IntType","SNP","Ref","Alt"))
  tfGenesImp <- tfGenesImp[,c("Source","Target","IntType","Sign","SNP","Ref","Alt","FC","RSID")]
}

#Bind
if(nrow(tfGenesImp)>0 & nrow(tfEnhImp)>0){
  finalDf <- rbind(tfEnhImp,tfGenesImp)
  finalDf <- rbind(finalDf,enhNetImp)
}else if(nrow(tfGenesImp)==0 & nrow(tfEnhImp)>0){
  finalDf <- rbind(tfEnhImp,enhNetImp)
}else{
  finalDf <- tfGenesImp
}

#Write
write.table(finalDf,args[3],sep="\t",row.names = F,col.names = T,quote = F)
