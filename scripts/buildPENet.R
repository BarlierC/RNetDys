#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(reshape2)
library(doParallel)
library(stats)
library(propagate)
library(sparseMatrixStats)
library(utils)
library(slam)

### Identify active promoters

#TFGene Network
net = fread(args[1])
expg <- unique(c(net$Source,net$Target))
rm(net)

#ChIP coordinates file
chip = fread(args[2])
#ChIP with expressed genes/TFs
chip <- chip[which(chip[,4] %in% expg | chip[,5] %in% expg),]
#ATAC coordinates that intersected with ChIP
chipatac_net <- data.table("Chr"=chip[,1],"Start"=chip[,2],"Stop"=chip[,3],"Target"=chip[,4])
chipatac_net$coord <- paste(chipatac_net$Chr,paste(chipatac_net$Start,chipatac_net$Stop,sep="-"),sep=":")
#Unique target & coordinates
chipatac_targetcoord <- chipatac_net[,c(4,5)]
chipatac_targetcoord <- chipatac_targetcoord[!duplicated(chipatac_targetcoord[,c(1,2)])]
rm(chip)
rm(expg)

#Active enhancers (with ATAC coordinates that intersected with enhancer regions)
activEnh = fread(args[3])
activEnh$enh <- paste(activEnh[,1],paste(activEnh[,2],activEnh[,3],sep="-"),sep=":")
activEnh$coord <- paste(activEnh[,5],paste(activEnh[,6],activEnh[,7],sep="-"),sep=":")

#scATACseq matrix
mtx_atac = read.table(args[4],header=T,row.names = 1)
if(colnames(mtx_atac)[1] == "X"){
  rownames(mtx_atac) <- mtx_atac$X
  mtx_atac$X <- NULL
}

#Retrieve Enhancer-Promoter interactions

#Correlation matrix
mtx_enh <- mtx_atac[which(rownames(mtx_atac) %in% c(activEnh$coord)),]
mtx_enh <- t(mtx_enh)
mtx_prom <- mtx_atac[which(rownames(mtx_atac) %in% c(chipatac_net$coord)),]
mtx_prom <- t(mtx_prom)

#Open enhancers ATAC peaks
enhvec <- colnames(mtx_enh)
#Open promoters ATAC peaks
promot <- colnames(mtx_prom)

rm(mtx_atac)
rm(chipatac_net)
#Free memory
gc()

#N times M matrix M with N = number of columns in x and M = number of columns in y
#x=promoters (will be splitted in blocks)
#y=enhancers
ffcorMat <- propagate::bigcor(x=mtx_prom,y=mtx_enh,fun="cor",method="pearson")
#Select enhancers in cols and promoters in rows
corMat <- as(ffcorMat[seq(1,ncol(mtx_prom)),seq(1,ncol(mtx_enh))],"dgCMatrix")
rm(mtx_enh)
rm(mtx_prom)
#Free memory
gc()

#Perform statistics to identify high probability active enhancer & significant interactions

#Compute matrix of z-score
meanMtx <-  mean(corMat)
sdMtx <- sd(corMat)
corMatSig <- (corMat - meanMtx) / sdMtx
print(class(corMatSig))
rm(meanMtx)
rm(sdMtx)
gc()

print("Computation of the statistics")

stm <- slam::as.simple_triplet_matrix(x=corMatSig)
#Free memory
rm(corMatSig)
rm(corMat)
#Compute statistics
enhPro <- slam::colapply_simple_triplet_matrix(test2,function(x){
    #Compute pvalues for each enhancer (cols)
    pv <- 2*pnorm(q=x, lower.tail=FALSE)
    #Adjust p-values for each enhancer
    padj <- p.adjust(pv,method="BH")
    #Get significant enhancer-promoter interactions
    id <- which(padj<as.numeric(args[5]))
    return(id)
  })
#Free memory
rm(stm)
gc()
names(enhPro) <- enhvec
enhPro <- enhPro[lapply(enhPro,length)>0]

#Network construction
print("Network construction")

#Unique enhancers
uenh <- activEnh[,c("coord","enh")]
uenh <- uenh[!duplicated(uenh[,1:2]),]
uenh <- uenh[which(uenh$coord %in% names(enhPro)),]
ue <- unique(uenh$enh)

# create progress bar
pb <- txtProgressBar(min = 0, max = length(ue), style = 3)

#For each enhancer 
epNet <- do.call("rbind",lapply(X = 1:length(ue),FUN = function(a){
  
  Sys.sleep(0.1)
  
  #Chromosome enhancer
  chren <- strsplit(ue[a],":")[[1]][1]
  
  # update progress bar
  setTxtProgressBar(pb, a)

  #ATAC peaks enhancer
  atacenh <- uenh$coord[which(uenh$enh == ue[a])]
  #Id coresp enhPro
  ida <- which(names(enhPro) %in% atacenh)
  
  #Promoter names
  if(length(ida)>1){
    pc <- data.frame("id"=enhPro[[ida[1]]],"coord"=promot[enhPro[[ida[1]]]],"ida"=rep(ida[1],length(promot[enhPro[[ida[1]]]])))
    for(x in seq(2,length(ida))){
      pc <- rbind(pc,data.frame("id"=enhPro[[ida[x]]],
                                "coord"=promot[enhPro[[ida[x]]]],
                                "ida"=rep(ida[x],length(enhPro[[ida[x]]]))))
    }
  }else{
    pc <- data.frame("id"=enhPro[[ida]],"coord"=promot[enhPro[[ida]]],"ida"=rep(ida,length(enhPro[[ida]])))
  }
  
  #Keep promoter of the same chromosome
  pc <- pc[str_detect(pc$coord,paste0("^",chren)),]
  
  if(nrow(pc)>0){
    #Get promoter gene names
    prom <- chipatac_targetcoord[which(chipatac_targetcoord$coord %in% pc$coord),]
    pcprom <- merge(pc,prom,by="coord")
    #Net
    if(length(ida)>1){
      net <- data.table("Enhancer"=rep(ue[a],nrow(pcprom[which(pcprom$ida == ida[1]),])),
                        "Promoter"=pcprom$Target[which(pcprom$ida == ida[1])],
                        "CorVal"=ffcorMat[pcprom$id[which(pcprom$ida == ida[1])],ida[1]],
                        "Count"=rep(1,nrow(pcprom[which(pcprom$ida == ida[1]),])))
      for(y in seq(1,length(ida))){
        net <- rbind(net,data.table("Enhancer"=rep(ue[a],nrow(pcprom[which(pcprom$ida == ida[y]),])),
                                    "Promoter"=pcprom$Target[which(pcprom$ida == ida[y])],
                                    "CorVal"=ffcorMat[pcprom$id[which(pcprom$ida == ida[y])],ida[y]],
                                    "Count"=rep(1,nrow(pcprom[which(pcprom$ida == ida[y]),]))))
      }
    }else{
      net <- data.table("Enhancer"=rep(ue[a],nrow(pcprom)),
                        "Promoter"=pcprom$Target,
                        "CorVal"=ffcorMat[pcprom$id,ida],
                        "Count"=rep(1,nrow(pcprom)))
    }
    
    #Compute score: sum of atac peaks correlated with same promoter
    net <- aggregate(net[,c("CorVal","Count")],by=list(net$Enhancer,net$Promoter),sum)
    colnames(net) <- c("Source","Target","CorVal","Count")

  }else{
    net <- data.table("Source"=NA,
                      "Target"=NA,
                      "CorVal"=NA,
                      "Count"=NA)
  }

  return(net)          
}))

close(pb)

#Free memory
rm(chipatac_targetcoord)
rm(enhPro)
rm(ffcorMat)
gc()

#Remove potential NAs
epNet <- epNet[!is.na(epNet$Source),]

#Replace correlation value by actual sign
epNet$Sign <- rep(NA,nrow(epNet))
epNet$Sign[which(epNet$CorVal>0)] <- "+"
epNet$Sign[which(epNet$CorVal<0)] <- "-"

#Intersect with GeneHancer backbone
genehancerf <- fread(args[6])
genehancerf$Source <- paste(genehancerf$Chr,paste(genehancerf$Start,genehancerf$End,sep="-"),sep=":")
genehancerf$Target <- genehancerf$ConnectedGene
genehancerf <- genehancerf[,c("Source","Target")]
genehancerf <- genehancerf[!duplicated(genehancerf[,1:2]),]
epNet <- merge(epNet,genehancerf,by=c("Source","Target"))

#Network file
write.table(epNet,args[7],sep="\t",row.names=F,col.names=T,quote=F)

enh_df <- data.frame("Chr"=rep(NA,nrow(epNet)),"Start"=rep(NA,nrow(epNet)),"Stop"=rep(NA,nrow(epNet)))
enh_df$Chr <- sapply(epNet$Source,function(x){strsplit(x,split=":")[[1]][1]})
enh_df$Start <- sapply(epNet$Source,function(x){strsplit(strsplit(x,split=":")[[1]][2],split="-")[[1]][1]})
enh_df$Stop <- sapply(epNet$Source,function(x){strsplit(strsplit(x,split=":")[[1]][2],split="-")[[1]][2]})
enh_df <- enh_df[!duplicated(enh_df[,c("Chr","Start","Stop")]),]

#Bed file of enhancers
write.table(enh_df,args[8],sep="\t",row.names=F,col.names=T,quote=F)
