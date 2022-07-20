#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(reshape2)

### Impaired interactions
impNet = fread(args[1],header=T)

# TF-promoter
impNetTfPro <- impNet[which(impNet$IntType == "TF-Promoter"),]

# Init final df
if(ncol(impNetTfPro) == 9){
  affinityNetTfPro <- data.frame("Source"=NA,"SNP"=NA,"Target"=NA,"Ref"=NA,"Alt"=NA,"RSID"=NA,"Motif"=NA,"FC"=NA)
}else{
  affinityNetTfPro <- data.frame("Source"=NA,"SNP"=NA,"Target"=NA,"Ref"=NA,"Alt"=NA,"Motif"=NA,"FC"=NA)
}

print("Promoter regions")
if(nrow(impNetTfPro)>0){
  ugp <- unique(impNetTfPro$Target)
  # For each gene
  for(i in seq(1,length(ugp))){
    
    #TF interactions
    impNetTfPro_spe <- impNetTfPro[which(impNetTfPro$Target == ugp[i]),]
    if(ncol(impNetTfPro_spe) == 9){
      impNetTfPro_spe <- impNetTfPro_spe[,c("Source","Target","SNP","Ref","Alt","RSID")]
      impNetTfPro_spe <- impNetTfPro_spe[!duplicated(impNetTfPro_spe[,1:6]),]
    }else{
      impNetTfPro_spe <- impNetTfPro_spe[,c("Source","Target","SNP","Ref","Alt")]
      impNetTfPro_spe <- impNetTfPro_spe[!duplicated(impNetTfPro_spe[,1:5]),]
    }

    #For each SNP
    snps <- unique(impNetTfPro_spe$SNP)
    for (j in seq(1,length(snps))) {
      #Retrieve fasta sequences for promoter with SNPs
      # 25 dp [SNP] 25 bp
      snpPos <- snps[j]
      snpPosSplit <- c(strsplit(snpPos,":")[[1]][1], #chr
                       strsplit(strsplit(snpPos,":")[[1]][2],"-")[[1]][1], #start
                       strsplit(strsplit(snpPos,":")[[1]][2],"-")[[1]][2]) #stop 
      startReg <- as.numeric(snpPosSplit[2]) - 50
      stopReg <- as.numeric(snpPosSplit[3]) + 50
      pos <- paste(snpPosSplit[1],paste(startReg,stopReg,sep=","),sep=":")
      
      cmd <- c(paste(paste0("curl -s 'http://genome.ucsc.edu/cgi-bin/das/",args[2],"/dna?segment="),paste(pos,paste0("' | xmllint --xpath '/DASDNA/SEQUENCE/DNA/text()' - | tr -d '\n' > ",args[4],"seqSNPreg.fasta"),sep=""),sep=""))
      system(cmd)
      
      if(file.exists(paste0(args[4],"seqSNPreg.fasta"))){
        if(file.info(paste0(args[4],"seqSNPreg.fasta"))$size > 0){
          fasf <- suppressWarnings(read.delim(paste0(args[4],"seqSNPreg.fasta"),sep="\t",header=F))
          #Rm file
          system(paste0("rm ",args[4],"seqSNPreg.fasta"))
          if(nrow(fasf)>0){
            seq <- fasf$V1
            seq <- strsplit(seq,split="")[[1]]
            seq <- toupper(seq)
            seqstart <- paste(seq[1:50],collapse = "")
            seqStop <- paste(seq[52:101],collapse = "")
            seqsnp <- paste("[",paste((impNetTfPro_spe$Ref[j]),paste("/",paste((impNetTfPro_spe$Alt[j]),"]",sep=""),sep=""),sep=""),sep="") #[Ref/Alternate]
            seqfinale <- paste(seqstart,paste(seqsnp,seqStop,sep=""),sep="")
            
            dfsnp <- data.frame("rs"=snps[j],"Seq"=seqfinale)
            #Remove duplicates
            dfsnp <- dfsnp[!duplicated(dfsnp[,1:2]),]
            #Remove NAs
            dfsnp <- dfsnp[!is.na(dfsnp$Seq),]
            if(nrow(dfsnp)>0){
              write.table(dfsnp,paste0(args[4],"SNPinfo.txt"),sep=" ",row.names=F,col.names=F,quote = F)
              
              #Run Java command
              setwd(args[5])
              cmdjava <- paste0("java -cp ape.jar ru.autosome.perfectosape.SNPScan ",args[3]," ",paste0(args[4],"SNPinfo.txt --pvalue-cutoff 0.05 --fold-change-cutoff 2"),"  > ",args[4],"BindingAffinity.txt")
              system(cmdjava)
              system(paste0("rm ",paste0(args[4],"SNPinfo.txt")))
              
              #Process binding affinity info
              bi <- fread(paste0(args[4],"BindingAffinity.txt"),header=T)
              
              if(nrow(bi)>0){
                bi <- bi[,c(1,2,12)]
                colnames(bi) <- c("SNP","Motif","FC")
                bi$Source <- sapply(bi$Motif,function(x){strsplit(x,split="_")[[1]][1]})
                #Filter based on FC cutoff
                #bi <- bi[which(bi$FC>=args[6]),]
                #Filter TF in interactions by merging
                if(nrow(bi)>0){
                  intAff <- merge(impNetTfPro_spe,bi,by=c("Source","SNP"))
                  print(head(intAff))
                  if(ncol(impNetTfPro) == 9){
                    intAff <- intAff[,c("Source","SNP","Target","Ref","Alt","RSID","Motif","FC")]
                  }else{
                    intAff <- intAff[,c("Source","SNP","Target","Ref","Alt","Motif","FC")]
                  }
                  #Bind to final df
                  affinityNetTfPro <- rbind(affinityNetTfPro,intAff)
                  rm(intAff)
                }
              }
              #Clean
              rm(bi)
              system(paste0("rm ",paste0(args[4],"BindingAffinity.txt")))
            }
            rm(seqfinale)
          }
        }
      }
    }
  }
  #Remove NA init
  if(nrow(affinityNetTfPro)>0){
    affinityNetTfPro <- affinityNetTfPro[!is.na(affinityNetTfPro$Source),]
  }
}

# TF-enhancer
print("Enhancer regions")
impNetTFEnh <- impNet[which(impNet$IntType == "TF-Enhancer"),]

# Init final df
if(ncol(impNetTFEnh) == 9){
  affinityNetTFEnh <- data.frame("Source"=NA,"SNP"=NA,"Target"=NA,"Ref"=NA,"Alt"=NA,"RSID"=NA,"Motif"=NA,"FC"=NA)
}else{
  affinityNetTFEnh <- data.frame("Source"=NA,"SNP"=NA,"Target"=NA,"Ref"=NA,"Alt"=NA,"Motif"=NA,"FC"=NA)
}

if(nrow(impNetTFEnh)>0){
  uenh <- unique(impNetTFEnh$Target)

  # For each enhancer
  for(i in seq(1,length(uenh))){
    
    #TF interactions
    impNetTfEnh_spe <- impNetTFEnh[which(impNetTFEnh$Target == uenh[i]),]
    if(ncol(impNetTfEnh_spe) == 9){
      impNetTfEnh_spe <- impNetTfEnh_spe[,c("Source","Target","SNP","Ref","Alt","RSID")]
      impNetTfEnh_spe <- impNetTfEnh_spe[!duplicated(impNetTfEnh_spe[,1:6]),]
    }else{
      impNetTfEnh_spe <- impNetTfEnh_spe[,c("Source","Target","SNP","Ref","Alt")]
      impNetTfEnh_spe <- impNetTfEnh_spe[!duplicated(impNetTfEnh_spe[,1:5]),]
    }
    
    snps <- unique(impNetTfEnh_spe$SNP)
    #For each SNP
    for (j in seq(1,length(snps))) {
      #Retrieve fasta sequences for promoter with SNPs
      # 25 dp [SNP] 25 bp
      snpPos <- snps[j]
      snpPosSplit <- c(strsplit(snpPos,":")[[1]][1], #chr
                       strsplit(strsplit(snpPos,":")[[1]][2],"-")[[1]][1], #start
                       strsplit(strsplit(snpPos,":")[[1]][2],"-")[[1]][2]) #stop 
      startReg <- as.numeric(snpPosSplit[2]) - 50
      stopReg <- as.numeric(snpPosSplit[3]) + 50
      pos <- paste(snpPosSplit[1],paste(startReg,stopReg,sep=","),sep=":")
      
      cmd <- c(paste(paste0("curl -s 'http://genome.ucsc.edu/cgi-bin/das/",args[2],"/dna?segment="),paste(pos,paste0("' | xmllint --xpath '/DASDNA/SEQUENCE/DNA/text()' - | tr -d '\n' > ",args[4],"seqSNPreg.fasta"),sep=""),sep=""))
      system(cmd)
      
      if(file.exists(paste0(args[4],"seqSNPreg.fasta"))){
        if(file.info(paste0(args[4],"seqSNPreg.fasta"))$size > 0){
          fasf <- read.table(paste0(args[4],"seqSNPreg.fasta"),sep="\t",header=F)
          #Rm file
          system(paste0("rm ",args[4],"seqSNPreg.fasta"))
          if(nrow(fasf)>0){
            seq <- fasf$V1
            seq <- strsplit(seq,split="")[[1]]
            seq <- toupper(seq)
            seqstart <- paste(seq[1:50],collapse = "")
            seqStop <- paste(seq[52:101],collapse = "")
            seqsnp <- paste("[",paste((impNetTfEnh_spe$Ref[j]),paste("/",paste((impNetTfEnh_spe$Alt[j]),"]",sep=""),sep=""),sep=""),sep="") #[Ref/Alternate]
            seqfinale <- paste(seqstart,paste(seqsnp,seqStop,sep=""),sep="")
            
            dfsnp <- data.frame("rs"=snps[j],"Seq"=seqfinale)
            dfsnp <- dfsnp[!is.na(dfsnp$Seq),]
            if(nrow(dfsnp)>0){
              write.table(dfsnp,paste0(args[4],"SNPinfo.txt"),sep=" ",row.names=F,col.names=F)
              
              #Run Java command
              cmdjava <- paste0("java -cp ",args[5],"ape.jar ru.autosome.perfectosape.SNPScan ",args[3]," ",paste0(args[4],"SNPinfo.txt --pvalue-cutoff 0.05 --fold-change-cutoff 2")," > ",args[4],"BindingAffinity.txt")
              system(cmdjava)
              system(paste0("rm ",paste0(args[4],"SNPinfo.txt")))
              
              #Process binding affinity info
              bi <- fread(paste0(args[4],"BindingAffinity.txt"),header=T)
              if(nrow(bi)>0){
                bi <- bi[,c(1,2,12)]
                colnames(bi) <- c("SNP","Motif","FC")
                bi$Source <- sapply(bi$Motif,function(x){strsplit(x,split="_")[[1]][1]})
                #Filter based on FC cutoff
                #bi <- bi[which(bi$FC>=args[6]),]
                if(nrow(bi)>0){
                  #Filter TF in interactions by merging
                  intAff <- merge(impNetTfEnh_spe,bi,by=c("Source","SNP"))
                  print(head(intAff))
                  if(ncol(impNetTFEnh) == 9){
                    intAff <- intAff[,c("Source","SNP","Target","Ref","Alt","RSID","Motif","FC")]
                  }else{
                    intAff <- intAff[,c("Source","SNP","Target","Ref","Alt","Motif","FC")]
                  }
                  #Bind to final df
                  affinityNetTFEnh <- rbind(affinityNetTFEnh,intAff)
                  rm(intAff)
                }
              }
              #Clean
              rm(bi)
              system(paste0("rm ",paste0(args[4],"BindingAffinity.txt")))
            }
            rm(seqfinale)
          }
        }
      }
    }
  }
  if(nrow(affinityNetTFEnh)>0){
    #Remove NA init
    affinityNetTFEnh <- affinityNetTFEnh[!is.na(affinityNetTFEnh$Source),]
  }
}

#Bind both analyses
if(nrow(affinityNetTfPro)>0 & nrow(affinityNetTFEnh)>0){
  affinityRes <- rbind(affinityNetTfPro,affinityNetTFEnh)
}else if(nrow(affinityNetTfPro)>0 & nrow(affinityNetTFEnh)==0){
  affinityRes <- affinityNetTfPro
}else if(nrow(affinityNetTfPro)==0 & nrow(affinityNetTFEnh)>0){
  affinityRes <- affinityNetTFEnh
}else{
  #Empty df
  affinityRes <- affinityNetEnhPro
}
affinityRes <- affinityRes[!is.na(affinityRes$Source),]

#TF binding affinity analysis
write.table(affinityRes,args[6],sep="\t",row.names=F,col.names=T,quote=F)

