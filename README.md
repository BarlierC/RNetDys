# RNetDys

:warning: **This repository contains the version of the method presented in my PhD thesis (University of Luxembourg): https://orbilu.uni.lu/handle/10993/52240**

:white_check_mark: **For the official and latest version of RNetDys please go to: https://gitlab.lcsb.uni.lu/CBG/RNetDys**

RNetDys is a systematic pipeline to decipher cell (sub)type specific regulatory interactions impaired by disease-related SNPs. This pipeline relies on multi-OMICS data to provide mechanistic insights into impaired regulatory interactions due to SNPs by leveraging the information obtained from the GRN inference. 

RNetDys contains two main parts composed of (1) the cell (sub)type specific GRN inference and (2) the capture of impaired regulatory interactions due to disease-related SNPs to gain regulatory mechanistic insights for the disease.

## Guidelines

### Content
The two parts implemented in RNetDys are accessible with two different scripts:

-	S2EPReg.sh for the GRN inference (part 1)
-	Contextualize.sh for the identification of impaired interactions (Part 2)

Please note that RNetDys requires samtools, bedtools, R (implemented with version 4) and Java.

-	Samtools: http://www.htslib.org/download/
-	BedTools: https://bedtools.readthedocs.io/en/latest/content/installation.html
-	R: https://www.r-project.org
-	Java: https://www.java.com/fr/download/help/download_options_fr.html

### Prior-knowledge data

RNetDys relies on ChIP-seq and enhancer-promoter prior-knowledge data:

- The ChIP-seq data should contain five columns: Chr, Start, Stop, Target (Promoter/Gene) - Source (TF that binds)

*ChIP-seq data can be downloaded and processed from the ChIP-ATLAS: https://www.chip-atlas.org or Cistrome database: http://cistrome.org*

**The ChIP-seq files provided in RNetDys needs to be sorted**
```console
cb:RNetDys@celine.barlier$ sort -V -k1,1 -k2,2 chipfile.bed
```

- The enhancer-promoter data should contain four columns: EnhChr, EnhStart, EnhStop, connected gene

*GeneHancer data (Human) can be requested at: https://www.genecards.org/Guide/GeneCard*

- Hocomoco PWM matrices can be downloaded here: https://hocomoco11.autosome.org/downloads_v11

### Preparing R session

The following R packages are required to run RNetDys:

```r
library(data.table)
library(stringr)
library(reshape2)
library(doParallel)
library(stats)
library(propagate)
library(sparseMatrixStats)
library(utils)
library(slam)
```

### Part 1: Cell (sub)type specific GRN inference

#### Parameters

-r scRNA-seq matrix file (genes in rows, cells in columns)
-a scATAC-seq matrix file (peaks in rows, cells in columns)
-c sorted ChIP-seq file (bed format, no header, sorted)
-cp sorted ChIP-seq file containing only promoter regions (bed format, no header, sorted)
-ep enhancer-promoter prior-knowledge data file (e.g., GeneHancer - bed format, no header)
-o output path 
-n name of the project 
-p padj cutoff for the peak correlation matrix

#### Command line

```console
cb:RNetDys@celine.barlier$ RNETDYSPATH="/Users/celine.barlier/Desktop/RNetDys/"
cb:RNetDys@celine.barlier$ OUTPUTPATH="/Users/celine.barlier/Desktop/"
cb:RNetDys@celine.barlier$ sh $RNETDYSPATH"S2EPreg.sh" -r "scRNAseq_mtx.txt" -a "scATACseq_mtx.txt" -c "sorted.all.chip.bed" -cp "sorted.prom.chip.bed" -ep "genehancer.bed" -o $OUTPUTPATH -n "Astro" -p 0.05
```

### Part 2: Contextualization towards the disease state

#### Parameters

-g GRN file 
-s SNPs file (should contain at least 6 columns: Chr, Start, Stop, Ref, Alt and MAF score; one additional column with RSID is optional) 
-cp sorted ChIP-seq file containing only promoter regions
-h HOCOMOCO PWM matrices path 
-o output path 
-n name of the project 

#### Command line

```console
cb:RNetDys@celine.barlier$ sh $RNETDYSPATH"Contextualize.sh" -g "Astro_S2PEreg_net.txt" -s "AD_SNPs.txt" -cp " sorted.prom.chip.bed" -h “~/Desktop/HOCOMOCO/PWM/” -o $OUTPUTPATH -n "AD_Astro"
```


