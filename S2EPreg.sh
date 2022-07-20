#!/bin/sh

#If required export Samtoola & bedtools in your PATH
#PATH=$PATH:/Users/celine.barlier/Desktop/Tools/samtools-1.11/
#PATH=$PATH:/Users/celine.barlier/Desktop/Tools/bedtools2/bin/
#export PATH

CUR_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#Call: sh Method.sh -r RNAmat.txt -a ATACmat.txt -c sorted.all.chip.bed -cp sorted.prom.chip.bed -ep genehancer.bed -o my/output/path/ -n myProject -p 0.05
#scRNAseq & scATACseq needs to be QC before using this pipeline !
while getopts r:a:s:o:n:c: option
do
case "${option}"
in
r) SCRNASEQ=${OPTARG};;
a) SCATACSEQ=${OPTARG};;
c) CHIP_ALL=${OPTARG};;
cp) CHIP_PROMOTER=${OPTARG};;
ep) GENEHANCER=${OPTARG};;
o) OUTPUTPATH=${OPTARG};;
n) PROJECTNAME=${OPTARG};;
p) CUTOFF=${OPTARG};;
esac
done

echo "Performing analysis using assembly ${ASSEMBLY}"
echo "Project name: ${PROJECTNAME}"
echo "Output path: ${OUTPUTPATH}"
echo "scRNAseq data: ${SCRNASEQ}"
echo "scATACseq data: ${SCATACSEQ}"
echo "ChIP-seq data: ${CHIP_ALL}"
echo "ChIP-seq promoter data: ${CHIP_PROMOTER}"
echo "Enhancer-Promoter prior-knowledge: ${GENEHANCER}"
echo "Cutoff: ${CUTOFF}"

###################### First step: inference of TF-Gene interactions ######################

echo 'Starting building of the regulatory interactions'

#1/ Extract bed coordinates from scATACseq matrix & sort bed file
echo '1/ Processing scATACseq'
Rscript "$CUR_DIR/"scripts/processATACseq.R $SCATACSEQ $OUTPUTPATH$PROJECTNAME".scATACseq.bed"
sort -k1,1 -k2,2n $OUTPUTPATH$PROJECTNAME".scATACseq.bed" > $OUTPUTPATH$PROJECTNAME".filtered.sorted.bed"

#2/ intersect with CHIP_PROMOTER (aspecific & containing only promoter regions e.g. [-1500,500])
echo '2/ Identifying open binding regions'
bedtools intersect -a $OUTPUTPATH$PROJECTNAME".filtered.sorted.bed" -b $CHIP_PROMOTER -f 0.48 -r -wo -sorted > $OUTPUTPATH$PROJECTNAME"_ChIP_ATAC_backbone.bed"

#3/ identify gene expressed with scRNAseq & filter out any gene not expressed in the network
echo '3/ TF-Gene regulatory interaction inference'
Rscript "$CUR_DIR/"scripts/getExpressedGenes.R $SCRNASEQ $OUTPUTPATH$PROJECTNAME"_ChIP_ATAC_backbone.bed" $OUTPUTPATH$PROJECTNAME"_TFGeneNetToSign.txt"

#4/ correlation of scRNAseq to sign TF-Gene interactions
echo '4/ Sign TF-Gene regulatory interactions'
Rscript "$CUR_DIR/"scripts/getSignTFGeneInts.R $SCRNASEQ $OUTPUTPATH$PROJECTNAME"_TFGeneNetToSign.txt" $OUTPUTPATH$PROJECTNAME"_TFGeneNet.txt"
##################### Output: TF-Gene network ##################################


############ Second step: identify enhancer - promoter interactions ###############

#1/ Identify open enhancer regions: intersect enhancers from GENEHANCER with ATACseq: 100% of ATAC peak fall inside an enhancer region
echo '4/ Identification of active enhancer regions'
bedtools intersect -a $GENEHANCER -b $OUTPUTPATH$PROJECTNAME".filtered.sorted.bed" -F 1.0 -wo -sorted > $OUTPUTPATH$PROJECTNAME"_activeEnhancers.bed"

#2/ Peaks correlation
echo '5/ Enhancer-Promoter regulatory interaction inference'
Rscript "$CUR_DIR/"scripts/buildPENet.R $OUTPUTPATH$PROJECTNAME"_TFGeneNet.txt" $OUTPUTPATH$PROJECTNAME"_ChIP_ATAC_backbone.bed" $OUTPUTPATH$PROJECTNAME"_activeEnhancers.bed" $SCATACSEQ $CUTOFF $GENEHANCER $OUTPUTPATH$PROJECTNAME"_PENet.txt" $OUTPUTPATH$PROJECTNAME"_EP.bed"

#3/ Process PENet file for further analyses
sort -k1,1 -k2,2n $OUTPUTPATH$PROJECTNAME"_EP.bed" > $OUTPUTPATH$PROJECTNAME"sorted.EP.bed"

echo '5/ TF-Enhancer regulatory interaction inference'
#4/ identify TF binding to enhancers identified in 2/ : 100% of TF peak fall within Open Enhancer region.
bedtools intersect -a $OUTPUTPATH$PROJECTNAME"sorted.EP.bed" -b $CHIP_ALL -F 1.0 -wo -sorted > $OUTPUTPATH$PROJECTNAME"TF_EP.bed"

#5/ process final TF-Enhancer-Promoter network 
Rscript "$CUR_DIR/"scripts/getTFPEnet.R $OUTPUTPATH$PROJECTNAME"TF_EP.bed" $OUTPUTPATH$PROJECTNAME"_TFGeneNet.txt" $OUTPUTPATH$PROJECTNAME"TFEPNet.txt"
############ Output: Enhancer - Promoter network ###################


######## Third step: combine all the inferred subnetworks of the two steps ########
echo 'Finalizing regulatory network'
#1/remove enhancer from enhancer-promoter if no TF is binding to it
Rscript "$CUR_DIR/"scripts/getFinalNet.R $OUTPUTPATH$PROJECTNAME"_TFGeneNet.txt" $OUTPUTPATH$PROJECTNAME"_PENet.txt" $OUTPUTPATH$PROJECTNAME"TFEPNet.txt" $OUTPUTPATH$PROJECTNAME"_S2PEreg_net_toSign.txt"

#2/ sign the remaining interactions of the network
Rscript "$CUR_DIR/"scripts/signFinalNet.R $SCRNASEQ  $OUTPUTPATH$PROJECTNAME"_S2PEreg_net_toSign.txt" $OUTPUTPATH$PROJECTNAME"_S2PEreg_net.txt"

#Clean folder
rm $OUTPUTPATH$PROJECTNAME".scATACseq.bed"
rm $OUTPUTPATH$PROJECTNAME"_ChIP_ATAC_backbone.bed"
rm $OUTPUTPATH$PROJECTNAME"_activeEnhancers.bed"
rm $OUTPUTPATH$PROJECTNAME".filtered.sorted.bed"
rm $OUTPUTPATH$PROJECTNAME"_PENet.txt"
rm $OUTPUTPATH$PROJECTNAME"_EP.bed"
rm $OUTPUTPATH$PROJECTNAME"TF_EP.bed"
rm $OUTPUTPATH$PROJECTNAME"TFEPNet.txt"
rm $OUTPUTPATH$PROJECTNAME"_TFGeneNet.txt"
rm $OUTPUTPATH$PROJECTNAME"sorted.EP.bed"
rm $OUTPUTPATH$PROJECTNAME"_TFGeneNetToSign.txt"
rm $OUTPUTPATH$PROJECTNAME"_S2PEreg_net_toSign.txt"
echo 'End'
####### FINAL OUTPUT: Enhancer-Promoter regulatory interaction specific to the subpopulation given as an input ########
