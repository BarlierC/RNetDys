#!/bin/sh

#If required export Samtoola & bedtools in your PATH
#PATH=$PATH:/Users/celine.barlier/Desktop/Tools/samtools-1.11/
#PATH=$PATH:/Users/celine.barlier/Desktop/Tools/bedtools2/bin/
#export PATH

CUR_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#Call: sh Method.sh -g network.txt -s SNPs -cp sorted.prom.chip.bed -h HOCOMOCO/pwm/  -o my/output/path/ -n myProject
#SNPs should have 6 columns: chrx, Start, Stop, Ref, Alt alleles, MAF, RS name (optional) but no header !
while getopts g:s:a:o:n: option
do
case "${option}"
in
g) NETWORK=${OPTARG};;
s) SNPS=${OPTARG};;
cp) CHIP_PROMOTER=${OPTARG};;
h) HOCOMOCOPATH=${OPTARG};;
o) OUTPUTPATH=${OPTARG};;
n) PROJECTNAME=${OPTARG};;
esac
done

echo "PART 2: network contextualization towards the disease state"
echo "Performing analysis using network ${NETWORK}"
echo "SNPs to contextualize network ${SNPS}"
echo "ChIP-seq promoter data: ${CHIP_PROMOTER}"
echo "PMW data: ${CHIP_PROMOTER}"
echo "Project name: ${PROJECTNAME}"
echo "Output path: ${OUTPUTPATH}"

###################### First step: mapping SNPs to regions #############################

echo 'Mapping SNPs to network'

#1/ Retrieve promoter coordinates based on assembly (/!\ network & SNPs should be from the same assembly)
echo '1/ Retrieving promoter & enhancers coordinates'
#Promoter & enhancers coordinates
/usr/local/bin/Rscript "$CUR_DIR/"scripts/retrieveCoord.R $NETWORK $CHIP_PROMOTER $OUTPUTPATH$PROJECTNAME".netCoord.bed"
#Sort prior to intersection with SNPs
sort -k1,1 -k2,2n $OUTPUTPATH$PROJECTNAME".netCoord.bed" > $OUTPUTPATH$PROJECTNAME".filtered.netCoord.bed"
#Sort SNPs
sort -k1,1 -k2,2n $SNPS > $SNPS".filtered.bed"
#2/ Map SNPs to regions
echo '2/ Mapping SNPs to regions'
/Users/celine.barlier/Documents/bedtools2/bin/bedtools intersect -a $OUTPUTPATH$PROJECTNAME".filtered.netCoord.bed" -b $SNPS".filtered.bed" -F 1.0 -wo -sorted > $OUTPUTPATH$PROJECTNAME"mappedSNPs.bed"
#3/ Process impaired interactions
echo '3/ Process impaired interactions'
/usr/local/bin/Rscript "$CUR_DIR/"scripts/impairedInt.R $OUTPUTPATH$PROJECTNAME"mappedSNPs.bed" $NETWORK $OUTPUTPATH$PROJECTNAME"_II.txt" $OUTPUTPATH$PROJECTNAME"_impaired_interactions.txt"

############ Output: impaired interactions & SNPs ######################


############ Third step: TF binding impaired & regulator ranking ###############

#Identify impaired TF binding sites
echo '3/ TF binding affinity analysis'
/usr/local/bin/Rscript "$CUR_DIR/"scripts/TFbindingAffinity.R $OUTPUTPATH$PROJECTNAME"_II.txt" $ASSEMBLY $HOCOMOCOPATH $OUTPUTPATH $JAVAFS $OUTPUTPATH$PROJECTNAME"_impaired_TFs.txt"

echo '4/ Refine impaired regulatory interactions'
/usr/local/bin/Rscript "$CUR_DIR/"scripts/bindImpRegInfo.R $OUTPUTPATH$PROJECTNAME"_impaired_interactions.txt" $OUTPUTPATH$PROJECTNAME"_impaired_TFs.txt" $OUTPUTPATH$PROJECTNAME"_impaired_regInt.txt"

echo '5/ TF regulator ranking'
/usr/local/bin/Rscript "$CUR_DIR/"scripts/TFregRanks.R $OUTPUTPATH$PROJECTNAME"_impaired_regInt.txt" $SNPS $OUTPUTPATH$PROJECTNAME"_ranked_TFreg.txt"
############ Output: TF affinity analysis results ######################

#Clean
rm $OUTPUTPATH$PROJECTNAME".netCoord.bed"
rm $OUTPUTPATH$PROJECTNAME".filtered.netCoord.bed"
rm $OUTPUTPATH$PROJECTNAME"mappedSNPs.bed"
rm $OUTPUTPATH$PROJECTNAME"_II.txt"
rm $OUTPUTPATH$PROJECTNAME"_impaired_interactions.txt"
rm $OUTPUTPATH$PROJECTNAME"_impaired_TFs.txt"
rm $SNPS".filtered.bed"

#End
echo 'end'
