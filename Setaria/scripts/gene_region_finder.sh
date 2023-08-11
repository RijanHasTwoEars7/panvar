#!/bin/bash
INPUT_FILE=$1
PROJECT_FOLDER=${2:-"default_project"}
LD_VALUE=${3:-.01}


##Make folders
mkdir -p ${PROJECT_FOLDER}/
mkdir -p ${PROJECT_FOLDER}/regions

##loop through regions
while read -r SNP_CHR SNP_LOCATION; do

echo "Chr	Snp" > ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp
echo "$SNP_CHR	$SNP_LOCATION" >> ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp

#creat temp VCF that is large window around SNP target correct for less than 1


##Create large window based of SNP position
SNP_TEMP_REGION_START=$((SNP_LOCATION-500000))
SNP_TEMP_REGION_STOP=$((SNP_LOCATION+500000))

if (( ${SNP_TEMP_REGION_START} < 1 )); then
	SNP_TEMP_REGION_START=1
fi

##CREATE A TEMPORARY VCF TO LOOK FOR LOCAL LD
tabix -h resources/BAP_vcf_files/TERRA_final_${SNP_CHR}.snpEff.vcf.gz ${SNP_CHR}:${SNP_TEMP_REGION_START}-${SNP_TEMP_REGION_STOP} > ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_temp_region.vcf

##CALCULATE LD FOR REGION
vcftools --vcf ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_temp_region.vcf --geno-r2-positions ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp --ld-window 500 --out ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD

##Remove SNPs that have LD > LD_VALUE, remove errored SNPs, remove header, Sort based on position and output to temp file

awk -v LD_VAR="$LD_VALUE" '{ if ( $6 >= LD_VAR ) {print} }' ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld | grep -v "nan" | grep -v 'N_INDV' | sort  -k4,4n > ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld.sorted.temp

REGION_CHROM=$( awk 'NR==1{print $3}' ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld.sorted.temp )
REGION_START=$( awk 'NR==1{print $4}' ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld.sorted.temp )
REGION_STOP=$( tail ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld.sorted.temp | awk 'NR==1{print $4}' - )


echo "$REGION_CHROM $REGION_START $REGION_STOP" > ${PROJECT_FOLDER}/regions/${SNP_CHR}_${SNP_LOCATION}.region

#CLEAN UP
rm ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_temp_region.vcf
rm ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp
rm ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld.sorted.temp
rm ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_SNP_LD.list.geno.ld

mkdir -p ${PROJECT_FOLDER}/logs
mv ${PROJECT_FOLDER}/*.log ${PROJECT_FOLDER}/logs/.

done <$INPUT_FILE
