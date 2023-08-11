#!/bin/bash
SNP_CHR=$1
SNP_LOCATION=$2
PROJECT_FOLDER=${3:-"default_project"}

##Make folders
mkdir -p ${PROJECT_FOLDER}/
mkdir -p ${PROJECT_FOLDER}/regions

##Create region files through regions
echo "Chr	Snp" > ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp
echo "$SNP_CHR	$SNP_LOCATION" >> ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp


##Create large window based of SNP position
REGION_START=$((SNP_LOCATION-30000))
REGION_STOP=$((SNP_LOCATION+30000))

if (( ${REGION_START} < 1 )); then
	REGION_START=1
fi

REGION_CHROM=${SNP_CHR}

echo "$REGION_CHROM $REGION_START $REGION_STOP" > ${PROJECT_FOLDER}/regions/${SNP_CHR}_${SNP_LOCATION}.region

#CLEAN UP
rm ${PROJECT_FOLDER}/${SNP_CHR}_${SNP_LOCATION}_snp_of_interest.temp
