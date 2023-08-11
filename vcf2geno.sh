#!/bin/bash

SNP_CHR=$1
SNP_LOCATION=$2

SPECIES="Sbicolor"
BASE_DIR="/shares/tmockler_share/private/Data/G2g/${SPECIES}"
CURRENT_DIR=$PWD


SNP_TEMP_REGION_START=$((SNP_LOCATION-5000))
SNP_TEMP_REGION_STOP=$((SNP_LOCATION+5000))

if (( ${SNP_TEMP_REGION_START} < 1 )); then
        SNP_TEMP_REGION_START=1

fi

gene_cord=$SNP_CHR":"$SNP_TEMP_REGION_START"-"$SNP_TEMP_REGION_STOP

tabix -h ${BASE_DIR}/resources/BAP_vcf_files/TERRA_final_${SNP_CHR}.snpEff.vcf.gz ${gene_cord} > ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_region.vcf


grep ^Chr ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_region.vcf | ${BASE_DIR}/resources/scripts/vcfEffOnePerLine.pl | java -jar ${BASE_DIR}/resources/scripts/SnpSift.jar extractFields - CHROM POS "ANN[*].FEATUREID" REF ALT "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" "GEN[*].GT" >> ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_geno.vcf


sed -i 's|0\/0|0|g; s|0\/1|1|g; s|1\/1|2|g; s|.\/.|NA|g' ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_geno.vcf


grep ${SNP_LOCATION} ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_geno.vcf | head -1 > ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_2_geno.vcf

#NA=grep -wc "NA" ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_2_geno.vcf
NAs=`awk -F 'NA' '{print NF-1}' ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_2_geno.vcf`

#if ((awk 'NR==1{print NF}' ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno.vcf) = "371") then cut -f 9- ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno.vcf >  ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno_only.vcf

cut -f 9- ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_2_geno.vcf >  ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_geno_only.vcf


REF=`awk -F '0' '{print NF-1}' ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_geno_only.vcf`
ALT=`awk -F '2' '{print NF-1}' ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_geno_only.vcf`



echo -e "SNP\tNA\tREF\tALT\n${SNP_CHR}_${SNP_LOCATION}\t${NAs}\t${REF}\t${ALT}" > ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno.count


#echo -e "${NAs}\t${REF}\t${ALT}" >> ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno.count
#printf '%s,%s,%s\n' "$NAs" "$REF" "$ALT"  >> ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno.countq


cat ${BASE_DIR}/resources/vcf.header ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp_2_geno.vcf > ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_geno.vcf

rm -rf ${CURRENT_DIR}/${SNP_CHR}_${SNP_LOCATION}_temp*.vcf
