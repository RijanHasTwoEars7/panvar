#!/bin/bash

##Import user generated list of genes / get number of lines to be tested
GENE_LIST=$1
FOLDER_NAME=${2:-"batch_results"}
SNP_REGION_NAME=${GENE_LIST%_just_genes.region}
OUTPUT_NAME=${SNP_REGION_NAME##*/}
SPECIES=$3

BASE_DIR="/shares/tmockler_share/private/Data/G2g/${SPECIES}"
CURRENT_DIR=$PWD

##Make results directory

mkdir -p ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits
mkdir -p ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME
mkdir -p ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_impactful
mkdir -p ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_full

##Read through list of genes, check if gene exists / if yes, process, if no skip.

while read GENE_ID; do

if grep -Fq "$GENE_ID" ${BASE_DIR}/resources/Zmays_gene_locations_v4/Zmays_v4_gene_loc.txt
then

##Get relevant info for gene of interest

gene_info=$(grep ${GENE_ID} ${BASE_DIR}/resources/Zmays_gene_locations_v4/Zmays_v4_gene_loc.txt)

gene_cord=$(echo ${gene_info} | awk '{print $1":"$2"-"$3}')

chr=$(echo ${gene_info} | awk '{print $1}')

##Extract VCF info using tabix

tabix -h ${BASE_DIR}/resources/snpeff_vcfs/hmp321_agpv4_chr${chr}.snpeff.vcf.gz ${gene_cord} | grep -E "#|${GENE_ID}" > ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${GENE_ID}_temp.vcf

##Convert to single line file

cat ${BASE_DIR}/resources/vcf.header > ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/temp_${GENE_ID}_diversity_output.txt

cat  ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${GENE_ID}_temp.vcf | ${BASE_DIR}/resources/scripts/vcfEffOnePerLine.pl | grep -E "#|${GENE_ID}" | \
java -jar ${BASE_DIR}/resources/scripts/SnpSift.jar extractFields - CHROM POS "ANN[*].FEATUREID" REF ALT "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" "GEN[*].GT" | \
grep Chr >> ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/temp_${GENE_ID}_diversity_output.txt

##Remove multiple/non-relvant gene snps

grep -v "-" ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/temp_${GENE_ID}_diversity_output.txt > ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_full/${GENE_ID}_diversity_output.txt


#Clean up
rm  ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${GENE_ID}_temp.vcf
rm  ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/temp_${GENE_ID}_diversity_output.txt

#Exhange allele scores for numeric values

sed -i 's|0\/0|0|g; s|0\/1|1|g; s|1\/1|2|g; s|.\/.|NA|g' ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_full/${GENE_ID}_diversity_output.txt 

##For filtering only on meaningful variants

grep -E 'CHROM|LOW|MODERATE|HIGH' ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_full/${GENE_ID}_diversity_output.txt > ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_impactful/${GENE_ID}_diversity_output_impactful.txt
 
fi
done < $GENE_LIST

##Merge the generated lists into one full file for each time
cat ${BASE_DIR}/resources/vcf.header > ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${OUTPUT_NAME}_all_gene_hits.txt

for file in ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_full/*.txt; do
grep -v CHROM ${file} >> ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/unsorted_${OUTPUT_NAME}.txt
done
sort -k1,1 -k2,2 >> ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${OUTPUT_NAME}_all_gene_hits.txt < ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/unsorted_${OUTPUT_NAME}.txt
sed -i 's/.v3.1//g' ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${OUTPUT_NAME}_all_gene_hits.txt

#Write file for merged meaningful mutations

grep -E 'CHROM|LOW|MODERATE|HIGH' ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${OUTPUT_NAME}_all_gene_hits.txt > ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${OUTPUT_NAME}_impactful_gene_hits.txt

#Clean up 
rm ${CURRENT_DIR}/${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/unsorted_${OUTPUT_NAME}.txt

