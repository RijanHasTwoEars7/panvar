#!/bin/bash

#import user generated list of genes / get number of lines to be tested
GENE_LIST=$1
FOLDER_NAME=${2:-"batch_results"}
SNP_REGION_NAME=${GENE_LIST%_just_genes.region}
OUTPUT_NAME=${SNP_REGION_NAME##*/}

##Make results directory

mkdir -p ${FOLDER_NAME}/results/gene_hits
mkdir -p ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME
mkdir -p ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_impactful
mkdir -p ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/individual_files_full

##Read through list of genes, check if gene exists / if yes, process, if no skip.

while read GENE_ID; do

if grep -Fq "$GENE_ID" resources/Zmays_gene_locations_v4/Zmays_v4_gene_loc.txt
then

##Get relevant info for gene of interest

gene_info=$(grep ${GENE_ID} resources/Zmays_gene_locations_v4/Zmays_v4_gene_loc.txt)

gene_cord=$(echo ${gene_info} | awk '{print $1":"$2"-"$3}')

chr=$(echo ${gene_info} | awk '{print $1}')

##Extract VCF info using tabix

tabix -h resources/snpeff_vcfs/hmp321_agpv4_chr${chr}.snpeff.vcf.gz ${gene_cord} | grep -E "#|${GENE_ID}" > ${FOLDER_NAME}/results/gene_hits/$OUTPUT_NAME/${GENE_ID}_temp.vcf

fi
done < $GENE_LIST
