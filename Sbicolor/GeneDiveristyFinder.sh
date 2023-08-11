#!/bin/bash

##Import user generated list of genes / get number of lines to be tested
GENE_LIST=$1
FOLDER_NAME=${2:-"batch_results"}
OUTPUT_NAME=${3:-"full_gene_diveristy_output"}
PHENOTYPE_FILE=${4:-""}
NUMBER_OF_LINES=$( wc -l < ${GENE_LIST})
CURRENT_NUMBER=1


##Read through list of genes, check if gene exists / if yes, process, if no skip.

while read GENE_ID; do
echo "Processing ${GENE_ID}, search query ${CURRENT_NUMBER} of ${NUMBER_OF_LINES}."
if grep -Fq "$GENE_ID" resources/BAP_gene_locations/Sbicolor_313_v3_1_gene_locations_extended.txt
then

##Add one to completed total.
CURRENT_NUMBER=$((CURRENT_NUMBER+1))

##Get relevant info for gene of interest

gene_info=$(grep ${GENE_ID} resources/BAP_gene_locations/Sbicolor_313_v3_1_gene_locations_extended.txt)

gene_cord=$(echo ${gene_info} | awk '{print $1":"$2"-"$3}')

chr=$(echo ${gene_info} | awk '{print $1}')

##Extract VCF info using tabix

tabix -h resources/BAP_vcf_files/TERRA_final_${chr}.snpEff.vcf.gz ${gene_cord} | grep -E "#|${GENE_ID}" >  ${GENE_ID}_temp.vcf

##Convert to single line file

cat resources/vcf.header > temp_${GENE_ID}_diversity_output.txt

cat ${GENE_ID}_temp.vcf | resources/scripts/vcfEffOnePerLine.pl | grep -E "#|${GENE_ID}" | \
java -jar resources/scripts/SnpSift.jar extractFields - CHROM POS "ANN[*].FEATUREID" REF ALT "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" "GEN[*].GT" | \
grep Chr >>  temp_${GENE_ID}_diversity_output.txt

##Make results directory

mkdir -p results/${FOLDER_NAME}
mkdir -p results/${FOLDER_NAME}/individual_files_impactful
mkdir -p results/${FOLDER_NAME}/individual_files_full

##Remove multiple/non-relvant gene snps

grep -v "-" temp_${GENE_ID}_diversity_output.txt > results/${FOLDER_NAME}/individual_files_full/${GENE_ID}_diversity_output.txt 

#Clean up
rm ${GENE_ID}_temp.vcf
rm temp_${GENE_ID}_diversity_output.txt

#Exhange allele scores for numeric values

sed -i 's|0\/0|0|g; s|0\/1|1|g; s|1\/1|2|g; s|.\/.|NA|g' results/${FOLDER_NAME}/individual_files_full/${GENE_ID}_diversity_output.txt 

##For filtering only on meaningful variants

grep -E 'CHROM|LOW|MODERATE|HIGH' results/${FOLDER_NAME}/individual_files_full/${GENE_ID}_diversity_output.txt > results/${FOLDER_NAME}/individual_files_impactful/${GENE_ID}_diversity_output_impactful.txt
 
echo "Search completed for ${GENE_ID}."

else

echo "${GENE_ID} not found."
##Add one to completed.
CURRENT_NUMBER=$((CURRENT_NUMBER+1))

fi
done < $GENE_LIST

##Merge the generated lists into one full file for each time
echo "Merging results files to final output files."

cat resources/vcf.header > results/${FOLDER_NAME}/${OUTPUT_NAME}.txt
for file in results/${FOLDER_NAME}/individual_files_full/*.txt; do
grep -v CHROM ${file} >> results/${FOLDER_NAME}/unsorted_${OUTPUT_NAME}.txt
done
sort -k1,1 -k2,2 >> results/${FOLDER_NAME}/${OUTPUT_NAME}.txt < results/${FOLDER_NAME}/unsorted_${OUTPUT_NAME}.txt

#Write file for merged meaningful mutations

grep -E 'CHROM|LOW|MODERATE|HIGH' results/${FOLDER_NAME}/${OUTPUT_NAME}.txt > results/${FOLDER_NAME}/${OUTPUT_NAME}_impactful.txt

#Clean up 
rm results/${FOLDER_NAME}/unsorted_${OUTPUT_NAME}.txt


cd results/${FOLDER_NAME}
Rscript --vanilla /shares/tmockler_share/private/Data/G2g/Sbicolor/resources/scripts/GeneDiversityFinder_Routputs.R $PHENOTYPE_FILE
cd ../..

#Write Out.

echo "Analysis completed."

