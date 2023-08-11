#!/bin/bash


INPUT_FILE=$1
PROJECT_NAME=${2:-"default_project"}
PHENO=$3
RESULTS_FOLDER=${4:-"results_folder"}
JOBS=${5:-4}


while read -r SNP_START SNP_LOCATION; do
bash resources/scripts/1.gene_region_finder.sh ${SNP_START} ${SNP_LOCATION} ${PROJECT_NAME} .001 &
done < ${INPUT_FILE}
wait

Rscript resources/scripts/2.gene_region_to_gene_lists.R $PROJECT_NAME --save
##Clean up outputs
sed -i 's/"//g' $PROJECT_NAME/working/gene_lists_only_by_region/*.region
sed -i 's/"//g' $PROJECT_NAME/results/snp_region_full_gene_lists/*.region

##launch gene finders


for region_list in  $PROJECT_NAME/working/gene_lists_only_by_region/*.region; do
bash resources/scripts/3.gene_diversity_finder_g2g.sh $region_list $PROJECT_NAME &
done
wait

##Generate outputs
cd  ${PROJECT_NAME}/results/gene_hits/
Rscript ../../../resources/scripts/Universal_G2G_Routputs_Zmays.R $PHENO
cd ../../..
 
##Clean up

mkdir -p ${PROJECT_NAME}/logs
mv ${PROJECT_NAME}/*log ${PROJECT_NAME}/logs/.

mkdir -p $RESULTS_FOLDER
mv ${PROJECT_NAME} ${RESULTS_FOLDER}/.

echo "Run complete."

