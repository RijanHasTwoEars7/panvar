#!/bin/bash


INPUT_FILE=$1
PROJECT_NAME=${2:-"default_project"}
PHENO=$3
#JOBS=${4:-4}


while read -r SNP_START SNP_LOCATION; do
bash resources/scripts/1.gene_region_finder.sh ${SNP_START} ${SNP_LOCATION} ${PROJECT_NAME} &
done < ${INPUT_FILE}
wait

Rscript resources/scripts/2.gene_region_to_gene_lists.R $PROJECT_NAME --save
##Clean up outputs
sed -i 's/"//g' $PROJECT_NAME/working/gene_lists_only_by_region/*.region


##launch gene finders


for region_list in  $PROJECT_NAME/working/gene_lists_only_by_region/*.region; do
bash resources/scripts/3.gene_diversity_finder_g2g.sh $region_list $PROJECT_NAME &
done
wait

##Clean up

mkdir -p ${PROJECT_NAME}/logs
mv ${PROJECT_NAME}/*log ${PROJECT_NAME}/logs/.

##Generate outputs
cd  ${PROJECT_NAME}/results/gene_hits/
Rscript ../../../resources/scripts/4.G2g_create_all_outputs.R
cd ../../..

echo "Run complete."

