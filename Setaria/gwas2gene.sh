#!/bin/bash


INPUT_FILE=$1
PHENO=$2
PROJECT_NAME=${3:-"default_project"}
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

for file in $PROJECT_NAME/results/gene_hits/*/*_impactful_gene_hits.txt; do
sed -i -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|1/2/g' $file
done

##Clean up

mkdir -p ${PROJECT_NAME}/logs
mv ${PROJECT_NAME}/*log ${PROJECT_NAME}/logs/.



##Generate outputs
cd  ${PROJECT_NAME}/results/gene_hits/

##Remove empty SNP folders before processing

for file in */*impactful_gene_hits.txt; do
LINECOUNT=`wc -l ${file} | cut -f1 -d' '`

if [[ $LINECOUNT == 1 ]]; then
   echo "${file%%/*} is empty."
        mkdir -p ../no_impactful_hits
        mv "${file%%/*}" ../no_impactful_hits/.
fi
done

##Output script
Rscript --vanilla ../../../resources/scripts/Universal_G2G_Routputs_SViridis_v2.R ${PHENO}
cd ../../..

echo "Run complete."

