#!/bin/bash
while getopts i:p:s:n opt; do
  case ${opt} in
    i ) INPUT_FILE=${OPTARG};;
    p ) PROJECT_NAME=${OPTARG};;
    s ) SPECIES=${OPTARG};;
    n ) PHENO=${OPTARG};;
    \? ) echo "Usage: DirectGeneSearch.sh [-i] INPUT_FILE [-p] PROJECT_NAME [-s] SPECIES [-n] Phenotype file"
      ;;
  esac
done
#INPUT_FILE=$1
#PROJECT_NAME=${2:-"default_project"}
#PHENO=$3
#JOBS=${4:-4}
#SPECIES=$4
#echo 'Which species would your like to run G2g? Sbicolor, Setaria or Zmays'
#read SPECIES
if  [ "${SPECIES}" = "SViridis" ];
then
SPECIES="Setaria";
fi

if [ "${SPECIES}" = "Sbicolor" ] || [ "${SPECIES}" = "Setaria" ] || [ "${SPECIES}" = "Zmays" ]; 
then
echo "Thank you, G2g is now running on ${SPECIES}!!"; 
else
echo "You did not input	the species correctly." && exit 0;
fi


BASE_DIR="/shares/tmockler_share/private/Data/G2g/${SPECIES}"
CURRENT_DIR=$PWD

#echo "${CURRENT_DIR}"
#echo "${INPUT_FILE}"
#echo "${PROJECT_NAME}"
#echo "${BASE_DIR}"


while read -r SNP_START SNP_LOCATION; do
bash ${CURRENT_DIR}/scripts/1.gene_region_finder.sh ${SNP_START} ${SNP_LOCATION} ${PROJECT_NAME} ${SPECIES} & done < ${INPUT_FILE}
wait


Rscript ${CURRENT_DIR}/scripts/2.gene_region_to_gene_lists_${SPECIES}.R $PROJECT_NAME --save
##Clean up outputs
sed -i 's/"//g' ${CURRENT_DIR}/$PROJECT_NAME/working/gene_lists_only_by_region/*.region



##launch gene finders


for region_list in  ${CURRENT_DIR}/$PROJECT_NAME/working/gene_lists_only_by_region/*.region; do
bash ${CURRENT_DIR}/scripts/3.gene_diversity_finder_g2g_${SPECIES}.sh $region_list $PROJECT_NAME $SPECIES&
done
wait




##Clean up

mkdir -p ${CURRENT_DIR}/${PROJECT_NAME}/logs
mv ${CURRENT_DIR}/${PROJECT_NAME}/*log ${CURRENT_DIR}/${PROJECT_NAME}/logs/.

##Generate outputs
ln -s ${BASE_DIR}/resources ${CURRENT_DIR}/${PROJECT_NAME}
cd ${CURRENT_DIR}/${PROJECT_NAME}/results/gene_hits/
Rscript ${CURRENT_DIR}/scripts/4.G2g_create_all_outputs.R 
#Rscript ${CURRENT_DIR}/scripts/Universal_DirectGeneSearch_Routputs_Sbicolor_v2.R
cd ../../..

echo "Run complete."

