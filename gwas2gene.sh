#!/bin/bash
while getopts i:p:s:n: opt; do
  case ${opt} in
    i ) INPUT_FILE=${OPTARG};;
    p ) PROJECT_NAME=${OPTARG};;
    s ) SPECIES=${OPTARG};;
    n ) PHENO=${OPTARG};;
    \? ) echo "Usage: gwas2gene.sh [-i] INPUT_FILE [-p] PROJECT_NAME [-s] SPECIES [-n] Phenotype file"
      ;;
  esac
done

#INPUT_FILE=$1
#PHENO=$2
#PROJECT_NAME=${3:-"default_project"}
#JOBS=${4:-4}

# Comments by Rijan: The following lines simply re-iterate the variables and does nothing else. This is useless.
echo "${SPECIES}"
echo "${PHENO}"
echo "${INPUT_FILE}"
echo "${PROJECT_NAME}"




if  [ "${SPECIES}" = "SViridis" ];
then
SPECIES="Setaria";
fi

if [ "${SPECIES}" = "Sbicolor" ] || [ "${SPECIES}" = "Setaria" ] || [ "${SPECIES}" = "Zmays" ]; 
then
echo "Thank you, G2g is now running on ${SPECIES}!!"; 
else
echo "You did not input the species correctly." && exit 0;
fi


# The following lines make no sense whatsoever from a design perspective. These are hardlinked paths that should not be here.
BASE_DIR="~/panvar/${SPECIES}"
PROGRAM_DIR="~/panvar/universal_species"
CURRENT_DIR=$PWD


while read -r SNP_START SNP_LOCATION; do
bash ${PROGRAM_DIR}/scripts/1.gene_region_finder.sh ${SNP_START} ${SNP_LOCATION} ${PROJECT_NAME} ${SPECIES} &
done < ${INPUT_FILE}
wait

Rscript ${PROGRAM_DIR}/scripts/2.gene_region_to_gene_lists_${SPECIES}.R ${PROJECT_NAME} --save

##Clean up outputs
sed -i 's/"g' ${CURRENT_DIR}/${PROJECT_NAME}/working/gene_lists_only_by_region/*.region # these region files are made by the rscript above.


##launch gene finders


for region_list in  ${CURRENT_DIR}/${PROJECT_NAME}/working/gene_lists_only_by_region/*.region; do
bash ${PROGRAM_DIR}/scripts/3.gene_diversity_finder_g2g_${SPECIES}.sh $region_list ${PROJECT_NAME} ${SPECIES}&
done
wait

for file in ${CURRENT_DIR}/${PROJECT_NAME}/results/gene_hits/*/*_impactful_gene_hits.txt; do
sed -i -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|1/2/g' $file
done

##Clean up

mkdir -p ${CURRENT_DIR}/${PROJECT_NAME}/logs
mv ${CURRENT_DIR}/${PROJECT_NAME}/*log ${CURRENT_DIR}/${PROJECT_NAME}/logs/.


##Generate outputs
cd  ${CURRENT_DIR}/${PROJECT_NAME}/results/gene_hits/

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
ln -s ${BASE_DIR}/resources ${CURRENT_DIR}/${PROJECT_NAME}
echo $PHENO
Rscript --vanilla ${PROGRAM_DIR}/scripts/GWAS2Gene_Pheno_Routputs_${SPECIES}_v3.R ${PHENO} ${CURRENT_DIR} ${PROJECT_NAME}
cd $BASE_DIR

echo "Run complete."

