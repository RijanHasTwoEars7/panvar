#!/bin/bash


INPUT_FILE=$1
PROJECT_NAME=${2:-"default_project"}
JOBS=${3:-4}

bash resources/scripts/gene_region_finder.sh ${INPUT_FILE} ${PROJECT_NAME}

 
Rscript resources/scripts/2.gene_region_to_gene_lists.R -f ${PROJECT_NAME}

echo "Run complete."

