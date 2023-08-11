#!/bin/bash
FILE_NAME=$1
PROJECT_NAME=$2
REGION_NAME_LABEL=${REGION_NAME##regions/}


FILE="${PROJECT_NAME}/gene_locations_with_regions/${FILE_NAME%.region}_gene_hits.region"


##If an old file exists, remove it
if [ -f $FILE ] ; then
    rm $FILE
fi


##Create output directory
mkdir -p ${PROJECT_NAME}/gene_locations_with_regions

#while read -r CHR START STOP; do
#	while read -r GENE_CHR GENE_START GENE_STOP GENE_ID; do 
#		if [[ "$CHR" == "$GENE_CHR"  ]] && (( $GENE_START >= $START )) && (($GENE_START <= $STOP )); then
#			echo "${REGION_NAME_LABEL%.region}	$GENE_ID	$GENE_CHR	$GENE_START	$GENE_STOP"	>> ${PROJECT_NAME}/gene_locations_with_regions/${REGION_NAME_LABEL%.region}_gene_hits.region
#		fi
#	done < resources/BAP_gene_locations/Sbicolor_313_v3_1_gene_locations_extended.txt

#if [[ $(wc -l ${PROJECT_NAME}/gene_locations_with_regions/${REGION_NAME_LABEL%.region}_gene_hits.region) = 0 ]]; then
#	rm ${PROJECT_NAME}/gene_locations_with_regions/${REGION_NAME_LABEL%.region}_gene_hits.region
#fi

#done < $REGION_NAME
