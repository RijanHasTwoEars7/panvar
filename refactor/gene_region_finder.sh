#!/bin/bash

# use the following setup to define longer flags in bash

if command -v tabix >/dev/null; then
    :
else
    echo "This shell does not have acess to tabix. This software requires tabix in path"
    exit
fi

if command -v vcftools >/dev/null; then
    :
else
    echo "This shell does not have acess to vcftools. This software requires vcftools in path"
    exit
fi


while test $# -gt 0; do
  case "$1" in
    --chromosome)
        shift
        chromosome=$1 # the chromosome you want to filter for
        shift
        ;;
    --vcf_file)
        shift
        vcf_file=$1 # this should be the vcf file
        shift
        ;;
    --output)
        shift
        output=$1 # this should be the output path to the output file
        shift
        ;;
    --distance) # this should be the ld distance that is calculated
        shift
        distance=$1
        shift
        ;;
  esac
done

# Extract the base name of the input vcf file for later use 

base_name="$(basename -- $vcf_file)"
base_name="${base_name%.*}" # this should be the base of the file without the extension

# Crunching the numbers for the linkage distance range from the input 
snp_start_ld=$((distance-500000))
snp_stop_ld=$((distance+500000))

## make sure that the start LD is not below zero

if (( ${snp_start_ld} < 1)); then
    snp_start_ld=0
    echo "the given distance was too close to the start and needed to be reset to 0"    
fi

# generate a base name for the output file using the base name of the input file and the ld range
output_base_name="${base_name}_${snp_start_ld}_${snp_stop_ld}.txt"

output_path="${output}/${output_base_name}"

tabix -h ${vcf_file} ${chromosome}:${snp_start_ld}-${snp_stop_ld} > ${output_path}

