#!/usr/bin/env python

import os
import pandas as pd
import argparse
import glob

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input_folder", help="Input folder name")
parser.add_argument("--gene_location_file", help="Gene location file")
parser.add_argument("--region_file", help="Region file")
parser.add_argument("--start_point", type=int, help="Start point")
parser.add_argument("--stop_point", type=int, help="Stop point")
parser.add_argument("--chromosome", help="Chromosome")
parser.add_argument("--output_dir", default=os.getcwd(), help="Output directory (default is current directory)")
parser.add_argument("--bulk_input_file:", help="If you want to supply a bulk set of inputs as a tsv file. Supply the same thing as single input but each arg should be a named column and each row is then a set.")
args = parser.parse_args()

# Load gene location data
gene_location = pd.read_csv(args.gene_location_file, sep="\t")

# Get all region file names
full_path = os.path.join(args.output_dir, args.input_folder)
region_path = os.path.join(full_path, "regions")
region_files = glob.glob(os.path.join(region_path, "*.region"))

# Make output directories
os.makedirs(os.path.join(full_path, "results", "snp_region_full_gene_lists"), exist_ok=True)
os.makedirs(os.path.join(full_path, "working", "gene_lists_only_by_region"), exist_ok=True)

# Creating a new DataFrame to hold data for later
ld_tracking_df = pd.DataFrame(columns=['region', 'LD'])

# Subset gene list based on region boundaries
# Load region of interest

def region_file_parser(file:str, chromosome:str, start:int, stop:int)
    
    gene_location = pd.read_csv(file, sep="\t")

    file_base_name = os.path.splitext(os.path.basename(file))[0]

    # This is only filtering the region file for rows and SNPs of interest.
    region_gene_subset = gene_location[(gene_location['Chrom'] == chromosome) & (gene_location['Ext_Start'] >= start) & (gene_location['Ext_Start'] <= stop)]

    if (len(region_gene_subset) <= 2):

        start = start + 1000
        stop = stop  + 1000

        if (start > min(gene_location['Ext_start']) & stop < max(gene_location['Ext_Start'])):

            region_file_parser(file:str, chromosome:str, start:int, stop:int)
        else:
            print("No genes found in the region.")



for current_region_file in region_files:
    region_of_interest = pd.read_csv(current_region_file, sep="\t", header=None)

    # Set region boundaries
    region_chr = args.chromosome
    region_start = args.start_point
    region_stop = args.stop_point
    region_name = os.path.splitext(os.path.basename(current_region_file))[0]

    # This is only filtering the region file for rows and SNPs of interest.
    region_gene_subset = gene_location[(gene_location['Chrom'] == region_chr) & (gene_location['Ext_Start'] >= region_start) & (gene_location['Ext_Start'] <= region_stop)]

    # Filter for the genes of interest.
    region_gene_hits_only = region_gene_subset['Gene']

    if len(region_gene_subset) < 2:
        ld_tracking_df = ld_tracking_df.append({'region': region_name, 'LD': "No genes in LD, using expanded region"}, ignore_index=True)
    if len(region_gene_subset) >= 2:
        ld_tracking_df = ld_tracking_df.append({'region': region_name, 'LD': "Genes found within LD region"}, ignore_index=True)

    # If there are no genes within region, extend until you find two closest to region
    while len(region_gene_subset) < 2:
        region_start -= 1000
        region_stop += 1000

        # This filters the genomic data for the chromosome of interest and the range of interest as defined by the LD
        region_gene_subset = gene_location[(gene_location['Chrom'] == region_chr) & (gene_location['Ext_Start'] >= region_start) & (gene_location['Ext_Start'] <= region_stop)]

    # Output the tables to file
    output_file_name = f"{region_name}_{region_chr}_{region_start}_{region_stop}.gene_list"
    region_gene_subset.to_csv(os.path.join(full_path, "results", "snp_region_full_gene_lists", output_file_name), sep="\t", index=False)
    region_gene_hits_only.to_csv(os.path.join(full_path, "working", "gene_lists_only_by_region", output_file_name), sep="\t", index=False, header=False)

ld_tracking_df.to_csv(os.path.join(full_path, "results", "LD_tracking_by_region_list.txt"), sep="\t", index=False)
