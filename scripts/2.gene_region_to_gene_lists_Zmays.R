#!/usr/bin/env Rscript

##LOAD DEPENDENCIES
if (!require(optparse)) install.packages('optparse')
if (!require(tidyr)) install.packages('tidyr')

require(optparse)

##LOAD IN ARGUMENTS FOR FILES
args <- commandArgs()
input_folder <- args[6]

##Get Path
species_path <- setwd("/shares/tmockler_share/private/Data/G2g/Zmays/")

## DEPENDENCIES

require(tidyr)

##Load Gene Location Data
gene_location <- read.table(file = paste(species_path, 
                                         "resources/Zmays_gene_locations_v4/Zmays_v4_gene_loc.txt", 
                                         sep="/"), header = T, sep = "\t")

wd_path <- setwd("/shares/tmockler_share/private/Data/G2g/universal_species")

##get all region file names
region_files <- list.files(path = paste(wd_path, "/", input_folder, "/regions/", sep = ""), pattern = ".region")

full_path <- paste(wd_path, input_folder, sep="/")
region_path <- paste(full_path, "regions/", sep="/")

##Make output directories

dir.create(file.path(full_path, "results"), showWarnings = FALSE)
dir.create(file.path(full_path, "results/snp_region_full_gene_lists"), showWarnings = FALSE)
dir.create(file.path(full_path, "working"), showWarnings = FALSE)
dir.create(file.path(full_path, "working/gene_lists_only_by_region"), showWarnings = FALSE)


ld_tracking_df <- data.frame(region=character(),
                 LD=character(),  
                 stringsAsFactors=FALSE) 

##Subset gene list based on region boundries
##Load region of interest

for( current_region_file in region_files ){

region_of_interest <- read.table(file=paste(region_path, current_region_file,
                                            sep = ""), header = FALSE)

#Set region boundries
region_chr <- region_of_interest$V1
region_start <- region_of_interest$V2
region_stop <- region_of_interest$V3
region_name <- gsub(x = current_region_file, pattern = ".region", replacement = "")

region_gene_subset <- gene_location %>% subset(Chrom == as.character(region_chr) & Ext_Start >= region_start &
                                                 Ext_Start <= region_stop)

region_gene_hits_only <- region_gene_subset$Gene

##annotate LD tracking DF

if(nrow(region_gene_subset) < 2){
ld_tracking_df_temp <- as.data.frame(rbind(c(region=region_name, LD=paste("No genes in LD, using expanded region"))))
ld_tracking_df <- rbind(ld_tracking_df, ld_tracking_df_temp)
}
if(nrow(region_gene_subset) >= 2){
ld_tracking_df_temp <- as.data.frame(rbind(c(region=region_name, LD=paste("Genes found within LD region"))))
ld_tracking_df <- rbind(ld_tracking_df, ld_tracking_df_temp)
}

#if there are no genes within region, extend until you find two closest to region
while(nrow(region_gene_subset) < 2){
region_start <- region_start - 1000
region_stop <- region_stop + 1000

region_gene_subset <- gene_location %>% subset(Chrom == as.character(region_chr) & Ext_Start >= region_start &
                                                 Ext_Start <= region_stop)
}


write.table(x = region_gene_subset, file = paste(full_path, "/results/", "snp_region_full_gene_lists/", region_name,
                                                 "_all_genes.region", sep=""), row.names = F, sep = "\t")
write.table(x = region_gene_hits_only, file = paste(full_path, "/working/", "gene_lists_only_by_region/",  
                                                    region_name, "_just_genes.region", sep=""), row.names = F, col.names = F, sep = "\t")                                                     
}

write.table(x = ld_tracking_df, file = paste(full_path, "/results/LD_tracking_by_region_list.txt", sep=""), row.names = F, sep = "\t")

