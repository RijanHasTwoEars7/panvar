#!/usr/bin/env Rscript

library(data.table)

##LOAD DEPENDENCIES
if (!require(optparse)) install.packages('optparse')
if (!require(tidyr)) install.packages('tidyr')

require(optparse)

##LOAD IN ARGUMENTS FOR FILES
args <- commandArgs()
input_folder <- args[6] ## comments by Rijan: The sixth argument is the name of the input folder. I am not sure where the SIXTH argument is coming from...

##Get Path
wd_path <- getwd()

## DEPENDENCIES

require(tidyr)

##Load Gene Location Data
gene_location = fread("~/panvar/Sbicolor/resources/BAP_gene_locations/Sbicolor_313_v3_1_gene_locations_extended.txt")

##get all region file names
region_files <- list.files(path = paste(wd_path, "/", input_folder, "/regions/", sep = ""), pattern = ".region")

full_path <- paste(wd_path, input_folder, sep="/")
region_path <- paste(full_path, "regions/", sep="/")

##Make output directories

dir.create(file.path(full_path, "results"), showWarnings = FALSE)
dir.create(file.path(full_path, "results/snp_region_full_gene_lists"), showWarnings = FALSE)
dir.create(file.path(full_path, "working"), showWarnings = FALSE)
dir.create(file.path(full_path, "working/gene_lists_only_by_region"), showWarnings = FALSE)

# creating a new data frame to hold data for later
ld_tracking_df <- data.frame(region=character(),
                 LD=character(),  
                 stringsAsFactors=FALSE) 

##Subset gene list based on region boundries
##Load region of interest

for( current_region_file in region_files ){

region_of_interest <- read.table(file=paste(region_path, current_region_file,
                                            sep = ""), header = FALSE)

#Set region boundries
# Comments by Rijan: Get the boundries of the regions from the columns. It probably makes sense to just keep these in headers and query them later?

region_chr <- region_of_interest$V1
region_start <- region_of_interest$V2
region_stop <- region_of_interest$V3
region_name <- gsub(x = current_region_file, pattern = ".region", replacement = "")

# This is only filtering the region file for rows and SNPs of interest.

region_gene_subset <- gene_location %>% subset(Chrom == as.character(region_chr) & Ext_Start >= region_start &
                                                 Ext_Start <= region_stop)

# filter for the genes of interest.

region_gene_hits_only <- region_gene_subset$Gene

## comments by Rijan: Just testing for the size of the gene list produced to not process something that is too small to process.

##annotate LD tracking DF


## comments by Rijan: ld_tracking_df was created above.
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

# This filters the genomic data for the chromsome of interest and the range of interest as defined by the LD
region_gene_subset <- gene_location %>% subset(Chrom == as.character(region_chr) & Ext_Start >= region_start &
                                                 Ext_Start <= region_stop)
}


# output the tables to file
write.table(x = region_gene_subset, file = paste(full_path, "/results/", "snp_region_full_gene_lists/", region_name,
                                                 "_all_genes.region", sep=""), row.names = F, sep = "\t")
write.table(x = region_gene_hits_only, file = paste(full_path, "/working/", "gene_lists_only_by_region/",  
                                                    region_name, "_just_genes.region", sep=""), row.names = F, col.names = F, sep = "\t")                                                     
}

write.table(x = ld_tracking_df, file = paste(full_path, "/results/LD_tracking_by_region_list.txt", sep=""), row.names = F, sep = "\t")