#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##Load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(gtable)
library(patchwork)
library(stringr)
library(gridExtra)
library(grid)

##Make folders

dir.create("../../outputs")
dir.create("../../outputs/summary_tables")
dir.create("../../outputs/diversity_pdfs")

##Load resource files
gene_db <- read.table("../../../resources/R_resources/Sviridis_311_v1.1.gene.gff_for_gene_models.txt", sep = "\t", header = T)

##models

##process transcript ID to text
transcript_to_geneid <- function(transcript){
  print(transcript)
}

##gene model functions
single_gene_to_model <- function(gene_id) {
  ##define target gene
  target_gene <- gene_id
  base_gene <- gsub(target_gene,pattern = "\\.[1-9]$",replacement = "")
  ##define the gene bounds
  gene_bounds <- subset(gene_db, GeneID == base_gene)
  gene_bounds$gene_number <- as.factor(1:nrow(gene_bounds))
  expand_factor <- (gene_bounds$Stop - gene_bounds$Start) * .025
  gene_bounds$Start1 <- round(gene_bounds$Start - expand_factor)
  gene_bounds$Stop1 <- round(gene_bounds$Stop + expand_factor)
  gene_bounds$Start <- NULL
  gene_bounds$Stop <- NULL
  gene_bounds$Strand <- lapply(gene_bounds$Strand, function(x) { gsub("\\-", "first", x) })
  gene_bounds$Strand <- lapply(gene_bounds$Strand, function(x) { gsub("\\+", "last", x) })
  
  ##Melt gene info for ggplot
  gene_space <- melt(gene_bounds, id.vars = c("GeneID", "Chrom", "Type", "Strand", "gene_number"))
  
  ##Define exon boundries
  exons_bounds <- gene_db %>% filter(Type == "CDS" & GeneID == target_gene)
  exons_bounds$observation <- as.factor(1:nrow(exons_bounds))
  exons_bounds$exon_label <- paste(exons_bounds$Type,exons_bounds$observation)
  exons_bounds$y1 <- .98
  exons_bounds$y2 <- 1.02
  
  ##Define 5'UTR
  five_prime_bounds <- gene_db %>% filter(Type == "five_prime_UTR" & GeneID == target_gene)
  
  if(nrow(five_prime_bounds)>=1){
    five_prime_bounds$y1 <- .98 
    five_prime_bounds$y2 <- 1.02 
    five_prime_status <- TRUE
  } else{
    five_prime_status <- FALSE
  }
  
  ##Define 3'UTR
  three_prime_bounds <- gene_db %>% filter(Type == "three_prime_UTR" & GeneID == target_gene)
  
  if(nrow(three_prime_bounds)>=1){
    three_prime_bounds$y1 <- .98 
    three_prime_bounds$y2 <- 1.02 
    three_prime_status <- TRUE
  } else{
    three_prime_status <- FALSE
  }
  
  ##Define X-axis
  
  gene_xlabel <- paste("Position on ", gene_bounds$Chrom,sep = "")
  
  ##convery gene_number to numeric
  gene_space$gene_number=as.numeric(levels(gene_space$gene_number))[gene_space$gene_number]
  ##Draw graph
  gene_model <- ggplot(data=gene_space, aes(x=value, y=gene_number, group=GeneID)) +
    geom_line(arrow = arrow(length=unit(0.30,"cm"), ends=paste(gene_bounds$Strand), type = "closed"), size = 1)+
    geom_point(alpha = 0.05) + scale_y_continuous(limits = c(.93, 1.07))
  
  gene_model <- gene_model + geom_rect(data=exons_bounds, inherit.aes=FALSE,
                                       aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2,
                                           group=exon_label), fill = "dodgerblue1") 
  
  ##Add 5' UTR if there is one
  if(five_prime_status==TRUE){
    gene_model <- gene_model + geom_rect(data=five_prime_bounds, inherit.aes=FALSE,
                                         aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2), fill = "chartreuse4") + xlab(gene_xlabel) +
      ggtitle(paste(target_gene)) +  theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
            axis.line.y=element_blank())
    
  }
  
  ##Add 3' UTR if there is one
  if(three_prime_status==TRUE){
    gene_model <- gene_model + geom_rect(data=three_prime_bounds, inherit.aes=FALSE,
                                         aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2), fill = "orange2") + xlab(gene_xlabel) +
      ggtitle(paste(target_gene)) +  theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
            axis.line.y=element_blank())+ 
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    
  }  
  
  if(three_prime_status == FALSE && five_prime_status == FALSE){
    gene_model <- gene_model + xlab(gene_xlabel) +
      ggtitle(paste(target_gene)) +  theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x = element_line(colour = "black"), 
            axis.line.y=element_blank())+ 
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
  }
  
  gene_model <- gene_model + theme(axis.text=element_text(size=14),
                                   axis.title=element_text(size=20,face="bold")) + 
  theme(plot.title = element_text(size=25, face ="bold"))
  
  return(gene_model)
}
single_gene_to_model_no_axis <- function(gene_id) {
  ##define target gene
  target_gene <- gene_id
  base_gene <- gsub(target_gene,pattern = "\\.[1-9]",replacement = "")
  
  ##define the gene bounds
  gene_bounds <- subset(gene_db, GeneID == base_gene)
  gene_bounds$gene_number <- as.factor(1:nrow(gene_bounds))
  expand_factor <- (gene_bounds$Stop - gene_bounds$Start) * .025
  gene_bounds$Start1 <- round(gene_bounds$Start - expand_factor)
  gene_bounds$Stop1 <- round(gene_bounds$Stop + expand_factor)
  gene_bounds$Start <- NULL
  gene_bounds$Stop <- NULL
  gene_bounds$Strand <- lapply(gene_bounds$Strand, function(x) { gsub("\\-", "first", x) })
  gene_bounds$Strand <- lapply(gene_bounds$Strand, function(x) { gsub("\\+", "last", x) })
  
  ##Melt gene info for ggplot
  gene_space <- melt(gene_bounds, id.vars = c("GeneID", "Chrom", "Type", "Strand", "gene_number"))
  
  ##Define exon boundries
  exons_bounds <- gene_db %>% filter(Type == "CDS" & GeneID == target_gene)
  exons_bounds$observation <- as.factor(1:nrow(exons_bounds))
  exons_bounds$exon_label <- paste(exons_bounds$Type,exons_bounds$observation)
  exons_bounds$y1 <- .98
  exons_bounds$y2 <- 1.02

  ##Define 5'UTR
  five_prime_bounds <- gene_db %>% filter(Type == "five_prime_UTR" & GeneID == target_gene)
  
  if(nrow(five_prime_bounds)>=1){
    five_prime_bounds$y1 <- .98 
    five_prime_bounds$y2 <- 1.02 
    five_prime_status <- TRUE
  } else{
    five_prime_status <- FALSE
  }
  
  ##Define 3'UTR
  three_prime_bounds <- gene_db %>% filter(Type == "three_prime_UTR" & GeneID == target_gene)
  
  if(nrow(three_prime_bounds)>=1){
    three_prime_bounds$y1 <- .98 
    three_prime_bounds$y2 <- 1.02 
    three_prime_status <- TRUE
  } else{
    three_prime_status <- FALSE
  }
  
  ##Define X-axis
  
  gene_xlabel <- paste("Position on ", gene_bounds$Chrom,sep = "")
  
  ##convery gene_number to numeric
  gene_space$gene_number=as.numeric(levels(gene_space$gene_number))[gene_space$gene_number]
  
  ##Draw graph
  gene_model <- ggplot(data=gene_space, aes(x=value, y=gene_number, group=GeneID)) +
    geom_line(arrow = arrow(length=unit(0.30,"cm"), ends=paste(gene_bounds$Strand), type = "closed"), size = 1)+
    geom_point(alpha = 0.05) + scale_y_continuous(limits = c(.93, 1.07))
  
  gene_model <- gene_model + geom_rect(data=exons_bounds, inherit.aes=FALSE,
                                       aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2,
                                           group=exon_label), fill = "steelblue") 
  
  ##Add 5' UTR if there is one
  if(five_prime_status==TRUE){
    gene_model <- gene_model + geom_rect(data=five_prime_bounds, inherit.aes=FALSE,
                                         aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2), fill = "green") + xlab(gene_xlabel) +
      ggtitle(paste(target_gene)) +  theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x=element_blank(), 
            axis.line.y=element_blank())
  }
  
  ##Add 3' UTR if there is one
  if(three_prime_status==TRUE){
    gene_model <- gene_model + geom_rect(data=three_prime_bounds, inherit.aes=FALSE,
                                         aes(xmin=Start,xmax=Stop,ymin=y1,ymax=y2), fill = "red") + xlab(gene_xlabel) +
      ggtitle(paste(target_gene)) +  theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x=element_blank(), 
            axis.line.y=element_blank())
  }  
  
  if(three_prime_status == FALSE && five_prime_status == FALSE){
    gene_model <- gene_model + xlab(gene_xlabel) +
      ggtitle(paste(target_gene)) +  theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.x=element_blank(), 
            axis.line.y=element_blank())
  }
  return(gene_model)
}
single_gene_model_mutations <- function(gene_id){
  
  ##Create base model
  gene_model_diagram <- single_gene_to_model(gene_id)
  
  ##set graph limits for alternating mutations
  y_limits <- as.matrix(rep(c(.98, .96),10000))
  yend_limits <- as.matrix(rep(c(1.04, 1.02),10000))
  point_end_limits <- as.matrix(rep(c(1.04, .96),10000))
  
  ##Extract mutation data for gene of interest
  line_mutations <- mutation_db_filt %>% filter(ID == gene_id)
  
  ##Add locations for plotting mutations if available
  if(nrow(line_mutations)>=1){
    line_mutations <- cbind(line_mutations, y=y_limits[1:nrow(line_mutations),]) 
    line_mutations <- cbind(line_mutations, yend=yend_limits[1:nrow(line_mutations),])
    line_mutations <- cbind(line_mutations, point_end=point_end_limits[1:nrow(line_mutations),])
    line_mutations_status <- TRUE
  } else{
    line_mutations_status <- FALSE
  }
  
  #Create color scheme for mutation effect
  line_mutations$colors <- line_mutations$EFFECT
  line_mutations$colors <- gsub(line_mutations$colors, pattern = "HIGH", replacement = "orangered")
  line_mutations$colors <- gsub(line_mutations$colors, pattern = "MODERATE", replacement = "springgreen2")
  line_mutations$colors <- gsub(line_mutations$colors, pattern = "LOW", replacement = "gray60")
  
  ##Modify graph if there are mutations
  if(line_mutations_status == TRUE){
    gene_model_diagram <- gene_model_diagram + 
      geom_segment(data = line_mutations,  aes(x = POS, y = y, xend = POS, yend = yend), 
                   inherit.aes=FALSE, color = line_mutations$colors, size = .8) +
      geom_point(data = line_mutations, aes(x = POS, y = point_end), 
                 inherit.aes=FALSE, fill = line_mutations$colors, size=line_mutations$size, 
                 color = "black", alpha = .7, shape = line_mutations$shape)
  }
  
  gene_model_diagram <- gene_model_diagram + theme(axis.text=element_text(size=14),
                                                   axis.title=element_text(size=20,face="bold")) + theme(plot.title = element_text(size=25, face ="bold"))
  return(gene_model_diagram)
  
}


## tables

make_gene_mutation_table <- function(gene_id){

  mutation_db_filt_tables <- mutation_db_filt %>% filter(ID == transcript_to_geneid(gene_id))
  mutation_table_db <- data.frame(Chromosome=mutation_db_filt_tables$CHROM, Position=mutation_db_filt_tables$POS,
                                  Reference=as.character(mutation_db_filt_tables$REF), Alternate=as.character(mutation_db_filt_tables$ALT),
                                  "Mutation_type"=as.character(mutation_db_filt_tables$TYPE), 
                                  "Predicted_effect" = as.factor(mutation_db_filt_tables$EFFECT),
                                  "Amino_acid_change" = as.character(mutation_db_filt_tables$AA_SUB),
                                  "Percent_Alternate" = paste(round(mutation_db_filt_tables$PER_ALT), "%", sep = ""))
  
  
  table <- tableGrob(mutation_table_db, theme = ttheme_default(base_size = 8))
  title <- textGrob(paste(gene_id, "genic mutations within sorghum BAP"),gp=gpar(fontsize=16, fontface = "bold"))
  footnote <- textGrob("", x=0, hjust=0,
                       gp=gpar( fontface="italic"))
  padding <- unit(0.5,"line")
  table <- gtable_add_rows(table, 
                           heights = grobHeight(title) + padding,
                           pos = 0)
  table <- gtable_add_rows(table, 
                           heights = grobHeight(footnote)+ padding)
  table <- gtable_add_grob(table, list(title, footnote),
                           t=c(1, nrow(table)), l=c(1,2), 
                           r=ncol(table))
  grid.newpage()
  return(grid.draw(table))
  
}

##Calculate N samples for geno/pheno plots
give.n <- function(x){
  return(c(y = median(x) * 1.02, label = length(x)))
}

#create geno/pheno plots
create_geno_pheno_boxplot <- function(gene){
  
  effects_of_interest <- c("HIGH", "MODERATE")
  pheno_geno_table <- pre_mutation_db %>% filter(EFFECT %in% effects_of_interest & ID == gene)
  pheno_geno_table$PER_ALT <- NULL
  
  if(dim(pheno_geno_table)[1] >= 1){
    genos_table <- as.data.frame(t(pheno_geno_table[,9:length(pheno_geno_table)]))
    genos_table2 <- data.frame(lapply(genos_table, function(x) {gsub(" ", "", x)}))
    genos_table2 <- data.frame(lapply(genos_table2, function(x) {gsub("0", "Ref/Ref", x)}))
    genos_table2 <- data.frame(lapply(genos_table2, function(x) {gsub("1", "Ref/Alt", x)}))
    genos_table2 <- data.frame(lapply(genos_table2, function(x) {gsub("2", "Alt/Alt", x)}))
    pheno_geno_table$names <- paste(pheno_geno_table$CHROM, pheno_geno_table$POS, pheno_geno_table$TYPE, sep = "_")
    colnames(genos_table2) <- t(pheno_geno_table$names)
    traits_to_graph <- colnames(genos_table2)
    
    genos_table_final <- cbind(PlantID=rownames(genos_table),genos_table2) 
    geno_pheno_full_table <- merge(phenotype, genos_table_final, by.x = "Genotype", by.y="PlantID", all=F)
    geno_pheno_full_table$Genotype <- NULL
    
    
    
    for(current_trait in traits_to_graph){
      
      ##pull data for current phenotype
      current_pheno_geno_table <- as.data.frame(cbind(geno_pheno_full_table[1], geno_pheno_full_table[current_trait]))
      current_pheno_geno_table <- current_pheno_geno_table[complete.cases(current_pheno_geno_table),]
      current_pheno_geno_table[[current_trait]] <- factor(current_pheno_geno_table[[current_trait]], levels = c("Ref/Ref",
                                                                                                                "Ref/Alt",
                                                                                                                "Alt/Alt"))
      current_pheno_geno_table <- current_pheno_geno_table[complete.cases(current_pheno_geno_table),]
      current_pheno_geno_table <- droplevels(current_pheno_geno_table)
      
      ##Run Anova if all three genotype states exist
      if(length(levels(current_pheno_geno_table[,2])) == 3){
        res.aov <- aov(current_pheno_geno_table[,1] ~ current_pheno_geno_table[,2])
        
        
        #make ANovA table
        Tukey_table <- TukeyHSD(res.aov, ordered = TRUE)
        pg_dt_tukey <- as.data.frame(Tukey_table$`current_pheno_geno_table[, 2]`)
        
        ##Make title for Tukey table
        title <- textGrob("TukeyHSD Multiple Pairwise-Comparions Summary",gp=gpar(fontsize=10), hjust=.4)
        padding <- unit(5,"mm")
        
        ##Make Tukey Table
        
        t1 <- tableGrob(pg_dt_tukey, theme = ttheme_default(base_size = 8)  )
        
        table <- gtable_add_rows(
          t1, 
          heights = grobHeight(title) + padding,
          pos = 0)
        table <- gtable_add_grob(
          table, 
          title, 
          1, 1, 1, ncol(table))
        
        ##residuals plot
        
        Model <- data.frame(Fitted = fitted(res.aov), 
                            Residuals = resid(res.aov), 
                            Treatment = current_pheno_geno_table[,2])
        
        res.model_plot <- ggplot(Model, aes(Fitted, Residuals, colour = Treatment)) + 
          geom_point() + ggtitle("Residuals vs Fitted") + theme_classic() + 
          theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) 
        
        qq_plot <- ggplot(current_pheno_geno_table, aes(sample = current_pheno_geno_table[,1])) + 
          stat_qq() + stat_qq_line() + theme_classic() +  ggtitle("QQ-Plot") +  
          theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) 
        
        ##Make pheno/geno boxplot
        x_name <- colnames(current_pheno_geno_table)[2]
        y_name <- colnames(current_pheno_geno_table)[1]
        pheno_stats <- table(current_pheno_geno_table[2])
        p1 <- ggplot(data = current_pheno_geno_table, 
                     aes(x=current_pheno_geno_table[[2]], y=current_pheno_geno_table[[1]], 
                         col = current_pheno_geno_table[[2]])) + geom_boxplot() + theme_classic() + 
          theme(legend.position = "none") + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
          xlab(paste(x_name)) + ylab(paste(y_name)) + 
          ggtitle(paste(gene)) + stat_summary(fun.data = give.n, geom = "text", size = 6, colour = "black", position = "identity")
        
        
        print(grid.arrange(p1,table, res.model_plot, qq_plot, ncol = 2, clip = FALSE, 
                           widths = unit(c( .5, .5), "npc")))
      }
      
      ##Run t-test if only two allele states
      if(length(levels(current_pheno_geno_table[,2])) == 2){
        res.aov <- aov(current_pheno_geno_table[,1] ~ current_pheno_geno_table[,2])
        
        
        #make ANovA table
        Tukey_table <- TukeyHSD(res.aov, ordered = TRUE)
        pg_dt_tukey <- as.data.frame(Tukey_table$`current_pheno_geno_table[, 2]`)
        
        ##Make title for Tukey table
        title <- textGrob("TukeyHSD Multiple Pairwise-Comparions Summary",gp=gpar(fontsize=10), hjust=.4)
        padding <- unit(5,"mm")
        
        ##Make Tukey Table
        
        t1 <- tableGrob(pg_dt_tukey, theme = ttheme_default(base_size = 8))
        
        table <- gtable_add_rows(
          t1, 
          heights = grobHeight(title) + padding,
          pos = 0)
        table <- gtable_add_grob(
          table, 
          title, 
          1, 1, 1, ncol(table))
        
        ##residuals plot
        
        Model <- data.frame(Fitted = fitted(res.aov), 
                            Residuals = resid(res.aov), 
                            Treatment = current_pheno_geno_table[,2])
        
        res.model_plot <- ggplot(Model, aes(Fitted, Residuals, colour = Treatment)) + 
          geom_point() + ggtitle("Residuals vs Fitted") + theme_classic() + 
          theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) 
        
        qq_plot <- ggplot(current_pheno_geno_table, aes(sample = current_pheno_geno_table[,1])) + 
          stat_qq() + stat_qq_line() + theme_classic() +  ggtitle("QQ-Plot") +  
          theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) 
        
        ##Make pheno/geno boxplot
        x_name <- colnames(current_pheno_geno_table)[2]
        y_name <- colnames(current_pheno_geno_table)[1]
        pheno_stats <- table(current_pheno_geno_table[2])
        p1 <- ggplot(data = current_pheno_geno_table, 
                     aes(x=current_pheno_geno_table[[2]], y=current_pheno_geno_table[[1]], 
                         col = current_pheno_geno_table[[2]])) + geom_boxplot() + theme_classic() + 
          theme(legend.position = "none") + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
          xlab(paste(x_name)) + ylab(paste(y_name)) + 
          ggtitle(paste(gene)) + stat_summary(fun.data = give.n, geom = "text", size = 6, colour = "black", position = "identity")
        
        
        print(grid.arrange(p1,table, res.model_plot, qq_plot, ncol = 2, clip = FALSE, 
                           widths = unit(c( .5, .5), "npc")))
      }
      
      ##If only one allele state exists in the geno/phenos
      if(length(levels(current_pheno_geno_table[,2])) == 1){
        ##Make pheno/geno boxplot
        x_name <- colnames(current_pheno_geno_table)[2]
        y_name <- colnames(current_pheno_geno_table)[1]
        pheno_stats <- table(current_pheno_geno_table[2])
        
        p1 <- ggplot(data = current_pheno_geno_table, 
                     aes(x=current_pheno_geno_table[[2]], y=current_pheno_geno_table[[1]], 
                         col = current_pheno_geno_table[[2]])) + geom_boxplot() + theme_classic() + 
          theme(legend.position = "none") + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
          xlab(paste(x_name)) + ylab(paste(y_name)) + 
          ggtitle(paste(gene)) + stat_summary(fun.data = give.n, geom = "text", size = 6, colour = "black", position = "identity")
        
        
        print(p1)
        
    }
    
}
  }
}

##get list of snp_targets

snp_targets <- list.files(path=".")

paste("Processing: ", snp_targets)

##read and process if phenotype data available
# test if there is at least one argument: if not, return an error
if (length(args)==1) {
  
  ##read in phenotype file
  phenotype_file <- args[1]
  path_to_pheno_file <- paste("../../../", phenotype_file, sep="")
  phenotype <- read.table(path_to_pheno_file, header = T, sep = "\t")
  
  #process all snps with phenotype info
  for(single_snp_target in snp_targets){
    
    path_to_target <- paste(single_snp_target,"/",single_snp_target,"_impactful_gene_hits.txt", sep="")
    
    ##make folders for each snp_target
    diversity_pdfs_specific <- paste("../../outputs/diversity_pdfs/", single_snp_target, sep = "")
    dir.create(diversity_pdfs_specific)
    
    
    ##Import mutations of interest and calculate mutation frequency
    pre_mutation_db <- read.table(path_to_target, sep = "\t", header = T)
    
    pre_mutation_db$ID <- gsub(x = pre_mutation_db$ID, pattern = ".v1.1", replacement = "")
    number_genos <- length(pre_mutation_db) - 9
    pre_mutation_db$PER_ALT = ((rowSums(pre_mutation_db[,9:length(pre_mutation_db)], na.rm = T) / 2 ) / number_genos ) * 100
    
    
    ##Create and add additional columns for graphing
    mutation_db <- cbind(pre_mutation_db[,1:8],PER_ALT=pre_mutation_db$PER_ALT)
    mutation_db$size <- mutation_db$PER_ALT/10
    mutation_db$shape <- ifelse(mutation_db$TYPE=='missense_variant',21, 
                                ifelse(mutation_db$TYPE=='synonymous_variant',22, 
                                       ifelse(mutation_db$TYPE=='disruptive_inframe_deletion',23,25)))
    
    mutation_db_filt <- droplevels(mutation_db %>% filter(size > 0))
    mutation_db_filt$ID <- as.factor(mutation_db_filt$ID)
    
    ##Make summary output file
    ##make empty dataframe 
    
    odat <- data.frame()
    
    output_db <- pre_mutation_db
    output_db$ID <- as.factor(output_db$ID)
    genes_of_interest <- levels(output_db$ID)
    
    ##loop through and generate summary tables
    for( ID_sample in genes_of_interest){
      
      sub_db <- output_db %>% filter(ID == ID_sample)
      gene <- t(as.data.frame(table(sub_db$EFFECT)))
      colnames(gene) <- gene[1,]
      gene <- gene[-c(1),]
      gene <- t(as.data.frame(gene))
      
      types <- t(as.data.frame(table(sub_db$TYPE)))
      colnames(types) <- types[1,]
      types <- types[-c(1),]
      types <- t(as.data.frame(types))
      
      temp_df <- cbind(ID_sample, gene, types)
      odat <- rbind(odat, temp_df)
    }
    
    output_summary_name_with_path <- paste("../../outputs/summary_tables/", single_snp_target, "_summary_output.txt", sep = "")
    write.table(odat, output_summary_name_with_path, sep = "\t", row.names = F)
    
    
    
    ##Create one file for gene of interest
    for(gene in levels(mutation_db_filt$ID)){
      pdf_name <- paste(gene,"diversity_summary.pdf", sep = "_")
      pdf(paste(diversity_pdfs_specific, "/", pdf_name, sep = ""), onefile = TRUE, width = 12)
      print(single_gene_model_mutations(gene))
      print(make_gene_mutation_table(gene))
      print(create_geno_pheno_boxplot(gene))
      dev.off()
    }
  }
  
  
  
  } else if (length(args)==0) {

    #process all snps without phenotype info
    for(single_snp_target in snp_targets){
      
      path_to_target <- paste(single_snp_target,"/",single_snp_target,"_impactful_gene_hits.txt", sep="")
      
      ##make folders for each snp_target
      diversity_pdfs_specific <- paste("../../outputs/diversity_pdfs/", single_snp_target, sep = "")
      dir.create(diversity_pdfs_specific)
      
      
      ##Import mutations of interest and calculate mutation frequency
      pre_mutation_db <- read.table(path_to_target, sep = "\t", header = T)
      
      pre_mutation_db$ID <- gsub(x = pre_mutation_db$ID, pattern = ".v1.1", replacement = "")
      number_genos <- length(pre_mutation_db) - 9
      pre_mutation_db$PER_ALT = ((rowSums(pre_mutation_db[,9:length(pre_mutation_db)], na.rm = T) / 2 ) / number_genos ) * 100
      
      
      ##Create and add additional columns for graphing
      mutation_db <- cbind(pre_mutation_db[,1:8],PER_ALT=pre_mutation_db$PER_ALT)
      mutation_db$size <- mutation_db$PER_ALT/10
      mutation_db$shape <- ifelse(mutation_db$TYPE=='missense_variant',21, 
                                  ifelse(mutation_db$TYPE=='synonymous_variant',22, 
                                         ifelse(mutation_db$TYPE=='disruptive_inframe_deletion',23,25)))
      
      mutation_db_filt <- droplevels(mutation_db %>% filter(size > 0))
      mutation_db_filt$ID <- as.factor(mutation_db_filt$ID)
      
      ##Make summary output file
      ##make empty dataframe 
      
      odat <- data.frame()
      
      output_db <- pre_mutation_db
      output_db$ID <- as.factor(output_db$ID)
      genes_of_interest <- levels(output_db$ID)
      
      ##loop through and generate summary tables
      for( ID_sample in genes_of_interest){
        
        sub_db <- output_db %>% filter(ID == ID_sample)
        gene <- t(as.data.frame(table(sub_db$EFFECT)))
        colnames(gene) <- gene[1,]
        gene <- gene[-c(1),]
        gene <- t(as.data.frame(gene))
        
        types <- t(as.data.frame(table(sub_db$TYPE)))
        colnames(types) <- types[1,]
        types <- types[-c(1),]
        types <- t(as.data.frame(types))
        
        temp_df <- cbind(ID_sample, gene, types)
        odat <- rbind(odat, temp_df)
      }
      
      output_summary_name_with_path <- paste("../../outputs/summary_tables/", single_snp_target, "_summary_output.txt", sep = "")
      write.table(odat, output_summary_name_with_path, sep = "\t", row.names = F)
      
      
      
      ##Create one file for gene of interest
      for(gene in levels(mutation_db_filt$ID)){
        pdf_name <- paste(gene,"diversity_summary.pdf", sep = "_")
        pdf(paste(diversity_pdfs_specific, "/", pdf_name, sep = ""), onefile = TRUE, width = 12)
        print(single_gene_model_mutations(gene))
        print(make_gene_mutation_table(gene))
        dev.off()
      }
    }
}

#merge output tables for master table
all_summary_tables_list <- list.files(path = "../../outputs/summary_tables/.", pattern = "summary_output.txt")
master_table <- data.frame()
for( summary_table_single in all_summary_tables_list){
  SNP_Hit_ID <- gsub(x = summary_table_single, pattern = "_summary_output.txt", replacement = "")
  current_table <- read.table(paste("../../outputs/summary_tables/", summary_table_single, sep=""), header = T, sep = "\t")
  Gene_ID <- str_replace(pattern = '.[1-9]$', replacement = "", string = current_table$ID_sample)
  current_table <- cbind(SNP_Hit_ID, Gene_ID, current_table)
  master_table <- bind_rows(master_table, current_table)
  }


setaria_anno <- read.csv("../../../resources/R_resources/Sviridis_311_v1.1.annotation_info.csv", header = T)
final_output <- merge(master_table, setaria_anno, by.x = "Gene_ID", by.y = "locusName", all.x = T, all.y = F)
write.csv(final_output, "../All_queries_gene_summary_table.csv", row.names = F)



