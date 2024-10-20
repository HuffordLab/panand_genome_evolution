library(tidyverse)
library(ggplot2)
library(ape)
library(reshape2)

# Load input files
div <- read.delim("../../2_multiple_alignment/data/pairwise_4dsites_divergence_from_B73.txt",sep = ",", header = TRUE)
gene_rates <- read.delim("../../5_syntenic_genes/data/singlecopy_gene_subgenome_orthologs_rate.txt",sep='\t',header=TRUE)
species_gene_rates <- read.delim("../../5_syntenic_genes/data/singlecopy_gene_species_orthologs_rate.txt",sep='\t',header=TRUE)

## Pre-calculate pairwise turnover
#tfbs <- read.delim("../output_files/Zm-B73-REFERENCE-NAM-5.0_1000bp_JASPARmotifs_translationSS_aln_denovo_withcoords_core_exp_chip_syn_counts_18122023.tsv",sep='\t',header=TRUE)
#turnover <- tfbs %>% group_by(query_subgenome) %>%
#  summarise(ref_motifs = sum(tot_ref),
#            hq_motifs = sum(ref_highqual),
#            lq_motifs = sum(ref_lowqual))
#write.csv(turnover,"../output/tfbs_turnover_statistics.csv",quote = FALSE,row.names = FALSE)
turnover <- read.csv("../output/tfbs_turnover_statistics.csv",header = TRUE)

turnover$turnover <- 1- ((turnover$hq_motifs+turnover$lq_motifs)/turnover$ref_motifs)
turnover_distg <- merge(turnover_dist,gene_rates,by.x ="query_subgenome", by.y = "subgenome")
turnover_distg <- turnover_distg[c("species",
                                   "query_subgenome",
                                   "distance",
                                   "turnover.y",
                                   "turnover.x",
                                   "gene",
                                   "zmays_genes")]

colnames(turnover_distg) <- c("species",
                              "query_subgenome",
                              "distance",
                              "turnover_gene",
                              "turnover_tfbs",
                              "genes",
                              "zmays_genes")

turnover_distg <- melt(turnover_distg,id=c("query_subgenome","distance","species","genes","zmays_genes"))
turnover_distg <- turnover_distg[turnover_distg$variable!="turnover_gene",]

# use only subgenome 1
# remove znicar due to missing gene annotation impact
single_subgenome <- c("tdacts","zdiplg","zdiplm","tdactn","zhuehu","zTIL18","zTIL25","zTIL01",
                      "zTIL11","vcuspi@1","snutan@1","sscopa@1","rrottb@1","telega@1","aburma@1",
                      "ttrian@1","etrips@1","hcompr@1","blagur@1","ppanic@1","agerar@1","ccitra@1",
                      "udigit@1","hconto@1","rtuber@1","achine@1","smicro@1","crefra@1","avirgi@1",
                      "irugos@1","sbicol@1")
trips <- c("tdacts","zdiplg","znicar","zdiplm","tdactn","zhuehu","zTIL18","zTIL25","zTIL01",
           "zTIL11")
turnover_distg_single <- turnover_distg[ turnover_distg$query_subgenome %in% single_subgenome, ]
turnover_distg_single <- turnover_distg_single %>% mutate(clade = ifelse(species %in% trips,"Tripsacineae","Other_Andropogoneae"))

# Add divergence data
turnover_distg_single <- merge(turnover_distg_single,div,by = c("species"))

# Include additional line excluding Tripsacineae
turnover_distg_single$clade_variable <- paste(turnover_distg_single$clade,turnover_distg_single$variable,sep='_')
turnover_distg_single_tmp <- turnover_distg_single[turnover_distg_single$variable=="turnover_tfbs",]
turnover_distg_single_tmp$clade_variable <- "All"
turnover_distg_single_big <- rbind(turnover_distg_single,turnover_distg_single_tmp)
turnover_distg_single_big$clade_variable <- gsub("Other_Andropogoneae_turnover_species_gene","All_turnover_species_gene",turnover_distg_single_big$clade_variable)
turnover_distg_single_big$clade_variable <- gsub("Tripsacineae_turnover_species_gene","All_turnover_species_gene",turnover_distg_single_big$clade_variable)
turnover_distg_single_big$clade_variable <- gsub("Other_Andropogoneae_turnover_tfbs_random","All_turnover_tfbs_random",turnover_distg_single_big$clade_variable)
turnover_distg_single_big$clade_variable <- gsub("Tripsacineae_turnover_tfbs_random","All_turnover_tfbs_random",turnover_distg_single_big$clade_variable)
turnover_distg_single_big <- turnover_distg_single_big[turnover_distg_single_big$clade_variable!="Tripsacineae_turnover_tfbs",]

turnover_distg_single_big_gene <- turnover_distg_single_big
turnover_distg_single_big_gene$variable <- "turnover_gene"
turnover_distg_single_big_gene <- merge(turnover_distg_single_big_gene,subset(species_gene_rates,select=-c(genes,zmays_genes)),by = c("species"))
turnover_distg_single_big_gene$value <- turnover_distg_single_big_gene$turnover
turnover_distg_single_big_gene$clade <- "Andropogoneae"
turnover_distg_single_big_gene$clade_variable <- "Andropogoneae"
turnover_distg_single_big_gene <- subset(turnover_distg_single_big_gene, select=-c(turnover))
turnover_bind <- rbind(turnover_distg_single_big,turnover_distg_single_big_gene)


turnover_bind$clade <- ifelse(turnover_bind$species %in% trips, 'Tripsacineae', 'Other_Andropogoneae')
turnover_bind$clade_variable <- gsub("All", "Andropogoneae TFBS", turnover_bind$clade_variable )
turnover_bind$clade_variable <- gsub("^Other_Andropogoneae_turnover_tfbs$", "Andropogoneae ex Tripsacineae TFBS", turnover_bind$clade_variable )
turnover_bind$clade_variable <- gsub("^Andropogoneae$", "Andropogoneae genes", turnover_bind$clade_variable )

my_x_title <- expression(paste("% divergence from ", italic("Z"),".",italic(" mays")))
ggplot(turnover_bind,aes(y=value,x=divergence,color=clade_variable)) + 
  geom_point(aes(shape=factor(clade),size=3),colour = "black", fill="grey", size = 5)+
  geom_smooth(method='lm', formula= y~x,fullrange=TRUE)+
  theme_bw()+
  xlim(0,20)+
  xlab(my_x_title)+
  ylab("Turnover rate")+
  theme(text = element_text(size = 18))+
  scale_shape_manual(values=c(21,24))+
  scale_colour_manual(values=c("#d1641a","#7770b2","#369d77"))
ggsave("Figure_TFBS_turnover.pdf",width = 9, height = 5)
