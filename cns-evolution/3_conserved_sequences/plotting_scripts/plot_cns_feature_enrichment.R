library(ggplot2)
library(tidyverse)
library(reshape2)
library(ape)
library(phytools)




bp_filt <- read.delim("../results/phastcons_filtered_enrichment_26082024.txt",sep='\t',header=TRUE)
cns_dist <- read.delim("../data/conservation/panand_raw_phastCons_gene_stats.bed", sep='\t',header=F)

# Feature enrichment within CNS

bp_filt <- bp_filt[bp_filt$Feature!="TE" & bp_filt$Feature!="TFBS",]
bp_filt$Feature <- gsub("ACR","open chromatin",bp_filt$Feature)
bp_filt <- bp_filt %>% mutate(point_size = ifelse(PhastCons_Query_intersect_bp/1000000 > 20, 15,
                                        ifelse(PhastCons_Query_intersect_bp/1000000 > 5,10,
                                        ))
                    )

ggplot(bp_filt,aes(reorder(Feature,-Enrichment),log(Enrichment),color="grey")) + 
  geom_point(size=bp_filt$point_size,color="grey")+
  geom_text(aes(label = round(Enrichment,1)),color="grey", position = position_nudge(0,0.4))+
  theme_bw()+
  ylim(0,3.5)+
  ylab("log(Fold enrichment)")+
  xlab("")+
  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45,vjust = 0.5))+
  theme(legend.position="none")
ggsave("../figures/Fig_conserved_sequence_feature_enrichment.pdf",width = 8)

# CNS distance to genes

cns_dist <- cns_dist %>%  mutate(updown = ifelse(V13 < 0,"upstream","downstream"))
ggplot(cns_dist[cns_dist$V13!=0,],aes(x=V13,fill=updown)) + geom_histogram(binwidth = 200, center = 100)+
  xlim(-10000,10000)+
  theme_bw()+
  xlab("Distance from nearest gene in bp")+
  ylab("CNS elements")+
  theme(text = element_text(size = 20))+
  scale_fill_manual(values = c("#1b9e77","#d95f02")) 
ggsave("../figures/Fig_conserved_sequence_distance_from_gene.pdf")



