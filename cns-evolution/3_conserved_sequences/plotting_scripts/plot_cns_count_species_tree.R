library(ggplot2)
library(tidyverse)
library(reshape2)
library(ape)
library(phytools)

df2matrix<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

setwd("C:\\Users\\Armin\\ws\\panand\\analysis\\genome_paper_analysis_30112023\\3_conserved_sequences")
pan_tree <- read.tree("C:\\Users\\Armin\\ws\\panand\\analysis\\panand_cns_analysis\\3_conserved_sequences\\data\\tree/panand_neutral_tree.nwk")
cns_pav <- read.delim("C:\\Users\\Armin\\ws\\panand\\analysis\\panand_cns_analysis\\3_conserved_sequences\\results\\all_species_cns_Nfilt_feat_bp_counts.txt",sep='\t',header = TRUE)
cns_pav_dcast <- dcast(cns_pav,species ~ feature, value.var = "bp")
#reorder columns
cns_pav_dcast <- cns_pav_dcast[,c(1,8,7,4,5,6,2,3)]
cns_pav_mat <- df2matrix(cns_pav_dcast)

cols<- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00")
plotTree.barplot(pan_tree,
                 cns_pav_mat/1000000,
                 add=TRUE,
                 args.plotTree=list(fsize=1,lwd=1,xlim=c(0,0.3)),
                 args.barplot=list(xlab="CNS elements (Mbp)",
                                   xlim=c(0.,40),
                                   col=cols))

