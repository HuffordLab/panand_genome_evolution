# R script to summarize phyloP results and GO enrichment
# author: Sheng-Kai Hsu
# date created: 2023.08.08
# data last edited: 2023.10.27
# update with results from orthofinder

#load package
rm(list=ls())
library(limma)
library(parallel)
library(topGO)

#read phyloP
phyloP_og = read.table("/workdir/sh2246/p_panAndOGASR/output/phyloP_res_og.txt",header = F)

phyloP_og[,1] = substr(phyloP_og[,1],1,14) #renaming 
phyloP_og$Score = -log10(1-pchisq(2*phyloP_og$V8,1))*ifelse(phyloP_og$V9>0,1,-1)

phyloP_og$Score[phyloP_og$Score==Inf]=20
phyloP_og$Score[phyloP_og$Score==-Inf]=-20

colnames(phyloP_og) = c("chr","start","end","name","null_scale","alt_scale","alt_subscale","lnlratio","pval","Score")
# write.table(phyloP_og,"/workdir/sh2246/p_panAndOGASR/output/phyloP_og.summaryTab.txt",quote = F,row.names = F,sep = "\t")

# candidate genes based on Bonferroni's threshold
bonferroni_cutoff = -log10(0.01/nrow(phyloP_og))
goi_acc_og = paste0(phyloP_og[phyloP_og[,10] < -bonferroni_cutoff ,1],".1")
goi_con_og = paste0(phyloP_og[phyloP_og[,10] > bonferroni_cutoff ,1],".1")

# chisq test for unproportional distribution of the accelerated and conserved genes
length(goi_acc_og)
length(goi_con_og)
chisq.test(c(length(goi_acc_og),length(goi_con_og)),p = c(.5,.5))

#GO enrichment analysis for the candidate genes
GODB <- readMappings("/workdir/sh2246/p_panAndOGASR/data/Pv_GOTable.txt",IDsep = ",")
# background = names(GODB)
background = intersect(names(GODB),paste0(phyloP_og[,1],".1"))
length(background)

tmp=factor(as.integer(background%in%goi_acc_og))
names(tmp)=background
tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.gene2GO, gene2GO = GODB)
resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
GO_res_table_acc_og=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)

tmp=factor(as.integer(background%in%goi_con_og))
names(tmp)=background
tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.gene2GO, gene2GO = GODB)
resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
GO_res_table_con_og=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)

write.table(GO_res_table_acc_og,"/workdir/sh2246/p_panAndOGASR/output/GO_phastAcc_og.txt",quote = F,sep = "\t",row.names = F)
write.table(GO_res_table_con_og,"/workdir/sh2246/p_panAndOGASR/output/GO_phastCon_og.txt",quote = F,sep = "\t",row.names = F)