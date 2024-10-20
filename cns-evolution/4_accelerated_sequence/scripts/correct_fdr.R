#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

options(scipen=999)
pvals <- read.delim(args[1],sep="\t",header=FALSE)
adj <- p.adjust(pvals$V1, method="fdr")
write.table(adj,"adj.tsv",sep = ",",row.names = FALSE, quote = FALSE)

