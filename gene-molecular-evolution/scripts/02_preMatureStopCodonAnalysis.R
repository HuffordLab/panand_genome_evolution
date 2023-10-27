# R script to identify and test for non-random occurrence of premature stop codon
# author: Sheng-Kai Hsu
# date created: 2023.08.01
# data last edited: 2023.10.27
# update with results from orthofinder

rm(list=ls())
library(MSA2dist)
library(Biostrings)
library(parallel)
library(limma)

# customize code to generate contingency table with classifier, query and background sets
cont_table=function(query,background,classifyer){
  p1=length(intersect(query,classifyer))
  q1=length(query)-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}

# read alignment
aln.flist = list.files("/workdir/sh2246/p_panAndOGASR/output/orthoFinderOGMSA_noSevir",full.names = T,
                       pattern = "*.gs.l0.fa$")
aln.list = mclapply(aln.flist,readDNAStringSet,mc.cores = 20)
names(aln.list) = substr(strsplit2(aln.flist,"/")[,7],1,9)

# get a contingency table of presence of premature stop codon per MSA
##            perennial  annual 
##  presence      a         b
##  absence       c         d
cont_tab_list = mclapply(aln.list,function(a){
  codonArray = sapply(a,function(x) substring(as.character(x),seq(1,nchar(x),3),seq(3,nchar(x),3)))
  
  stopSite = na.omit(unlist(apply(codonArray,2,function(x) grep("TAA|TAG|TGA",x)[1])))
  if (any(stopSite<nrow(codonArray)*0.9)){idx = names(which(stopSite<nrow(codonArray)*0.9))
  } else idx = NULL
  
  annualPattern = "Zm|Zh|Zl|Zn|Zv|Zx|Sobic|Te"
  cont_tab = cbind(table(grepl(annualPattern,idx))[c("TRUE","FALSE")],
                   table(grepl(annualPattern,setdiff(colnames(codonArray),idx)))[c("TRUE","FALSE")])
  colnames(cont_tab) = c("TRUE","FALSE")
  rownames(cont_tab) = c("annual","perennial")
  cont_tab[is.na(cont_tab)] = 0
  return(cont_tab)                               
},mc.cores = 20)

# specify the outgroup reference sequences ID, and name the list of contingency tables
ref = mclapply(aln.list,function(x) {
  tmp = names(x)[grep("Pavag",names(x))]
  if (length(tmp)>1) tmp = sample(tmp,1)
  if (length(tmp)==0) tmp = NULL
  return(tmp)
})
names(cont_tab_list) = gsub("...p| ","",ref)

# Fisher's exact test on the contingency tables for non-random occurrence of premature stop codon
FET.p = mclapply(cont_tab_list,function(x) fisher.test(x)$p.value,mc.cores = 20)
FET.or = mclapply(cont_tab_list,function(x) fisher.test(x)$estimate,mc.cores = 20)                 

sum(p.adjust(FET.p,"BH")<0.05&FET.or>1)
sum(p.adjust(FET.p,"BH")<0.05&FET.or<1)
goi_more = names(which(p.adjust(FET.p,"BH")<0.05&FET.or>1)) # BH correction
goi_less = names(which(p.adjust(FET.p,"BH")<0.05&FET.or<1))
chisq.test(c(length(goi_more),length(goi_less)),p = c(.5,.5))

# Output test results
outTab = data.frame(OGID = names(aln.list),pvID = names(FET.p),p = unlist(FET.p), padj = p.adjust(FET.p,"BH"),OR = unlist(FET.or))
write.table(outTab,"/workdir/sh2246/p_panAndOGASR/output/prematureStopCodonEnrichment_ogv9.2b.txt",col.names = T,row.names = F,quote = F,sep = "\t")


# GO enrichment analysis
library(topGO)
GODB <- readMappings("/workdir/sh2246/p_panAndOGASR/data/Pv_GOTable.txt",IDsep = ",")
background = intersect(names(GODB),paste0(names(FET.p),".1"))
# background = intersect(names(GODB),paste0(names(diffP),".1"))

tmp=factor(as.integer(background%in%paste0(goi_more,".1")))
names(tmp)=background
tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.gene2GO, gene2GO = GODB)
resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
GO_res_table_more=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)

tmp=factor(as.integer(background%in%paste0(goi_less,".1")))
names(tmp)=background
tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.gene2GO, gene2GO = GODB)
resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
GO_res_table_less=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)


# comparison to phyloP-based results
phyloP = read.table("/workdir/sh2246/p_panAndOGASR/output/phyloP_og.filtered.txt",header = T)

bonferroni_cutoff = -log10(0.01/nrow(phyloP))
bonferroni_cutoff
goi_acc = paste0(phyloP[phyloP[,10] < -bonferroni_cutoff ,1])
goi_con = paste0(phyloP[phyloP[,10] > bonferroni_cutoff ,1])

cont_table(goi_acc,phyloP[,1],goi_more)
fisher.test(cont_table(goi_acc,phyloP[,1],goi_more),alternative = "greater")

