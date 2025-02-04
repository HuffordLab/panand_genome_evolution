library(rGREAT)
library(BioMartGOGeneSets)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(ggplot2)
library(tibble)


split_table_columns <- function(data_table) {
  result_list <- list()
  
  for (i in 1:nrow(data_table)) {
    name <- as.character(data_table[i, 1])
    values <- strsplit(as.character(data_table[i, 2]), ",\\s*")[[1]]
    result_list[[name]] <- values
  }
  
  return(result_list)
}

###################
# Custom analysis #
###################

# 1) We need a gene set that corresponds to GO terms

custom_go <- read.delim(gzfile("../data/Zm00001eb.1.fulldata_goList.go2gene.txt.gz"),header = FALSE,sep = '\t')
go_list <- split_table_columns(custom_go)

go_terms <- read.delim(gzfile("../data/go_term_lookup.tsv.gz"),header=TRUE,sep='\t')


# 2) we need a our target regions as a granges object, ie CNS 

gr_obj =  import(gzfile("../results/panand_cns_phylop_ACC_LRT_adjusted_signif_merged_50bp.bed.gz"))
background = import(gzfile("../data/panand_cns_50bpmerge.bed.gz"))



# 3) we need the gene annotation as a granges object
gdf = read.delim(gzfile("../data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.tss.tsv.gz"),sep='\t',header = FALSE)
tss = GRanges(seqnames = gdf[, 1], ranges = IRanges(gdf[, 2], gdf[, 2]), 
              strand = gdf[, 3], gene_id = gdf[, 4])
chrom = c("1","2","3","4","5","6","7","8","9","10")
chromsize = c(308452471,243675191,238017767,250330460,226353449,181357234,185808916,182411202,163004744,152435371)
chromlen = setNames(chromsize,chrom)  
et = extendTSS(tss, chromlen)

# Merged CNS within 50bp
#background = import("data\\panand_cns_50bpmerge.bed")
resttest = great(background, go_list,extended_tss = et)

# recalculate p-values with higher resolution
tbtest = getEnrichmentTable(resttest)

tbtest_update = getEnrichmentTable(resttest) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(resttest,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      )
    )
  ) |>
  mutate(updated_pvalue_adj = p.adjust(updated_pvalue,method = "fdr"))


tb_deftest <- merge(x=tbtest_update,y=go_terms,by="id",all.x=TRUE)
tb_deftest <- tb_deftest[tb_deftest$fold_enrichment>2 & tb_deftest$p_adjust<0.1,]
write.table(tb_deftest,"great_all_cns_merged_50bp_go_enrichment.tsv",sep = '\t',row.names = FALSE,quote = FALSE)

reg_genetest <- getRegionGeneAssociations(resttest)
reg_gene_dftest <- data.frame(seqnames=seqnames(reg_genetest),
                              starts=start(reg_genetest)-1,
                              ends=end(reg_genetest)
)
reg_gene_dftest$annotated_genes <- data.frame(Genes = sapply(reg_genetest$annotated_genes, function(x) paste(x, collapse = ",")))$Genes
reg_gene_dftest$dist_to_TSS <- data.frame(Integers = sapply(reg_genetest$dist_to_TSS, function(x) paste(x, collapse = ",")))$Integers
write.table(reg_gene_dftest,"great_all_cns_merged_50bp_go_enrichment_region2geneassociation.tsv",sep = '\t',row.names = FALSE,quote = FALSE)

# Tripsacineae with Merged CNS within 50bp
#gr_obj =  import("data\\panand_cns_phylop_ACC_LRT_adjusted_signif_merged_50bp.bed")
#background = import("data\\panand_cns_50bpmerge.bed")
res = great(gr_obj, go_list,extended_tss = et, background = background)
tbtest = getEnrichmentTable(res) |>
  tibble::as_tibble() |>
  dplyr::mutate(
    updated_pvalue = exp(
      pbinom(
        q = observed_region_hits - 1,
        size = attr(res,which = "n_total"),
        prob = genome_fraction,
        lower.tail = F,
        log.p = T
      )
    )
  ) |>
  mutate(updated_pvalue_adj = p.adjust(updated_pvalue,method = "fdr"))

tb_deftest <- merge(x=tbtest,y=go_terms,by="id",all.x=TRUE)
tb_deftest <- tb_deftest[tb_deftest$fold_enrichment>2 & tb_deftest$p_adjust<0.1,]
write.table(tb_deftest,"great_trips_cns_merged_50bp_go_enrichment.tsv",sep = '\t',row.names = FALSE,quote = FALSE)

reg_genetest <- getRegionGeneAssociations(resttest)
reg_gene_dftest <- data.frame(seqnames=seqnames(reg_genetest),
                              starts=start(reg_genetest)-1,
                              ends=end(reg_genetest)
)
reg_gene_dftest$annotated_genes <- data.frame(Genes = sapply(reg_genetest$annotated_genes, function(x) paste(x, collapse = ",")))$Genes
reg_gene_dftest$dist_to_TSS <- data.frame(Integers = sapply(reg_genetest$dist_to_TSS, function(x) paste(x, collapse = ",")))$Integers
write.table(reg_gene_dftest,"great_trips_cns_merged_50bp_go_enrichment_region2geneassociation.tsv",sep = '\t',row.names = FALSE,quote = FALSE)

