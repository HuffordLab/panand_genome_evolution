## Calculating TFBS turnover rate between maize and other Andropogoneae

The overall approach to calculating the turnover of TFBS is as follows.

1. Define target genomic regions of 1kb upstream of translation start sites 
2. Predict TFBS locations de novo in all genomes
3. For each Z.mays gene target region, collect multiple alignments of TFBS as well as de novo annotated TFBS for all other Andropogoneae species
4. Motifs with a 100% sequence match without gaps to Z. mays and/or a denovo prediction of the same TFBS as Z.mays are defined as high-quality matching motifs
5. Target genomic regions for genes are retained only for genes that are: syntenic, core genes (Hufford et al. 2021), expressed at >0 TPM (Hufford et al., 2021), have a ChIP-seq peak (Tu et al., 2020) in their 1kb target region. The gene list is in `data/core_expr_chip_syn_genes.txt.gz`. 
6. The number of matching motifs and unmatched motifs is then counted for each species 

### Annotating TFBS

Plant consensus motifs from the JASPAR 2022 database were downloaded and truncated to the core motifs. Motifs were then predicted for each Andropogoneae genome using a custom kotlin script.

### Analysing TFBS

Sequences aligned to the maize motifs were extracted from the multiple alignment MAF file to generate a BED file as shown below. Note that the `species_list`,`motif_sequences`, and `motif_coordinates` are in the same order.
```
#chr,start,end,motif_id,motif_len,chr_gene_region,start_gene_region,end_gene_region,gene_transcript,strand,NA,species_list,motif_sequences,motif_coordinates
1       1002913 1002921 MA1232.1        8       1       1002403 1003403 Zm00001eb000210_T001    +       1       consensus,zB73v5,agerar@3,udigit@1,agerar@1,hconto@2,zTIL18,snutan@1,rrottb@1,zTIL11,smicro@1,sscopa@1,zluxur,zdiplm,sscopa@2,hcompr@1,ppanic@1,blagur@3,telega@2,ccitra@1,vcuspi@1,udigit@3,zhuehu,achine@1,hconto@4,hcompr@4,crefra@1,vcuspi@2,zTIL01,aburma@1,sbicol@1,ccitra@2,zTIL25,irugos@1,avirgi@1,hconto@3,hcompr@5,znicar,udigit@2,etrips@1,snutan@2,etrips@2,rrottb@2,blagur@2,zdiplg,blagur@1,rtuber@1,hcompr@6,telega@1,aburma@2,achine@2,hcompr@3,hconto@1,hcompr@2,ttrian@1,tdactn,agerar@2,tdacts      TCACCGTC,TCACCGTC,TCACCGTC,TCACCCTC,TCACCGTC,--------,TCACCGTC,--GTTGTC,--------,TCACCGTC,--------,TCATCGTC,TCACCGTC,TCACCGTC,--------,TC------,TCACCGTC,TCACCGTC,TCACCGTC,TCACCG--,TCACCGTC,TCACCCTC,TCACCGTC,TCAGCGTC,--------,TC------,TCACCGTC,TCACCGTC,TCACCGTC,TCACCGTC,TCACCGTC,TCACCGTC,TCACCGTC,TTACCGTC,TCATCGTC,--------,TC------,--------,TCACCCTC,TCACCGTC,--GTTGTC,TCACCGTC,--------,--------,TCACCGTC,TCACCGTC,TCACCGTC,TC------,TCACCGTC,TCACCGTC,TCAGCGTC,TC------,--------,--------,--------,TCACCGTC,TCACCGTC,TCACCGTC       NA,1@1002909@1002921@+,scaf_7@219379@219391@+,scaf_1@241706115@241706127@+,scaf_58@36142789@36142801@-,NA,chr1@1406813@1406825@+,scaf_6@71799671@71799677@-,NA,chr1@1405137@1405149@+,NA,scaf_39@177417@177429@+,chr1@302728@302740@+,chr1@2049770@2049782@+,NA,scaf_42@55777937@55777940@-,ctg_110892@379584@379596@-,scaf_30@1313366@1313378@+,ctg_61302@659278@659290@-,ctg_363678@1342655@1342665@+,scaf_5@5622113@5622125@-,scaf_7@157157981@157157993@-,chr1@687588@687600@+,ctg_309666@92027@92039@+,NA,scaf_26@69081178@69081181@-,scaf_1@52039454@52039466@-,scaf_113@227185@227197@-,chr1@1449621@1449633@+,ctg_354967@482836@482848@+,Chr01@80748100@80748112@-,ctg_8799@879718@879730@+,chr1@1501623@1501635@+,ctg_38502@315915@315927@+,chr2@112457260@112457272@-,NA,scaf_29@66789820@66789823@-,NA,scaf_2@931902@931914@-,scaf_5@85977@85989@+,scaf_1@115314@115320@+,scaf_3@178607681@178607693@-,NA,NA,chr1@1169519@1169531@+,scaf_8@60880@60892@+,ctg_84888@19361003@19361012@-,scaf_31@65508167@65508170@-,ctg_310389@876556@876568@-,ctg_119803@172199@172211@+,ctg_206679@294256@294268@+,scaf_32@360029@360032@+,NA,NA,NA,scaf_3@131213302@131213314@-,scaf_8@435670@435682@+,scaf_6@17106825@17106837@+
```

Mismatches and gaps between the maize reference TFBS sequence and each species' aligned sequence were then counted to generate a file classifying each TFBS alignment.

```
python count_conserved_tfbs.py tfbs_aligned.txt > tfbs_seq_comparison.txt
cut -f1-4,9,11,12- tfbs_seq_comparison.txt| awk '{print $1"_"$2"_"$3"_"$5".fa",$6,$9,$10,$7,$8,$4,$11}'|tr ' ' '\t' > tfbs_aligned_class.txt
```

In the next key step, all alignment and denovo TFBS prediction information `panand_1000bp_TSS_regions_JASPARmotifs.bed` for each site is integrated to rescue any conserved TFBS that could not be aligned. We also ensure that the compared motifs are in the putative regulatory regions of collinear syntenic orthologs using the reference information in `collinear_orthologous_syntenic_genes_panand.txt`.

```
python rescue_tfbs_with_denovo_alignment_w_synteny_input_withcoords.py maize_tfbs.bed panand_tfbs.bed tfbs_aligned_class.txt collinear_orthologous_syntenic_genes_panand.txt > tfbs_aligned_class_denovo.txt
```

Next, genes and their assigned TFBS are excluded if they are dispensable (rather than core) genes within maize, are not expressed (TPM=0), have no ChIP-seq signal in their putative regulatory region, or are not syntenic.

```
python filter_genes.py data/core_expr_chip_syn_genes.txt tfbs_aligned_class_denovo.txt > tfbs_aligned_class_denovo_filt.txt
```

A single edge case led to a collinear gene having no denovo TFBS, so this was manually set to a value of "NA" rather than being left blank.
```
awk '$11!=""' tfbs_aligned_class_denovo_filt.txt > tmp_file_without_null.txt
awk '$11==""' tfbs_aligned_class_denovo_filt.txt | sed 's/$/NA/' > tmp_file_with_null.txt
cat  tmp_file_without_null.txt tmp_file_with_null.txt > tfbs_aligned_class_denovo_filt.txt
```

Finally, turnover rates were calculated for each species using a custom python script.
```
python calculate_turnover.py  tfbs_aligned_class_denovo_filt.txt > tfbs_aligned_class_denovo_filt_turnover.txt 
```

The result is summarized in `output/tfbs_turnover_statistics.csv`.
