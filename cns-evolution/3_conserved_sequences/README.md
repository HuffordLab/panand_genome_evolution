## Calling conserved regions from the PanAnd alignment

A parameter search was conducted by varying the target coverage of conserved region in the alignment from 0.1-0.5 in increments of 0.05. For example, a model was fit given a target coverage of 0.45 using a neutral model based on 4d sites and outputting a model with branch lengths under a conserved model fitting the target coverage. 

```
phastCons --target-coverage 0.45 chr1.maf data/models/chr10_4d_single_nozea.mod --no-post-probs --estimate-rho chr1_cov0.45
```

Next conserved regions were called using all of the pre-estimated models with coverages 0.1-0.5, testing expected lengths of conserved regions of 8,12,20, and 40. 

```
phastCons --most-conserved chr1_cov0.45_len8.most-cons.bed --target-coverage 0.50 --expected-length 8 chr1.maf data/models/chr1_4dsingle_cov0.45_single.cons.mod,data/models/chr1_4dsingle_cov0.45_single.noncons.mod  > chr1_cov0.45_len8.wig
```

The resulting conserved regions are stored under `data/parametrization/`. The optimal parametrization was selected to maximize the enrichment of accessible chromatin based on regions made available as GEO_ID:GSE155178 by Marand et al., 2021 (doi: 10.1016/j.cell.2021.04.014), which led to an optimal expected length of 8 and target coverage of 0.45.

## Filtering conserved regions

Raw conserved regions were then filtered to generate a set of high-quality conserved non-coding elements using the following criteria:

* Regions intersecting B73 NAMv5 "CDS" features were removed
* Remaining regions were retained if they had a minimum length of 5bp
* Elements with a blastx match to any protein in the Uniprot TREMBL database with an e-value threshold of 0.01 were removed

The remaining 1,664,343 elements were considered high-quality CNS and used for further analyses. 

## Calculating enriched features in conserved regions

Enrichments of sequence features were calculated based on the ratio of proportion of the base pairs of a feature within the CNS and within the whole-genome background. If a feature covers 1% of the genome and 10% of the CNS, then this feature is 10-fold enriched. Gene-related features were extracted from the public `Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1` annotation. Additional features assessed were:

* Chromatin loops (GEO ID:GSE120304) 
* Accessible chromatin regions (GEO ID: GSE155178)
* MOA motifs from Savadel et al, 2021 found in file "S8 Zip File MOA-seq Motifs BED files for B73v5 mapped via BLAST-over from B73v3" (doi: 10.1371/journal.pgen.1009689)
* ChIP-seq peaks from Supplementary Table 5 in Tu et al., 2020 (doi: 10.1038/s41467-020-18832-8) 

All BED features in these files were merged to ensure no features overlapped before calculating enrichments.

## Counting conserved sequences for each species

Conserved regions are not always present in every single species included in the alignment. To therefore count the number of CNS per species, each CNS was extracted from the alignment, converted to FASTA format, and then checked using a custom python script `cns_maf2mfa.py`. MAF regions overlapping CNS were extracted using `mafsInRegion` (available from: [https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)). A simple bash wrapper performs the operation, given a set of MAF alignments named by chromosome or scaffold. 

```
get_sp_seq_from_region.sh chr1_34713_34721
```

By intersecting CNS with gene annotations and using `bedtools merge` to collapse the result and then `bedtools closest -D b -t first` to get the distance to the nearest gene, it was possible to generate a BED file with the following information:
```
chr1    34713   34721   chromosome,exon,five_prime_UTR,gene,mRNA        chr1    34616   40204   Zm00001eb000010       .       +       0
```
Here column four is a list of all the unique features that intersect the CNS element while column 11 is the distance to the nearest gene with negative values denoting that the feature is upstream of the gene (<0) and positive values that it is downstream (>0) of or within (=0) the gene. This file was used as input to the custom script `classify_cns.py` to assign a single feature (`downstream|downstream_distal|five_prime_UTR|intron|three_prime_UTR|upstream|upstream_distal`) to each CNS.

