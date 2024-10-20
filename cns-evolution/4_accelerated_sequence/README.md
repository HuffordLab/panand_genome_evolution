## Detecting accelerated sequences in the Tripsacineae lineage

All CNS on B73 chromosomes were used to test for acceleration on the Tripsacineae branch with the LRT method implemented in phyloP. It was important that the discovery of conserved elements was conducted with the reference and all Tripsacineae species masked, to prevent a reference bias.

```
phyloP --features 1 --method LRT --branch branch_tripsacineae --gff-scores --mode ACC data/allchr.4d.gt05.mod chr1.maf > chr1_phylop_acc_lrt_trips.gff
```

The resulting log p-values were converted to p-values and then corrected using the Benajmini-Hochberg FDR method with the two custom scripts `convert_pval.py` and `correct_fdr.R`, respectively. The result is stored in `results/panand_cns_phylop_ACC_LRT_adjusted.bed.gz`.

To test whether Tripsacineae-accelerated CNS or just the CNS in general was enriched for specific gene ontology features, an analysis was conducted with rGREAT. As the number of CNS was large and included many neighboring CNS, all CNS within 50bp were merged into a single element to avoid a bias towards more fragmented CNS.

```
cut -f1-3,6 panand_cns_phylop_ACC_LRT_adjusted.bed| bedtools merge -i - -o min -c 4 -d 50| awk '$4<0.05' > panand_cns_phylop_ACC_LRT_adjusted_signif_merged_50bp.bed
```

The rGREAT analysis relies on several pieces of input information:

* Transcription start sites of all B73 genes: `data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.tss.tsv.gz` 
* List of genes and corresponding GO terms: `data/Zm00001eb.1.fulldata_goList.go2gene.txt.gz`
* GO term lookup: `data/go_term_lookup.tsv.gz`

The analysis is conducted with an R script `scripts/rgreat_enrichment.R` and generates a list of enriched GO terms `results/great_trips_cns_merged_50bp_go_enrichment.tsv` and a list of regions associated with genes matched to the enriched GO terms `results/great_trips_cns_merged_50bp_go_enrichment_region2geneassociation.tsv`.

