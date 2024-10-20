
## Defining collinear syntenic genes across the PanAnd group

To ensure comparisons of TFBS turnover between PanAnd species are comparing orthologous gene regulatory regions, a set of target genes is defined based on the following criteria:

* Cluster in same orthogroup
* Are collinear
* Are syntenic

A custom script uses all of the data from files in `data/` to output this set of genes.

```
python filter_syntenic_gene_regions.py >collinear_orthologous_syntenic_genes_panand.txt
```

Tripsacum species `tdactn` and `tdacts` are special cases since their assembly and annotation was redone during this analysis, so orthology and collinearity information for these species is drawn from an earlier separate set of files.

The resulting set of genes can now be used to filter the TFBS in gene regulatory regions.

## Get statistics on single copy Z. mays gene turnover per subgenome

To compare the TFBS turnover to gene turnover, single copy genes were also assessed for turnover between species.
```
python calc_gene_turnover.py  > singlecopy_gene_subgenome_orthologs.txt
python count_turnover.py singlecopy_gene_subgenome_orthologs.txt > singlecopy_gene_subgenome_orthologs_stats.txt
```

In total 5951 single copy genes that overlapped the ~11k TFBS genes were assessed.

```
python count_turnover.py singlecopy_gene_subgenome_orthologs.txt| tr ' ' '\t'| grep -v unknown| grep -v maize2| egrep -v "maize1a|maize1b"| sort -k3,3n| awk -v OFS='\t' '{print $2,$1,$3,"5951",1-$3/5951}'
```

Because some subgenomes are incomplete, it makes more sense to compare the entire genome to check for gene loss:

```
cut -f2 -d' ' singlecopy_gene_subgenome_orthologs.txt| sort| uniq| while read l; do echo $l; grep $l singlecopy_gene_subgenome_orthologs.txt | cut -d' ' -f1| sort| uniq| wc -l;done| paste - -| sort -k2,2n| awk -v OFS='\t' '{print $1,$2,1-$2/5951}' > singlecopy_gene_species_orthologs_rate.txt
```

