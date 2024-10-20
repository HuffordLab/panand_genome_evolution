## Identifying collinear genomic regions using gene annotations

Collinear regions were identified with MCScanX using public B73 NAMv5 gene annotations and novel annotations from the PanAnd species. All-vs-All protein BLAST results of canonical gene isoforms were generated using diamond.

`diamond blastp --threads 48 --sensitive --db Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical.cds.dmnd --evalue 1e-10 --max-target-seqs 5 --outfmt 6 --query all.faa` 

After formatting the BLAST output and the gene annotations into files with the prefix zB73v5_public, MCScanX could be used to generate a custom output text file of collinear genes.

`MCScanX zB73v5_public`

The resulting `.collinearity` file was processed to generated a bed file of collinear regions between genomes.

```
grep -A 1 -B 1 "## Alignment" zB73v5_public.collinearity| grep -v "\-\-"| tail -n +2| paste - - - | sed 's/ //g' | sed '$ d' > zB73v5_public.collinearity.blocks.txt
python blocks2regions.py zmays.gff zB73v5_public.collinearity.blocks.txt > panand_collinearity.txt
``` 

Similarly, a collinearity analysis using Paspalum vaginatum rather than B73 as the reference can be conducted, revealing collinearity to the ancestral subgenome within the polyploid scaffold-level or ancient polyploid Tripsacineae genomes. By formatting the resulting collinearity output, a custom python script can be used to assign pseudosubgenomes to the polyploid species.

```
cut -f1 paspalum_collinearity.bed | sed 's/\..*//'| sort| uniq| while read l; do grep $l paspalum_collinearity.bed | grep -v pvagin.scaf| sort -k4,4 -k5,5n > $l.bed;done
python cluster_syntenic_paths.py achine.bed 2 > achine.2n.subgenomes.bed
```

By concatenating all of the resulting BED files, we can generate the `panand_subgenomes.txt` file.
