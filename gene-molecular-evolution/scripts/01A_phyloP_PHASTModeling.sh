#!/bin/bash

## phyloP scan for panand assemblies
export PATH=/programs/phast/bin:$PATH
cd /workdir/sh2246/p_panAndOGASR/

# orthofinder cds alignment
# identify Paspalum reference
mkdir output/ref_og
cat data/filteredOrthogroups_final.csv |parallel -j 30 "grep '>Pavag' data/og_v9/og-CDS-aln_with-alt_v9.2b/{}_mafft.aln|sed 's/>//g'|shuf -n 1 > output/ref_og/{}.ref"

# Step 1: process alignment (remove gap)
mkdir output/orthoFinderOGMSA/
cat data/filteredOrthogroups_final.csv |parallel -j 38 "Rscript paml_pipeline/Ungap_MSA_SKH.R data/og_v9/og-CDS-aln_with-alt_v9.2b/{}_mafft.aln output/ref_og/{}.ref output/orthoFinderOGMSA/{}.gs.fa"

find output/orthoFinderOGMSA/ -name "*.gs.fa" -type f |parallel -j 30 "cat {} | seqtk seq -l0 > {= s/gs.fa/gs.l0.fa/  =}"

## exclude outgroup setaria 
mkdir output/orthoFinderOGMSA_noSevir/
find output/orthoFinderOGMSA/ -name "*.gs.l0.fa" -type f| parallel -j 30 "seqkit grep -p "Sevir*" -r -v {} -w 0 -o {= s/orthoFinderOGMSA/orthoFinderOGMSA_noSevir/  =}"


# Step 2: create dummy feature files for CDS alignment
mkdir output/featureFiles_og/
find output/orthoFinderOGMSA_noSevir/ -name "*.gs.l0.fa" -type f | parallel -j 30 "head -n 2 {} |tail -n 1 |wc -c| awk '{print \"dummy\tCDS\t\"1\"\t\"\$0-1\"\t.\t+\t.\"\"\tID = gene1\"}' > output/featureFiles_og/{/.}.tmp"

cat data/filteredOrthogroups_final.csv | parallel -j 30 'paste -d "\t" output/ref_og/{}.ref output/featureFiles_og/{}.gs.l0.tmp > output/featureFiles_og/{}.gff3'

# step 3: get 4d MSA
mkdir output/og_aln_4d
cat data/filteredOrthogroups_final.csv | parallel -j 30 'Rscript src/03_run_get4d.msa.R --feature output/featureFiles_og/{}.gff3 --input output/orthoFinderOGMSA_noSevir/{}.gs.l0.fa  --output output/og_aln_4d/{}.4d.fa' 

# step 4: construct RAxML tree
mkdir output/tree_4d_og/
cd output/tree_4d_og/
find /workdir/sh2246/p_panAndOGASR/output/og_aln_4d -type f | parallel -j 30 '/workdir/sh2246/p_panAndOGASR/paml_pipeline/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s {} -# 1 -n {/.}.tree'
cd /workdir/sh2246/p_panAndOGASR/


# step 5: fit neutral model
mkdir output/neutralModel4d_og
cat data/filteredOrthogroups_final.csv | parallel -j 40 'phyloFit --tree output/tree_4d_og/RAxML_bestTree.{}.4d.tree -i FASTA -I 10 --EM --out-root output/neutralModel4d_og/{} output/og_aln_4d/{}.4d.fa' 

# step 6: label tree
cat data/filteredOrthogroups_final.csv | parallel -j 30 'tree_doctor --name-ancestors output/neutralModel4d_og/{}.mod > output/neutralModel4d_og/{}.labeled.mod'

# step 7: identify branch of interest
mkdir output/target_og
cat data/filteredOrthogroups_final.csv |parallel -j 40 "grep 'Zm\|Zh\|Zl\|Zn\|Zv\|Zx\|Sobic\|Te' data/og_v9/og-CDS-aln_with-alt_v9.2b/{}_mafft.aln |sed 's/>//g' > output/target_og/{}.target"

mkdir output/target_phast_og
cat data/filteredOrthogroups_final.csv | parallel -j 30 'Rscript src/LabelNodes_phast.R output/tree_4d_og/RAxML_bestTree.{}.4d.tree output/target_og/{}.target output/target_phast_og/{}'

# step 8: phyloP test
mkdir output/phyloP_og
cat data/filteredOrthogroups_final.csv | parallel -j 30 'phyloP --branch `cat output/target_phast_og/{}` --method LRT --mode CONACC -i FASTA --feature output/featureFiles_og/{}.gff3 output/neutralModel4d_og/{}.labeled.mod output/orthoFinderOGMSA_noSevir/{}.gs.l0.fa > output/phyloP_og/{}.feature'

cat output/phyloP_og/* |grep -v "#" > output/phyloP_res_og.txt
sed -i 's/0.00000/0.00001/g' output/phyloP_res_og.txt
