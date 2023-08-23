#!/bin/bash
# process mikado gff3 files to procuce a clean gff3 file
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus=${cpus:-36}
# cleaned up gff
awk '$3=="mRNA" {split($9,a,";"); gsub("ID=","",a[1]); gsub("Parent=","",a[2]); print a[1]"\t"a[2]}' mikado.loci.gff3 > mikado-mRNA-gene.ids
${scmd} mikado util grep mikado-mRNA-gene.ids mikado.loci.gff3 > mikado.genes.gff3
#convert to gtf
${scmd} mikado util convert -of gtf mikado.genes.gff3 mikado.genes.gtf

