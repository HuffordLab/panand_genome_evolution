#!/bin/bash
gff=${gff:-"mikado.genes.gff3"}
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus=${cpus:-36}
ml purge
awk '$NF>0 {print $1}' mikado.transcripts-counts > mikado.transcripts-to-reatain.txt
grep -Fw -f mikado.transcripts-to-reatain.txt mikado-mRNA-gene.ids |cut -f 2 > mikado.genes-to-reatain.txt
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado2
mikado util grep --genes mikado.genes-to-reatain.txt mikado.genes.gff3 > mikado.genes-filtered.gff3
