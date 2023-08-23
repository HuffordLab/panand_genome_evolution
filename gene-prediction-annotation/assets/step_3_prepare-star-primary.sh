#!/bin/bash
dir=$1
cd ${dir}
cwd=$(pwd)
genome="$cwd/00_raw-data/genome/${dir}_v1_p.fasta"
mv 00_raw-data/fastq-files/*.gz ./
ln -s $genome
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
${scripts}/step_1_make-star-slurm.sh
