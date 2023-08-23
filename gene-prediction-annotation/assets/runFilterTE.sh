#!/bin/bash
teTSV=$2
gff=$1
genome=$3
name=$4
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado
awk '$1 !~ /^#/ {print $1}' ${teTSV} |sort | uniq > te.ids
mikado util grep -v te.ids $gff > ${gff%.*}_noTE.gff
ml purge
ml cufflinks
gffread ${gff%.*}_noTE.gff -g ${genome} -x ${name}_cds.fasta -y ${name}_pep.fasta
