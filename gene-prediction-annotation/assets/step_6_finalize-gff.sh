#!/bin/bash
# should be in 06_annotation/finalFiles folder
dir=$1
version=1.0
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
# gff files
primary="${cwd}/00_raw-data/genome/${dir}_*_p.fasta"
bind=$(find ${cwd}/04_transcript-assembly/gemoma/gemoma_output_v1 -name "final_annotation.gff")

mkdir -p 06_annotation/finalFiles
cd 06_annotation/finalFiles
ml purge
ml genometools
gt gff3 -sort -tidy -o ${name}_BIND.v${version}.gff3 $bind &> ${name}.stdout;
ml purge
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate
agat_sp_statistics.pl --gff ${name}_BIND.v${version}.gff3 -o ${name}_BIND.v${version}.stats
ml purge
ml cufflinks
gffread ${name}_BIND.v${version}.gff3 -g ${primary} -x ${name}_cds.v${version}.fasta -y ${name}_pep.v${version}.fasta
ml purge
ml bioawk
for f in ${name}_cds.v${version}.fasta; do
bioawk -c fastx 'BEGIN{OFS="\t"}{print $name,length($seq),gc($seq),substr($seq,0,3),reverse(substr(reverse($seq),0,3))}' $f > ${f%.*}_len-gc-start-stop.tsv;
done
