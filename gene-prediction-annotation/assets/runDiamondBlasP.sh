#!/bin/bash
ml purge
#ml diamond
cpus=$SLURM_JOB_CPUS_PER_NODE
qid=NCBITAXAID
tid=$1
diamond blastp --threads $cpus --db $tid --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score --out $tid.out --query $qid.faa
echo -e "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tscore\tstaxid" > ${tid}.tab
awk -v x=${tid} 'BEGIN{OFS=FS="\t"}{print $1,$2,$7,$8,$9, $10, $11, $13, x}' $tid.out >> ${tid}.tab
rm $tid.out
