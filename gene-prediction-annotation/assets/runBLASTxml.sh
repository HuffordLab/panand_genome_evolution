#!/bin/bash
module load blast-plus
input=$1
cwd=$(pwd)
procs=4
blastx -max_target_seqs 5 -num_threads $procs -query ${input} -outfmt 5 -db uniprot-sprot_viridiplantae.fasta -evalue 0.000001 2> ${TMPDIR}/${input%.*}.log | sed '/^$/d' | gzip -c - > ${TMPDIR}/${input%.*}.blast.xml.gz
mv ${TMPDIR}/${input%.*}.blast.xml.gz $cwd/
