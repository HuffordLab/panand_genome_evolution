#!/bin/bash
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
genome=$(find . -maxdepth 1 -name "*.fasta")
index=$(basename $(pwd))
echo "${scripts}/runSTARmap_I.sh $genome" > ${index}.cmds
for f in *_1.fq.gz; do echo "${scripts}/runSTARmap_R1.sh $f $(echo $f |sed 's/_1.fq.gz/_2.fq.gz/g')"; done >> ${index}.cmds
echo "awk -f ${scripts}/sjCollapseSamples.awk *_SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > SJ.all" >> ${index}.cmds
for f in *_1.fq.gz; do echo "${scripts}/runSTARmap_R2.sh $f $(echo $f |sed 's/_1.fq.gz/_2.fq.gz/g')"; done >> ${index}.cmds
echo "${scripts}/mergeBAM.sh" >> ${index}.cmds
makeSLURMs.py 100 ${index}.cmds

