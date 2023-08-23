#!/bin/bash
# create dir and copy required files
mikado=$(find $(pwd) -name "mikado_prepared.fasta")
cd $(dirname ${mikado})
mkdir -p blastjobs
cd blastjobs
ln -s ${mikado}
for db in /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/uniprot-sprot_viridiplantae.fasta*; do
ln -s $db;
done
# spilit input
fasta-splitter.pl --n-parts 8 mikado_prepared.fasta
unlink mikado_prepared.fasta
# blast in parallel
ml purge
ml parallel
ml blast-plus
parallel \
  -j 8 \
  --joblog blast_progress.log \
  --workdir $PWD \
   "/work/triffid/arnstrm/PanAnd/current-versions-genomes/runBLASTxml.sh {}" ::: mikado_prepared.part-?.fasta
# extract XML
for xml in *xml.gz; do
   gunzip $xml;
done
# merge xml
ml purge
python2 ~/gitdirs/common_analyses/BlastXMLmerge.py ../mikado.blast.xml *.xml
# compress XML
for xml in *xml; do
   gzip $xml;
done
