#!/bin/bash
cpus=${SLURM_JOB_CPUS_PER_NODE:-36}
genome=$1
bam=$2
braker=$3
mikado=$4
dir="$5"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
sorghumFa="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//Sorghum-bicolor/Sbicolor_454_v3.0.1.fa"
sorghumGFF="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//Sorghum-bicolor/Sbicolor_454_v3.1.1.gene_exons.gff3"
maizeFa="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//Zea-mays/GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna"
maizeGFF="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//Zea-mays/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff"
# directory for assemblies
ml purge
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate gemoma
java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar ${scripts}/GeMoMa-1.8.jar CLI GeMoMaPipeline sc=false o=true p=true pc=true pgr=true threads=$cpus outdir="gemoma_output_v3" \
   r=MAPPED ERE.m=${bam} ERE.s=FR_FIRST_STRAND \
   t=$genome \
   i=sb a=${sorghumGFF} g=${sorghumFa} w=1 \
   ID=braker e=$braker weight=0.1 \
   ID=mikado e=$mikado weight=0.9 \
   GeMoMa.v=true \
   AnnotationFinalizer.u=YES AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=${name} \
#   i=zm a=${maizeGFF} g=${maizeFa} w=1 \
# GAF.f="start=='M' and stop=='*' and score/aa>=0.75 and (evidence>1 or tpc==0.5)" 

