#!/bin/bash
SIF="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//ltr_retriever_latest.sif"
module purge
ml singularity
genome=$1
intact=$(find $(pwd) -name "*fasta.mod.pass.list")
anno=$(find $(pwd) -name "*fasta.mod.EDTA.TEanno.out")
singularity exec --bind $PWD $SIF /LTR_retriever/LAI -genome $genome -intact $intact -all $anno -t 36

