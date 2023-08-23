#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado2
bam=$1
ref=$2
cpus=$SLURM_JOB_CPUS_PER_NODE
portcullis full \
   --threads=${cpus} \
   --verbose \
   --use_csi \
   --strandedness=UNKNOWN \
   --orientation=UNKNOWN \
   ${ref} \
   ${bam}
