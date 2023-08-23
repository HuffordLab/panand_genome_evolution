#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate stringtie
bam="$1"
out=$(basename ${bam%.*})
cpus=$SLURM_JOB_CPUS_PER_NODE
stringtie \
   ${bam} \
   --rf \
   -m 100 \
   -p $cpus \
   -v \
   -o ${out}_stringtie.gtf
