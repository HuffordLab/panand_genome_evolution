#!/bin/bash
ml purge
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate class
bam="$1"
out=$(basename ${bam%.*})
cpus=$SLURM_JOB_CPUS_PER_NODE
psiclass -b $bam \
         -o ${out}_psiclass.gtf \
         -p $cpus \
        --primaryParalog
