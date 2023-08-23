#!/bin/bash
ml purge
ml samtools
ls *.bam > bam.fofn
cpus=$SLURM_JOB_CPUS_PER_NODE
index=$(basename $(pwd))
samtools merge --threads $cpus -b bam.fofn merged_${index}.bam
