#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate multiqc
index=$(basename $(pwd))
read1=$1
read2=$2
cpus=$SLURM_JOB_CPUS_PER_NODE
out=$(basename ${read1} | sed 's/_1.fq.gz//g')
STAR \
  --genomeDir $index \
  --runThreadN $cpus \
  --runMode alignReads \
  --readFilesIn $read1 $read2 \
  --readFilesCommand zcat \
  --outFileNamePrefix ${out}_ \
  --outSAMtype None

