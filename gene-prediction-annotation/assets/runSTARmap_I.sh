#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate multiqc
genome=$1
cpus=$SLURM_JOB_CPUS_PER_NODE
index=$(basename $(pwd))
mkdir -p $index
STAR \
  --runMode genomeGenerate \
  --runThreadN $cpus \
  --genomeDir $index \
  --genomeFastaFiles ${genome}

