#!/bin/bash
cds=$1
module purge
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate edta
ml hmmer blast-plus

cpus=$SLURM_JOB_CPUS_PER_NODE
mkdir -p tmp
TEsorter -p $cpus -db rexdb-plant --tmp-dir tmp $cds

