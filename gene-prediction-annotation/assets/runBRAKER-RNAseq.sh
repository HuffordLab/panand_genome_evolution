#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
source /work/LAS/mhufford-lab/arnstrm/programs/sourceme
#cp /work/LAS/mhufford-lab/arnstrm/programs/gm_key_64 ~/.gm_key
conda activate braker
genome=$1
bam=$2
cpus=$SLURM_JOB_CPUS_PER_NODE
profile="${genome%.*}_$(date +"%Y%m%d")"
echo "profile generated for this genome is $profile"
GENEMARK_PATH="/work/LAS/mhufford-lab/arnstrm/programs/genemark-es-4.65"
braker.pl \
   --genome=${genome} \
   --bam=${bam} \
   --softmasking \
   --species=${profile} \
   --cores $cpus
