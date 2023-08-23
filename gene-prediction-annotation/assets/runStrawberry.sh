#!/bin/bash
#source /work/LAS/mhufford-lab/arnstrm/programs/strawberry_1.1.1/sourceme
#module use /work/GIF/software/modules
#module load GIF/strawberry/1.1.1
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/programs/strawberry_1.1.1
bam="$1"
cpus=$SLURM_JOB_CPUS_PER_NODE
strawberry \
   --output-gtf ${bam%.*}_strawberry.gtf \
   --logfile strawberry_assembly.log \
   --no-quant \
   --num-threads $cpus \
   --verbose \
   --fr \
   --min-transcript-size 100 \
     ${bam}
