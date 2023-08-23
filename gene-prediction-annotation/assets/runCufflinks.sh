#!/bin/bash
module load cufflinks
bam="$1"
out=$(basename ${bam%.*})
mkdir -p ${out}
cpus=$SLURM_JOB_CPUS_PER_NODE
cufflinks \
   --output-dir "${out}" \
   --num-threads $cpus \
   --frag-len-mean 100 \
   --library-type fr-firststrand \
   --multi-read-correct \
   --verbose \
   --no-update-check \
   ${bam}

