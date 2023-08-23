#!/bin/bash
ml purge
ml singularity
scmd="singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/trinity.sif"
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus="${cpus:-36}"
echo $cpus
bam="$1"
out=$(basename ${bam%.*})
$scmd Trinity \
   --genome_guided_bam ${bam} \
  --include_supertranscripts \
   --max_memory 300G \
   --min_contig_length 100 \
   --genome_guided_max_intron 10000 \
   --full_cleanup \
   --CPU $cpus
