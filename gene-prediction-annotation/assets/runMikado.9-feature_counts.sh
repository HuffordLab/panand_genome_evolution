#!/bin/bash
bam=$1
gff=$2
gff=${gff:-"mikado.genes.gff3"}
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus=${cpus:-36}
echo -e "\t\tbam = $bam\n\t\tGFF = $gff\n\t\tcpus = $cpus"
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate multiqc
mkdir -p tmp
ml purge
ml subread
featureCounts -T ${cpus} -a ${gff} -o mikado.transcripts-counts -t mRNA -g Name --tmpDir ./tmp ${bam}

