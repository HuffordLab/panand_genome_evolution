#!/bin/bash
# run transdecoder
# create dir and copy required files
mikado=$(find $(pwd) -name "mikado_prepared.fasta")
cd $(dirname ${mikado})
mkdir -p transdecoder
cd transdecoder
ln -s $mikado
# load modules
ml purge
ml r-seqlogo
ml r-devtools r-ggplot2
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado2
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus=${cpus:-36}
# run longorfs td
TransDecoder.LongOrfs \
   -t mikado_prepared.fasta
# run predict td
TransDecoder.Predict \
   -t mikado_prepared.fasta \
   --cpu $cpus

