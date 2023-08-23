#!/bin/bash
cpus=${SLURM_JOB_CPUS_PER_NODE:-36}
genome=$1
bam=$2
transcripts=$3
proteins=$4
stringtie=$5
dir="$6"
FUNANNOTATE_DB="/work/LAS/mhufford-lab/arnstrm/DATABASES/FUNANNOTATE"
ml purge
species=$(echo $dir | sed 's/-/ /g')
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate funannotate
funannotate predict \
    --input ${genome} \
    --out funannotate_out \
    --species "${species}" \
    --transcript_evidence $transcripts \
    --protein_evidence ${proteins} $FUNANNOTATE_DB/uniprot_sprot.fasta \
    --rna_bam ${bam} \
    --stringtie ${stringtie} \
    --busco_seed_species maize \
    --optimize_augustus \
    --organism other \
    --busco_db embryophyta \
    --max_intronlen 50000 \
    --force \
    --cpus ${cpus}

