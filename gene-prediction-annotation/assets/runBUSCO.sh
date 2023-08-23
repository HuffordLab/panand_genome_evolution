#!/bin/bash
ml purge
cwd=$(pwd)
fasta=$1
cpus="$SLURM_JOB_CPUS_PER_NODE"
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate busco5
cd /local/scratch/arnstrm/${SLURM_JOB_ID}
rsync -avP /work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2/busco_downloads.tar.gz .
rsync -avPL ${cwd}/${fasta} ./
tar xf busco_downloads.tar.gz
for spp in poales_odb10 liliopsida_odb10; do
out=$(basename ${fasta%.*})
busco \
   -i ${fasta} \
   -c ${cpus} \
   -o ${out}_${spp} \
   -m genome \
   -l ${spp} \
   --offline \
   -f
cp ${out}_${spp}/short* ./
tar -cvzf ${out}_${spp}.tar.gz  ${out}_${spp}/
rsync -avP ${out}_${spp}.tar.gz ${cwd}/
done
