#!/bin/bash
dir=$1
cd ${dir}
cwd=$(pwd)
genome=$(find ${cwd}/00_raw-data/genome -name "${dir}_v*_p.fasta")
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
cd ${cwd}/04_transcript-assembly/trinity
trinity=$(find $(pwd)/trinity_out_dir -maxdepth 1 -name "Trinity-GG.fasta")
mkdir -p ${cwd}/04_transcript-assembly/trinity_gmap
cd "${cwd}/04_transcript-assembly/trinity_gmap"
index="GMAP_${dir}_p"
ln -s ${genome}
ln -s ${trinity}
cat > ${name}-gmap.sub <<'Endofmessage1'
#!/bin/bash
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH -p amd
#SBATCH --time=4-00:00:00
#SBATCH --account=las
#SBATCH --qos=las
#SBATCH --job-name=JOBNAME
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
cwd=$SLURM_SUBMIT_DIR
cpus=$SLURM_JOB_CPUS_PER_NODE
ulimit -s unlimited
ml purge
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate multiqc
Endofmessage1
cat >> ${name}-gmap.sub <<Endofmessage2
name=${name}
trinity=$(basename ${trinity})
dir=${dir}
index=${index}
dbpath=$(pwd)
genome=$(basename ${genome})
Endofmessage2
cat >> ${name}-gmap.sub <<'Endofmessage3'
gmap_build \
   -d ${index} \
   -D ${dbpath} \
      $(basename $genome)
gmapl.sse42 \
   -D ${dbpath} \
   -d ${index} \
   -B 4 \
   -t ${cpus} \
   -f gff3_match_cdna \
      ${trinity} \
    > ${name}_TrinityGG-mapped.gff3  \
    2> ${name}_TrinityGG-mapped.stderr-out || {
echo >&2 "gmapl is not suitable for $dir, trying regular gmap now.."
gmap.sse42 \
   -D ${dbpath} \
   -d ${index} \
   -B 4 \
   -t ${cpus} \
   -f gff3_match_cdna \
      ${trinity} \
    > ${name}_TrinityGG-mapped.gff3 \
    2> ${name}_TrinityGG-mapped.stderr-out
}
Endofmessage3
pwd
sed -i 's/JOBNAME/'${name}'-gmap/g' ${name}-gmap.sub
sbatch  ${name}-gmap.sub




