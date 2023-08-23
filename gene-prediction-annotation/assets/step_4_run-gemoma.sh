#!/bin/bash
dir=$1
cd ${dir}
# variables
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2"
genome="${cwd}/00_raw-data/genome/${dir}_*_p.fasta"
bam=$(find $(pwd)/03_STAR -name "merged*.bam")
braker=$(find ${cwd}/04_transcript-assembly/braker/braker -maxdepth 1 -name "braker.gtf")
mikado=$(find ${cwd}/04_transcript-assembly/mikado_v2.3 -maxdepth 1 -name "mikado.genes.gff3")
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
# setup work dir and softlink required files
mkdir -p ${cwd}/04_transcript-assembly/gemoma
cd ${cwd}/04_transcript-assembly/gemoma
ln -s ${bam} rnaseq.bam
ln -s ${genome} genome.fa
ln -s ${braker} braker.gtf
ln -s ${mikado} mikado.gff3
# create jobfile
cat > ${name}-genmoma.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH -n 36
#SBATCH --time=1-00:00:00
#SBATCH --account=las
#SBATCH --qos=las
#SBATCH --job-name=JOBNAME
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
ulimit -s unlimited
ml purge
cd $SLURM_SUBMIT_DIR
cwd=$SLURM_SUBMIT_DIR
cpus=$SLURM_JOB_CPUS_PER_NODE
Endofmessage1
cat >> ${name}-genmoma.sub <<Endofmessage2
dir="$dir"
genome=genome.fa
bam=rnaseq.bam
scripts=$scripts
braker=braker.gtf
mikado=mikado.gff3
Endofmessage2
cat >> ${name}-genmoma.sub <<"Endofmessage3"
${scripts}/runGeMoMa.v1.sh ${genome} ${bam} ${braker} ${mikado} "${dir}"
${scripts}/runGeMoMa.v2.sh ${genome} ${bam} ${braker} ${mikado} "${dir}"
scontrol show job $SLURM_JOB_ID
Endofmessage3
sed -i 's/JOBNAME/'$name'_GMM/g' ${name}-genmoma.sub
# submit job
sbatch ${name}-genmoma.sub
