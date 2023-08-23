#!/bin/bash
dir=$1
cd ${dir}
# variables
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
genome="${cwd}/00_raw-data/genome/${dir}_v1_p.fasta"
bam=$(find $(pwd)/03_STAR -name "merged*.bam")
proteins=$(find ${cwd}/04_transcript-assembly/mikado_v2.3 -name "mikado.proteins.fasta")
transcripts=$(find ${cwd}/04_transcript-assembly/mikado_v2.3 -name "mikado.transcripts.fasta")
stringtie=$(find ${cwd}/04_transcript-assembly/stringtie -name "*_stringtie.gtf")
repeats=$(find ${cwd}/05_repeatmasker/primary -name "*_v1_p.fasta.out.gff")
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
# setup work dir and softlink required files
mkdir -p ${cwd}/04_transcript-assembly/funannotate
cd ${cwd}/04_transcript-assembly/funannotate
ln -s ${bam} rnaseq.bam
ln -s ${repeats} repeats.gff
ln -s ${genome} genome.fa
ln -s ${stringtie} stringtie.gtf
ln -s ${transcripts} transcripts.fa
ln -s ${proteins} proteins.fa
# create jobfile
cat > ${name}-funanot.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=14-00:00:00
#SBATCH --account=triffid
#SBATCH --qos=triffid
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
cat >> ${name}-funanot.sub <<Endofmessage2
dir="$dir"
genome=genome.fa
bam=rnaseq.bam
proteins=proteins.fa
transcripts=transcripts.fa
repeats=repeats.gff
scripts=$scripts
Endofmessage2
cat >> ${name}-funanot.sub <<"Endofmessage3"
${scripts}/runSoftmasker.sh ${genome} ${repeats}
${scripts}/runFunannotate.sh ${genome%.*}.softmasked.fasta ${bam} ${transcripts} ${proteins} ${stringtie} "${dir}"
scontrol show job $SLURM_JOB_ID
Endofmessage3
sed -i 's/JOBNAME/'$name'_fa/g' ${name}-funanot.sub
# submit job
sbatch ${name}-funanot.sub
