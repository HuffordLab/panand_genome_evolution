#!/bin/bash
dir=$1
cd $dir
# variables
cwd=$(pwd)
rundir="mikado_v2.3"
container="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//mikado-v2.3.3.sif"
scmd="singularity exec --bind $PWD $container"
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
genome="${cwd}/00_raw-data/genome/${dir}*_p.fasta"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
bam=$(find $(pwd)/03_STAR -name "merged*.bam")
stringtie=$(find "${cwd}/04_transcript-assembly/stringtie" -name "merged*stringtie.gtf")
class2=$(find "${cwd}/04_transcript-assembly/class2" -name "*_psiclass.gtf_vote.gtf")
cufflinks=$(find "${cwd}/04_transcript-assembly/cufflinks" -name "transcripts.gtf")
trinity=$(find "${cwd}/04_transcript-assembly/trinity_gmap"  -name "*TrinityGG-mapped.gff3")
strawberry=$(find "${cwd}/04_transcript-assembly/strawberry" -name "*_strawberry.gtf")
portcullis=$(find "${cwd}/04_transcript-assembly/portcullis" -name "portcullis_filtered.pass.junctions.bed")
# directory for running
mkdir -p ${cwd}/04_transcript-assembly/${rundir}
cd ${cwd}/04_transcript-assembly/${rundir}
# softlink files needed
ln -s $class2 class.gtf
ln -s $cufflinks cufflinks.gtf
ln -s $strawberry strawberry.gtf
ln -s $stringtie stringtie.gtf
ln -s $trinity trinity.gff3
ln -s $portcullis junctions.bed
ln -s $(realpath $genome) $(basename $genome)
ln -s $(realpath $bam) $(basename $bam)
# create the list file
read -r -d '' list<<"Endofmessage1"
class.gtf\tcl\tTrue\t1
cufflinks.gtf\tcf\tTrue\t1
strawberry.gtf\tsb\tTrue\t1
stringtie.gtf\tst\tTrue\t1
trinity.gff3\ttr\tFalse\t-0.5
Endofmessage1
echo -e "$list" > list.txt
# create job file
cat > ${name}-mikado.sub <<"Endofmessage2"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH -n 36
#SBATCH --time=2-00:00:00
#SBATCH --account=las
#SBATCH --qos=las
#SBATCH --job-name=JOBNAME
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
ulimit -s unlimited
cd $SLURM_SUBMIT_DIR
cwd=$SLURM_SUBMIT_DIR
cpus=$SLURM_JOB_CPUS_PER_NODE
Endofmessage2
cat >> ${name}-mikado.sub <<Endofmessage3
scmd="${scmd}"
genome="$(basename $genome)"
bam="$(basename $bam)"
name="${name}"
scripts="${scripts}"
list="list.txt"
Endofmessage3
cat >> ${name}-mikado.sub <<"Endofmessage4"
ml purge
ml singularity
if [ ! -f "all-steps.done" ]; then
if [ ! -f "configure.done" ]; then
${scmd} ${scripts}/runMikado.1-configure.sh $list $genome || {
echo >&2 "ERROR: Step 1: Configure failed, exiting."
exit 1
}
fi
touch configure.done
if [ ! -f "prepare.done" ]; then
${scmd} ${scripts}/runMikado.2-prepare.sh || {
echo >&2 "ERROR: Step 2 Prepare failed, exiting."
exit 1
}
touch prepare.done
fi
if [ ! -f "blast.done" ]; then
${scripts}/runMikado.3-blast.sh || {
echo >&2 "ERROR: Step 3 BLASTx failed, exiting."
exit 1
}
touch blast.done
fi
if [ ! -f "transdecoder.done" ]; then
${scripts}/runMikado.4-transdecoder.sh || {
echo >&2 "ERROR: Step 4 TransDecoder failed, exiting."
exit 1
}
touch transdecoder.done
fi
if [ ! -f "serialise.done" ]; then
${scmd} ${scripts}/runMikado.5-serialise.sh || {
echo >&2 "ERROR: Step 5 Serialise failed, exiting."
exit 1
}
touch serialise.done
fi
if [ ! -f "pick.done" ]; then
${scmd} ${scripts}/runMikado.6-pick.sh || {
echo >&2 "ERROR: Step 6 Pick failed, exiting."
exit 1
}
touch pick.done
fi
if [ ! -f "mikado.done" ]; then
${scmd} ${scripts}/runMikado.7-process.sh|| {
echo >&2 "ERROR: Step 7 Process failed, exiting."
exit 1
}
touch mikado.done
fi
if [ ! -f "seq-extract.done" ]; then
${scripts}/runMikado.8-extract.sh $genome || {
echo >&2 "ERROR: Step 8 Extract failed, exiting."
exit 1
}
touch seq-extract.done
fi
if [ ! -f "featurecounts.done" ]; then
${scripts}/runMikado.9-feature_counts.sh ${bam} || {
echo >&2 "ERROR: Step 9 feature counts failed, exiting."
exit 1
}
touch featurecounts.done
fi
if [ ! -f "exp-filtering.done" ]; then
${scripts}/runMikado.10-expression_filter.sh || {
echo >&2 "ERROR: Step 10 filtering failed, exiting."
exit 1
}
touch exp-filtering.done
fi
touch all-steps.done
fi
Endofmessage4
sed -i 's/JOBNAME/'$name'_mikado/g' ${name}-mikado.sub
#sbatch ${name}-mikado.sub
