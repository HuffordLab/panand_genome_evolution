#!/bin/bash
if [ "$#" -lt 3 ] ; then
echo "please provide:"
echo -e "\t\t(1) directory name"
echo -e "\t\t(2) genome version (eg: 0, 1, 2)"
echo -e "\t\t(3) annotation version (eg: 0.1, 1.0, 2.0)"
echo "";
echo "./step_7_run-edata-on-genome.sh <dir> <genome version> <annotation version>" ;
echo "";
exit 0;
fi



dir=$1
gver=$2
aver=$3
cd $dir
cwd=$(pwd)
gffdir=${cwd}/06_annotation/finalFiles
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2"
genome=$(find $(pwd)/00_raw-data/genome -name "*_v${gver}.fasta")
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"

# check if the step_6_finalize-gff.sh is run
if [[ -d "${gffdir}" ]]
then
cd ${gffdir}
if [[ -f ${name}_cds.v${aver}.fasta ]]
fcds=$(realpath ${name}_cds.v${aver}.fasta)
then
echo "GFF has been finalized, will use this CDS for running EDTA"
else
echo "no cds found"
exit 1
fi
else
echo "GFF has not been finalized, please run step_6_finalize-gff.sh first"
exit 1
fi
# create all files
ml purge
ml bioawk
mkdir -p ${cwd}/07_edta-run
cd ${cwd}/07_edta-run
cp ${fcds} cds.fasta
ln -s ${genome} $(basename ${genome})
cat > ${name}-edta.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=21-00:00:00
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
ml purge
Endofmessage1
cat >> ${name}-edta.sub <<Endofmessage2
genome=$(basename ${genome})
cds=cds.fasta
scripts=$scripts
Endofmessage2
cat >> ${name}-edta.sub <<"Endofmessage3"
${scripts}/runEDTA.sh $genome $cds
${scripts}/runLTRret.sh $genome
scontrol show job $SLURM_JOB_ID
Endofmessage3
sed -i 's/JOBNAME/'$name'_edta/g' ${name}-edta.sub
sbatch ${name}-edta.sub
