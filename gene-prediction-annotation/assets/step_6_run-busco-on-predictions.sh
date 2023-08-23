#!/bin/bash
dir=$1
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
bind=$(find ${cwd}/06_annotation/finalFiles -maxdepth 1 -name "*pep.v*.fasta")
ml purge
ml bioawk
mkdir -p ${cwd}/06_annotation/busco
cd ${cwd}/06_annotation/busco
ln -s ${bind}
echo "preparing job for ${dir}"
cat > ${name}-busco.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH -n 16
#SBATCH --partition=amd
#SBATCH --time=12:00:00
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
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate busco5
Endofmessage1
cat >> ${name}-busco.sub <<"Endofmessage2"
for prot in *.fasta; do
if [ ! -f "${prot%.*}_liliopsida_odb10.done" ]; then
busco -i $prot -c ${cpus} -o ${prot%.*}_liliopsida_odb10 -m prot -l liliopsida_odb10 -r || {
echo >&2 "ERROR: BUSCO on ${prot%.*} (liliopsida_odb10) failed. Exiting.."
exit 1
}
touch "${prot%.*}_liliopsida_odb10.done"
fi
if [ ! -f "${prot%.*}_poales_odb10.done" ]; then
busco -i $prot -c ${cpus} -o ${prot%.*}_poales_odb10 -m prot -l poales_odb10  -r || {
echo >&2 "ERROR: BUSCO on ${prot%.*} (poales_odb10) failed. Exiting.."
exit 1
}
touch "${prot%.*}_poales_odb10.done"
fi
done
scontrol show job $SLURM_JOB_ID
Endofmessage2
sed -i 's/JOBNAME/'$name'_busco/g' ${name}-busco.sub
sbatch ${name}-busco.sub
echo "Done.."
