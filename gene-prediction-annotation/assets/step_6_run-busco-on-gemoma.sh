#!/bin/bash
dir=$1
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
gemoma_v1=$(find ${cwd}/04_transcript-assembly/gemoma/gemoma_output_v1 -name "gemoma.v1_proteins.fasta")
gemoma_v2=$(find ${cwd}/04_transcript-assembly/gemoma/gemoma_output_v2 -name "gemoma.v2_proteins.fasta")

mkdir -p ${cwd}/06_annotation/busco
cd ${cwd}/06_annotation/busco
ln -s ${gemoma_v1} gemoma_pep_v1.fasta
ln -s ${gemoma_v2} gemoma_pep_v2.fasta

cat > ${name}-busco2.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --partition=scavenger
#SBATCH --time=12:00:00
#SBATCH --account=triffid
#SBATCH --qos=triffid
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
cat >> ${name}-busco2.sub <<"Endofmessage2"
for prot in gemoma_pep_v1.fasta gemoma_pep_v2.fasta; do
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
sed -i 's/JOBNAME/'$name'_busco2/g' ${name}-busco2.sub
sbatch ${name}-busco2.sub
