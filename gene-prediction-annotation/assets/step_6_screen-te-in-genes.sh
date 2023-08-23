#!/bin/bash
dir=$1
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
braker=$(find ${cwd}/04_transcript-assembly/braker/braker -maxdepth 1 -name "augustus.hints.codingseq")
mikado=$(find ${cwd}/04_transcript-assembly/mikado_v2.3 -name "mikado.transcripts.fasta")
bind=$(find ${cwd}/04_transcript-assembly/gemoma/gemoma_output_v1 -name "predicted_cds.fasta")
# job_header
read -r -d '' VAR <<EOF1
#!/bin/bash\n
#SBATCH --nodes=1\n
#SBATCH -n 16\n
#SBATCH --partition=amd\n
#SBATCH --time=12:00:00\n
#SBATCH --account=las\n
#SBATCH --qos=las\n
#SBATCH --job-name=JOBNAME\n
#SBATCH --output=nova-%x.%j.out\n
#SBATCH --error=nova-%x.%j.err\n
#SBATCH --mail-user=arnstrm@gmail.com\n
#SBATCH --mail-type=begin\n
#SBATCH --mail-type=end\n
cd \$SLURM_SUBMIT_DIR\n
cwd=\$SLURM_SUBMIT_DIR\n
ulimit -s unlimited\n
ml purge\n
EOF1
mkdir -p 06_annotation/tesorter/{braker,mikado,bind}
cd 06_annotation/tesorter/braker
ln -s ${braker}
echo -e ${VAR} > braker.sub
echo "$scripts/runTEsorter.sh $(basename $braker)" >> braker.sub
sed -i 's/JOBNAME/'${name}'-BrTE/g' braker.sub
echo "scontrol show job \$SLURM_JOB_ID" >> braker.sub
sed -i 's/^[[:space:]]*//g' braker.sub
##sbatch braker.sub
cd ${cwd}/06_annotation/tesorter/mikado
ln -s ${mikado}
echo -e ${VAR} > mikado.sub
echo "$scripts/runTEsorter.sh $(basename $mikado)" >> mikado.sub
echo "scontrol show job \$SLURM_JOB_ID" >> mikado.sub
sed -i 's/JOBNAME/'${name}'-MiTE/g' mikado.sub
sed -i 's/^[[:space:]]*//g' mikado.sub
##sbatch mikado.sub
cd ${cwd}/06_annotation/tesorter/bind
ln -s ${bind}
echo -e ${VAR} > bind.sub
echo "$scripts/runTEsorter.sh $(basename $bind)" >> bind.sub
echo "scontrol show job \$SLURM_JOB_ID" >> bind.sub
sed -i 's/^[[:space:]]*//g' bind.sub
##sbatch bind.sub
