#!/bin/bash
# should filter in 06_annotation/filterTE folder
dir=$1
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
# gff files
primary="${cwd}/00_raw-data/genome/${dir}_*_p.fasta"
braker=$(find ${cwd}/04_transcript-assembly/braker/braker -maxdepth 1 -name "braker.gtf")
mikado=$(find ${cwd}/04_transcript-assembly/mikado_v2.3 -name "mikado.genes-filtered.gff3")
bind=$(find ${cwd}/04_transcript-assembly/gemoma/gemoma_output_v1 -name "final_annotation.gff")
# tsv files

teBraker=$(find ${cwd}/06_annotation/tesorter/braker -maxdepth 1 -name "*rexdb-plant.cls.tsv")
teMikado=$(find ${cwd}/06_annotation/tesorter/mikado -maxdepth 1 -name "*rexdb-plant.cls.tsv")
teBind=$(find ${cwd}/06_annotation/tesorter/bind -maxdepth 1 -name "*rexdb-plant.cls.tsv")



# job_header
read -r -d '' VAR <<EOF1
#!/bin/bash\n
#SBATCH --nodes=1\n
#SBATCH -n 2\n
#SBATCH --partition=amd\n
#SBATCH --time=1:00:00\n
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
mkdir -p 06_annotation/filterTE/{braker,mikado,bind}

# for braker
cd 06_annotation/filterTE/braker
ln -s ${braker}
ln -s ${teBraker}
ln -s ${primary}
echo -e ${VAR} > braker.sub
echo "$scripts/runFilterTE.sh $(basename $braker) $(basename ${teBraker}) $(basename ${primary}) ${name}" >> braker.sub
sed -i 's/JOBNAME/'${name}'-BrTE/g' braker.sub
echo "scontrol show job \$SLURM_JOB_ID" >> braker.sub
sed -i 's/^[[:space:]]*//g' braker.sub
sbatch braker.sub

# for mikado
cd ${cwd}/06_annotation/filterTE/mikado
ln -s ${mikado}
ln -s ${teMikado}
ln -s ${primary}
echo -e ${VAR} > mikado.sub
echo "$scripts/runFilterTE.sh $(basename $mikado) $(basename ${teMikado}) $(basename ${primary}) ${name}" >> mikado.sub
echo "scontrol show job \$SLURM_JOB_ID" >> mikado.sub
sed -i 's/JOBNAME/'${name}'-MiTE/g' mikado.sub
sed -i 's/^[[:space:]]*//g' mikado.sub
sbatch mikado.sub

# for bind
cd ${cwd}/06_annotation/filterTE/bind
ln -s ${bind}
ln -s ${teBind}
ln -s ${primary}
echo -e ${VAR} > bind.sub
echo "$scripts/runFilterTE.sh $(basename $bind) $(basename ${teBind}) $(basename ${primary}) ${name}" >> bind.sub
echo "scontrol show job \$SLURM_JOB_ID" >> bind.sub
sed -i 's/JOBNAME/'${name}'-IdTE/g' bind.sub
sed -i 's/^[[:space:]]*//g' bind.sub
sbatch bind.sub
