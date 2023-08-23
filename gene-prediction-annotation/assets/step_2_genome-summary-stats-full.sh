#!/bin/bash
dir=$1
cd ${dir}
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2"
genome=$(find ${cwd}/00_raw-data/genome -name "${dir}_v[012].fasta")
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"

cd $cwd/00_raw-data/genome
cat > script1.sh <<Endofmessage1
#!/bin/bash
ml purge
ml bioawk
bioawk -c fastx '{print \$name}' $genome > .all-names.txt
grep -v "^alt" .all-names.txt > primary-scafs.txt
grep "^alt" .all-names.txt > alt-scafs.txt
seqtk subseq $genome primary-scafs.txt > ${genome%.*}_p.fasta
seqtk subseq $genome alt-scafs.txt > ${genome%.*}_a.fasta
rm .all-names.txt
Endofmessage1
chmod +x script1.sh

primary=$(realpath ${genome%.*}_p.fasta)
alt=$(realpath ${genome%.*}_a.fasta)
cd $cwd
mkdir -p 02_genome-stats
cd 02_genome-stats
size=$(grep -F "${dir}" /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//sizes-genome.tsv |awk '{print $2*1000000}')
mkdir summary-stats
cd summary-stats
cat > script2.sh <<Endofmessage2
#!/bin/bash
new_Assemblathon.pl -csv -genome_size $size $primary > primary.stats
new_Assemblathon.pl -csv -genome_size $size $alt > secondary.stats
new_Assemblathon.pl -csv -genome_size $size $genome > combined.stats
Endofmessage2
chmod +x script2.sh

mkdir -p $cwd/02_genome-stats/busco-stats
cd $cwd/02_genome-stats/busco-stats
ln -s $primary
ln -s $alt
ln -s $genome
ml purge
cat > script3.sh <<Endofmessage3
#!/bin/bash
${scripts}/runBUSCO.sh $(basename $genome)
${scripts}/runBUSCO.sh $(basename $primary)
${scripts}/runBUSCO.sh $(basename $alt)
Endofmessage3
chmod +x script3.sh
cd $cwd/02_genome-stats
cat > ${name}-stats.sub <<"Endofmessage4"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=4-00:00:00
#SBATCH --account=las
#SBATCH --qos=las
#SBATCH --partition=amd
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
Endofmessage4

cat >> ${name}-stats.sub <<Endofmessage5
cd $cwd/00_raw-data/genome
./script1.sh
cd $cwd/02_genome-stats/summary-stats
./script2.sh
cd $cwd/02_genome-stats/busco-stats
./script3.sh
Endofmessage5
sed -i 's/JOBNAME/'$name'_stats/g' ${name}-stats.sub
#sbatch ${name}-stats.sub

