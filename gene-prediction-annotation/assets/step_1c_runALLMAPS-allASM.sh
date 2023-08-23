#!/bin/bash
dir=$1
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"

hap1=$(find 0_data/assembly -name "*.bp.hap1.p_ctg.fasta")
hap2=$(find 0_data/assembly -name "*.bp.hap2.p_ctg.fasta")
primary=$(find 0_data/assembly -name "*.bp.p_ctg.fasta")

mkdir -p ${cwd}/1c_allmaps/{hap1,hap2,primary}
cp ${hap1} ${cwd}/1c_allmaps/hap1/hap1.ctg.fasta
cp ${hap2} ${cwd}/1c_allmaps/hap2/hap2.ctg.fasta
cp ${primary} ${cwd}/1c_allmaps/primary/primary.ctg.fasta

rm .temp &> /dev/null
cat > .temp <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=amd
#SBATCH --time=96:00:00
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
Endofmessage1
sed 's/JOBNAME/pg.hap1/g' .temp > ${cwd}/1c_allmaps/hap1/pg.sub
sed 's/JOBNAME/pg.hap2/g' .temp > ${cwd}/1c_allmaps/hap2/pg.sub
sed 's/JOBNAME/pg.prim/g' .temp > ${cwd}/1c_allmaps/primary/pg.sub
rm .temp &> /dev/null
cat >> ${cwd}/1c_allmaps/hap1/pg.sub <<"Endofmessage2a"
genome=hap1.ctg.fasta
/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2/runPanGenomeALLMAPS.sh ${genome} ${genome%.*}
for f in B73 TIL11 TIL25 TIL18 TIL01 Sorghum; do /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/make-scfs-per-chr-for-dotplots.sh pangenome.agp ${genome} ${f}; done
rm chr{1..10}_scafs.txt
rm chr{1..10}_scafs.fasta
for f in B73 TIL11 TIL25 TIL18 TIL01 Sorghum; do /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/make-chr-to-chr-dotplots.sh pangenome.chr.fasta $f; done
Endofmessage2a

cat >> ${cwd}/1c_allmaps/hap2/pg.sub <<"Endofmessage2b"
genome=hap2.ctg.fasta
/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2/runPanGenomeALLMAPS.sh ${genome} ${genome%.*}
for f in B73 TIL11 TIL25 TIL18 TIL01 Sorghum; do /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/make-scfs-per-chr-for-dotplots.sh pangenome.agp ${genome} ${f}; done
rm chr{1..10}_scafs.txt
rm chr{1..10}_scafs.fasta
for f in B73 TIL11 TIL25 TIL18 TIL01 Sorghum; do /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/make-chr-to-chr-dotplots.sh pangenome.chr.fasta $f; done
Endofmessage2b

cat >> ${cwd}/1c_allmaps/primary/pg.sub <<"Endofmessage2c"
genome=primary.ctg.fasta
/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2/runPanGenomeALLMAPS.sh ${genome} ${genome%.*}
for f in B73 TIL11 TIL25 TIL18 TIL01 Sorghum; do /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/make-scfs-per-chr-for-dotplots.sh pangenome.agp ${genome} ${f}; done
rm chr{1..10}_scafs.txt
rm chr{1..10}_scafs.fasta
for f in B73 TIL11 TIL25 TIL18 TIL01 Sorghum; do /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/make-chr-to-chr-dotplots.sh pangenome.chr.fasta $f; done
Endofmessage2c
cd ${cwd}/1c_allmaps/primary
sbatch pg.sub
cd ${cwd}/1c_allmaps/hap2
sbatch pg.sub
cd ${cwd}/1c_allmaps/hap1
sbatch pg.sub

