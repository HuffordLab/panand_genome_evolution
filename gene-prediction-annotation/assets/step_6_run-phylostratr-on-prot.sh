#!/bin/bash
set -x
dir=$1
taxaid=$2
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
prot=$(find ${cwd}/06_annotation/finalFiles -name "*_pep.v1.fasta" )
rundir=$(basename ${prot%.*})
mkdir -p ${cwd}/06_annotation/phylostratr/${rundir}
cd ${cwd}/06_annotation/phylostratr/${rundir}
mkdir -p blastdb uniprot-seqs
sed 's/NCBITAXAID/'${taxaid}'/g' ${scripts}/runDiamondBlasP.sh > uniprot-seqs/runDiamondBlasP.sh
sed 's/NCBITAXAID/'${taxaid}'/g' ${scripts}/runPhylostratR.sh > runPhylostratR.sh
head -n 21 runPhylostratR.sh > temp.sh
chmod +x temp.sh runPhylostratR.sh uniprot-seqs/runDiamondBlasP.sh
ml purge
ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' ${prot} | fold > uniprot-seqs/${taxaid}.faa
cat > ${rundir}-ps.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=amd
#SBATCH --time=2-00:00:00
#SBATCH --account=las
#SBATCH --qos=las
#SBATCH --job-name=JOBNAME
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
cwd=$SLURM_SUBMIT_DIR
SCMD="singularity exec --bind $cwd /work/LAS/mhufford-lab/arnstrm/rstudio_containers/PhylostratR_R.v3.6.2.sif"
ulimit -s unlimited
Endofmessage1
cat >> ${rundir}-ps.sub <<Endofmessage2
fasta=${taxaid}.faa
taxaid=${taxaid}
Endofmessage2
cat >> ${rundir}-ps.sub <<"Endofmessage3"
ml purge
ml singularity
$SCMD ./temp.sh
cd uniprot-seqs
ml purge
ml diamond
for f in *.faa; do diamond makedb --in ${f} --db ${f%.*}; done
for f in *.faa; do ./runDiamondBlasP.sh ${f%.*}; done
rm *.out *.dmnd
mv *.tab $cwd/
ml purge
ml blast-plus
cd ${cwd}/blastdb
for f in ${cwd}/uniprot-seqs/*.faa; do makeblastdb -in $f  -dbtype prot -parse_seqids -out ${f}; done
mv ${cwd}/uniprot-seqs/*.faa.p?? ./
cd ${cwd}
ml purge
ml singularity
$SCMD ./runPhylostratR.sh &> ps.stdout
scontrol show job $SLURM_JOB_ID
Endofmessage3
code=$(echo $rundir |cut -c 1-2)
sed -i 's/JOBNAME/'$name'-ps_'$code'/g' ${rundir}-ps.sub
sbatch ${rundir}-ps.sub
