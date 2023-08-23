#!/bin/bash
dir=$1

# job_header
read -r -d '' VAR <<EOF1
#!/bin/bash\n
#SBATCH --nodes=1\n
#SBATCH --exclusive\n
#SBATCH --mem=350GB\n
#SBATCH --time=4-00:00:00\n
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
#SBATCH --partition=scavenger\n

# some variables
cd ${dir}
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
genome="${cwd}/00_raw-data/genome/${dir}*_p.fasta"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
bam=$(find $(pwd)/03_STAR -name "merged*.bam")
# directory for assemblies
mkdir -p ${cwd}/04_transcript-assembly
cd ${cwd}/04_transcript-assembly
subfolders=(strawberry class2 stringtie trinity cufflinks portcullis braker)
for subfolder in "${subfolders[@]}"; do
    mkdir -p ${cwd}/04_transcript-assembly/${subfolder}
    cd ${cwd}/04_transcript-assembly/${subfolder}
    ln -s ${bam}
    ln -s ${genome}
done
#class2 job
cd ${cwd}/04_transcript-assembly/class2;
echo -e ${VAR} > calss2.sub
echo "$scripts/runPSIClass.sh $(basename $bam)" >> calss2.sub
sed -i 's/JOBNAME/'${name}'-cl/g' calss2.sub
# stringtie job
cd ${cwd}/04_transcript-assembly/stringtie;
echo -e ${VAR} > stringtie.sub
echo "$scripts/runStringtie.sh $(basename $bam)" >> stringtie.sub
sed -i 's/JOBNAME/'${name}'-st/g' stringtie.sub
#trinity job
cd ${cwd}/04_transcript-assembly/trinity;
echo -e ${VAR} > trinity.sub
echo "$scripts/runTrinity.sh $(basename $bam)" >> trinity.sub
sed -i 's/JOBNAME/'${name}'-tr/g' trinity.sub
#cufflinks job
cd ${cwd}/04_transcript-assembly/cufflinks;
echo -e ${VAR} > cufflinks.sub
echo "$scripts/runCufflinks.sh $(basename $bam)" >> cufflinks.sub
sed -i 's/JOBNAME/'${name}'-cf/g' cufflinks.sub
#portcullis job
cd ${cwd}/04_transcript-assembly/portcullis;
echo -e ${VAR} > portcullis.sub
echo "$scripts/runPortcullis.sh $(basename $bam) $(basename ${genome})" >> portcullis.sub
sed -i 's/JOBNAME/'${name}'-po/g' portcullis.sub
#braker job
cd ${cwd}/04_transcript-assembly/braker;
echo -e ${VAR} > braker.sub
echo "$scripts/runBRAKER-RNAseq.sh $(basename ${genome}) $(basename $bam)" >> braker.sub
sed -i 's/JOBNAME/'${name}'-br/g' braker.sub
#strawberry job
cd ${cwd}/04_transcript-assembly/strawberry;
echo -e ${VAR} > strawberry.sub
echo "$scripts/runStrawberry.sh $(basename $bam)" >> strawberry.sub
sed -i 's/JOBNAME/'${name}'-sw/g' strawberry.sub
#submit jobs
for subfolder in "${subfolders[@]}"; do
    cd ${cwd}/04_transcript-assembly/$subfolder
    sed -i 's/^[[:space:]]*//g' *.sub
    sbatch *.sub
    echo $(pwd)
done
