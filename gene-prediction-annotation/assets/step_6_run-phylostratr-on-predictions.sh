#!/bin/bash
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -p|--prediction)
      prediction=both
      shift # past argument
      shift # past value
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done
if [ "${#POSITIONAL[@]}" -lt 2 ] ; then
echo "needs a genus-species name as first argument and its taxa-id as second"
echo "./step_6_run-phylostratr-on-predictions.sh <genus-species> <taxa-id> [<prediction-type: braker, mikado, both>]"
echo "prediction-type: default is both"
exit 0
fi
dir=${POSITIONAL[0]}
taxaid=${POSITIONAL[1]}
cd $dir
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
braker=$(find $(pwd)/04_transcript-assembly/braker/braker -maxdepth 1 -name "augustus.hints.aa")
mikado=$(find $(pwd)/04_transcript-assembly/mikado_v2.3 -name "mikado.proteins.fasta")
if [ "${prediction}" == "both" ]; then
    echo "both"
    option=(barker mikado)
elif [ "${prediction}" == "braker" ]; then
    echo "braker"
    option=braker
elif [ "${prediction}" == "mikado" ]; then
    echo "mikado"
    option=mikado
else
    echo "default is both"
    option=(braker mikado)
fi
mkdir -p ${cwd}/06_annotation/phylostratr
cd ${cwd}/06_annotation/phylostratr
ln -s $braker braker.faa
ln -s $mikado mikado.faa
for f in ${option[@]}; do
mkdir -p ${cwd}/06_annotation/phylostratr/${f};
cd ${cwd}/06_annotation/phylostratr/${f};
mkdir -p blastdb uniprot-seqs
sed 's/NCBITAXAID/'${taxaid}'/g' ${scripts}/runDiamondBlasP.sh > uniprot-seqs/runDiamondBlasP.sh
sed 's/NCBITAXAID/'${taxaid}'/g' ${scripts}/runPhylostratR.sh > runPhylostratR.sh
head -n 20 runPhylostratR.sh > temp.sh
chmod +x temp.sh runPhylostratR.sh uniprot-seqs/runDiamondBlasP.sh
ml purge
ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' ${cwd}/06_annotation/phylostratr/${f}.faa | fold > uniprot-seqs/${taxaid}.faa
cat > ${f}-ps.sub <<"Endofmessage1"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=amd
#SBATCH --time=12:00:00
#SBATCH --account=triffid
#SBATCH --qos=triffid
#SBATCH --job-name=JOBNAME
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
cwd=$SLURM_SUBMIT_DIR
ulimit -s unlimited
Endofmessage1
cat >> ${f}-ps.sub <<Endofmessage2
fasta=${taxaid}.faa
taxaid=${taxaid}
Endofmessage2
cat >> ${f}-ps.sub <<"Endofmessage3"
ml purge
ml r-devtools blast-plus
./temp.sh
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
ml r-devtools blast-plus
./runPhylostratR.sh &> ps.stdout
scontrol show job $SLURM_JOB_ID
Endofmessage3
code=$(echo $f |cut -c 1-2)
sed -i 's/JOBNAME/'$name'-ps_'$code'/g' ${f}-ps.sub
sbatch ${f}-ps.sub
done
