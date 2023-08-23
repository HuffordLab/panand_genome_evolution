#!/bin/bash
if [ "$#" -ne 2 ] ; then
echo "please provide:"
echo -e "\t\t(1) list file"
echo -e "\t\t(2) genome file"
echo "";
echo "step_2-prepare.sh <list.txt> <genome.fasta>" ;
echo "";
exit 0;
fi
list=$1
genome=$2
junctions=junctions.bed
# run configure
mikado configure --list ${list} --reference ${genome} --mode permissive --scoring plant.yaml --copy-scoring plant.yaml --junctions ${junctions}  --toml > config.toml
cp /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid//plants.yaml ./
config=$(realpath plants.yaml)
sed -i 's+plant.yaml+'$config'+g' config.toml
