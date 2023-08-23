#!/bin/bash
if [ "$#" -ne 1 ] ; then
echo "please provide:"
echo -e "\t\t(1) genome fasta file (genome.fasta)"
echo "";
echo "step_8-extract.sh  <genome.fasta>";
echo "";
exit 0;
fi
genome=$1
ml purge
ml cufflinks bioawk
gffread mikado.genes.gff3 -g ${genome} -t mRNA -x mikado.transcripts.fasta -y temp_mikado.proteins.fasta
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' temp_mikado.proteins.fasta | fold > mikado.proteins.fasta
rm temp_mikado.proteins.fasta

