#!/bin/bash
genome=$1
gff=$2
ml purge
ml bioawk
bioawk -c fastx '{print ">"$name"\n"toupper($seq)}' ${genome} > temp_${genome}
ml bedtools2
if [ ! -f "${genome%.*}.softmasking.done" ]; then
bedtools maskfasta -soft -fi temp_${genome} -bed ${gff} -fo ${genome%.*}.softmasked.fasta || {
echo >&2 "ERROR: softmasking failed"
exit 1
}
touch ${genome%.*}.softmasking.done
rm temp_${genome}
fi
