#!/bin/bash
for f in *.fasta; do
g=$(echo $f |sed 's/_v..fasta//g');
mkdir -p $g;
mv $f ./$g;
done
#mkdir -p Elionorus-tripsacoides
#mkdir -p Hemarthria-compressa
#mkdir -p Themeda-triandra
