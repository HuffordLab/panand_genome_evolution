#!/bin/bash
# shell script to run TRaCE pipeline
# author: Sheng-Kai Hsu
# date created: 2022.08.10
# date last edited: 2022.09.15

cd /workdir/sh2246/p_panAndOGASR

## step0: manually number the files based on species from 01 to 36

## step1: run interProScan

parallel -j 3 'gffread -g data/assembly/final_8_19_22/{}_*.fasta -y output/aa/{}.aa.fa data/annotation/{}_*.gff3' ::: {34..36}


mkdir output/pfam
parallel -j 3 'interproscan-5.44-79.0/interproscan.sh -i output/aa/{}.aa.fa -f tsv -appl Pfam -d output/pfam -cpu 1' ::: {34..36}

#interproscan-5.44-79.0/interproscan.sh -i output/aa/Zmays.aa.modified.fa -f tsv -appl Pfam -d output/pfam -cpu 40


## step 2: stringtie assembly
# for bam in *.bam; do
# 	bash ../src/HuffordLab-NAM-genomes/gene-prediction/scripts-evidence/runStringtie.sh $bam
# done
parallel -j 3 'mkdir output/{}_out_stringtie' ::: {34..36}

parallel -j 3 'for bam in `ls data/RNASeq/{}*/*.bam`; do bash src/HuffordLab-NAM-genomes/gene-prediction/scripts-evidence/runStringtie.sh $bam output/{}_out_stringtie ;done' ::: {34..36}


## step 4: run TRaCE
#cd ..
# perl data/TRaCE/TRaCE.pl data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 output/pfam/Zmays.aa.modified.fa.tsv 0.5 0.5 0.5 0.5 output/*out.bam_stringtie.gtf > output/trace.out

# although there're default values for the four threshold parameters, you still need to provide it. Otherwise, there will be problem reading the input gtf

mkdir output/trace_output/
parallel -j 3 'perl data/TRaCE/TRaCE_SKH.pl data/annotation/{}_*.gff3 output/pfam/{}.aa.fa.tsv 0.5 0.5 0.5 0.5 output/{}_out_stringtie/*_stringtie.gtf > output/trace_output/{}_trace.out' ::: {34..36}
