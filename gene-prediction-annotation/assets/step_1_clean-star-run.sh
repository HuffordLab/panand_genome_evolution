#!/bin/bash
mkdir -p 00_raw-data/fastq-files 01_STAR/job-files 01_STAR/star-misc-files 01_STAR/bam-files/ 00_raw-data/genome
mv *.fq.gz 00_raw-data/fastq-files
mv *_0.sub *_0.o* *_0.e* *.cmds SJ.all bam.fofn *.out *.tab *_STARtmp *_STARgenome 01_STAR/star-misc-files/
index=$(basename $(pwd))
mv $index 01_STAR/STAR_$index
mv *sortedByCoord.out.bam 01_STAR/bam-files/
mv *.bam ./01_STAR/
mv *.fasta 00_raw-data/genome/
