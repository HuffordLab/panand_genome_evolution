#!/bin/bash
dir="03_STAR"
if [[ -d "$dir" ]]
then
    echo "dir exists.. skipping"
    exit 0;
fi


mkdir -p 00_raw-data/fastq-files 03_STAR/job-files 03_STAR/star-misc-files 03_STAR/bam-files/
mv *.fq.gz 00_raw-data/fastq-files
mv *_0.sub *_0.o* *_0.e* 03_STAR/job-files/
mv *.cmds SJ.all bam.fofn *.out *.tab *_STARtmp *_STARgenome 03_STAR/star-misc-files/
index=$(basename $(pwd))
mv $index 01_STAR/STAR_${index}_p
mv *sortedByCoord.out.bam 03_STAR/bam-files/
mv *.bam ./03_STAR/
