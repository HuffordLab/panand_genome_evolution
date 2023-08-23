#!/bin/bash
dir=$1
cd $dir
cwd=$(pwd)
mkdir -p 03_STAR/qcfiles
cd 03_STAR/qcfiles

for file in $(find ${cwd}/03_STAR/star-misc-files -name "*2Log.final.out"); do
ln -s $file;
done
multiqc .
mv multiqc_report.html ${dir}_mulitqc.html
