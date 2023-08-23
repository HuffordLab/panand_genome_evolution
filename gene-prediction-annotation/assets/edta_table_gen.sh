#!/bin/bash
dir=$1
stats=$(find $(pwd)/${dir}/07_edta-run -maxdepth 1 -type f -name "*mod.EDTA.TEanno.sum")
awk '/Repeat Classes/,/Repeat Stats/' ${stats} |awk 'NF>0 {$1=$1;print}' |sed 's/ /\t/g' |grep -v "^=" |grep -v "^-" |grep -v "Repeat" |sed 's/total\tinterspersed/total_interspersed/g' | awk '/^Class/,/Total/' > ${dir}_edta_table.tsv
