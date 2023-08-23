#!/bin/bash
dir=$1
stats=$(find $(pwd)/${dir}/06_annotation/finalFiles -type f -name "*.stats")
awk '/Compute/,/Re-compute/' $stats |sed 's/ /_/g' |sed 's/__/\t/g' |sed 's/\t_/\t/g' |tr -s "\t" |awk 'NF>1' > ${dir}_curated_table.tsv
