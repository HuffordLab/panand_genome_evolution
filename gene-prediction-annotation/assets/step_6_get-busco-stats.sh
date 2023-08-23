#!/bin/bash
dir=$1
busco=$(find "${dir}/06_annotation/busco/mikado-pep-exp-noTE_poales_odb10" -maxdepth 1 -name "short*")
#echo $busco
g=$(grep -h -A 5 "C:" ${busco} | grep -v "C:" | awk '{print $1}' | paste - - - - -)
echo -e "$dir\t$g"
