#!/bin/bash
dir=$1
# functions
gtf () {
  genes=$(grep -v "^#" $1 | awk '$3=="gene"' | wc -l)
  trans=$(grep -v "^#" $1 | awk '$3=="transcript"' | wc -l)
  echo -e "${dir}\t$genes\t$trans"
}
gff() {
  genes=$(grep -v "^#" $1 | awk '$3=="gene"' | wc -l)
  trans=$(grep -v "^#" $1 | awk '$3=="mRNA"' | wc -l)
  echo -e "${dir}\t$genes\t$trans"
}
fcount() {
  trans=$(grep -c ">" $1)
  genes=$(grep ">" $1 |rev |cut -f 2- -d "_" |sort |uniq |wc -l)
  echo -e "${dir}\t$genes\t$trans"
}
# some variables
cd ${dir}
cwd=$(pwd)
scripts="/work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.v2.triffid/"
genome="${cwd}/00_raw-data/genome/${dir}_v1_p.fasta"
name="$(echo $dir | cut -c 1)$(echo $dir | cut -f 2 -d "-" | cut -c 1-3)"
bam=$(find $(pwd)/03_STAR -name "merged*.bam")
# directory for assemblies
cd ${cwd}/04_transcript-assembly
subfolders=(strawberry class2 stringtie trinity cufflinks braker mikado_v2.3)
braker=$(find "$(pwd)/braker/braker" -maxdepth 1 -name "braker.gtf")
stringtie=$(find "$(pwd)/stringtie" -name "merged*stringtie.gtf")
class2=$(find "$(pwd)/class2" -name "*_psiclass.gtf_vote.gtf")
cufflinks=$(find "$(pwd)/cufflinks" -name "transcripts.gtf")
trinity=$(find "$(pwd)/trinity/trinity_out_dir"  -maxdepth 1 -name "Trinity-GG.fasta")
strawberry=$(find "$(pwd)/strawberry" -name "*_strawberry.gtf")
mikado=$(find "$(pwd)/mikado_v2.3" -name "mikado.genes.gff3")
mikadoFilt=$(find "$(pwd)/mikado_v2.3" -name "mikado.genes-filtered.gff3")
br=$(gtf ${braker})
cl=$(gtf ${class2})
st=$(gtf ${stringtie})
sb=$(gtf ${strawberry})
tr=$(fcount ${trinity})
cf=$(gtf ${cufflinks})
mk=$(gff ${mikado})
mf=$(gff ${mikadoFilt})
# print counts:
echo -e "method\tgenome\tgenes\ttranscripts"
echo -e "braker\t$br"
echo -e "cufflinks\t$cl"
echo -e "stringtie\t$st"
echo -e "strawberry\t$sb"
echo -e "trinity\t$tr"
echo -e "class2\t$cf"
echo -e "mikado\t$mk"
echo -e "mikado_filtered\t$mf"
