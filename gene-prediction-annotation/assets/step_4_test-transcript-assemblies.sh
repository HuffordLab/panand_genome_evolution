#!/bin/bash
dir=$1
cd ${dir}
cwd=$(pwd)
# directory for assemblies
cd ${cwd}/04_transcript-assembly
subfolders=(strawberry class2 stringtie trinity cufflinks braker portcullis)
braker=$(find "$(pwd)/braker/braker" -maxdepth 1 -name "braker.gtf")
stringtie=$(find "$(pwd)/stringtie" -name "merged*stringtie.gtf")
class2=$(find "$(pwd)/class2" -name "*_psiclass.gtf_vote.gtf")
cufflinks=$(find "$(pwd)/cufflinks" -name "transcripts.gtf")
trinity=$(find "$(pwd)/trinity/trinity_out_dir"  -maxdepth 1 -name "Trinity-GG.fasta")
strawberry=$(find "$(pwd)/strawberry" -name "*_strawberry.gtf")
portcullis=$(find "$(pwd)/portcullis/portcullis_out/3-filt" -name "portcullis_filtered.pass.junctions.bed")
for file in $braker $class2 $cufflinks $trinity $strawberry $portcullis; do
if [[ -s $file  ]]; then
echo "$(basename $file) exists"
else
echo "$(basename $file) not found"
fi
done


