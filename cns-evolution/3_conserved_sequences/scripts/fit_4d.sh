#$ -S /usr/bash
#$ -cwd
#$ -pe threads 1
#$ -N fit4d
#$ -o fit4d.out
#$ -j y
#$ -l m_mem_free=4G
## Run

source ~/.bashrc
conda activate batgene
phyloFit --tree astral_panand_tree.newick chr10_4d.maf --out-root panand_10_4d_single_nozea_mod
