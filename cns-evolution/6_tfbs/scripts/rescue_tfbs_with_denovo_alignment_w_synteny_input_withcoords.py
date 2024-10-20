import sys
import gzip
import random 
from collections import Counter

# Usage: python rescue_tfbs_with_denovo_alignment.py translationSS/Zm-B73-REFERENCE-NAM-5.0_1000bp_translationSS_JASPARmotifs.bed panand_1000bp_TSS_regions_JASPARmotifs.bed Zm-B73-REFERENCE-NAM-5.0_1000bp_translationSS_JASPARmotifs.classification.wzincarexta.11102023.txt panand_syntenic_orthogroup_genes_1to1_matched_translationSS.tsv

sp1_motifs = sys.argv[1]
sp2_motifs = sys.argv[2]
alntable = sys.argv[3]
gtable = sys.argv[4]

query_species = ["aburma@1","aburma@2","achine@1","achine@2","agerar@1","agerar@2","agerar@3","avirgi@1","blagur@1","blagur@2","blagur@3","ccitra@1","ccitra@2","crefra@1","etrips@1","etrips@2","hcompr@1","hcompr@2","hcompr@3","hcompr@4","hcompr@5","hcompr@6","hconto@1","hconto@2","hconto@3","hconto@4","irugos@1","ppanic@1","rrottb@1","rrottb@2","rtuber@1","sbicol@1","smicro@1","snutan@1","snutan@2","sscopa@1","sscopa@2","tdactn","tdacts","telega@1","telega@2","ttrian@1","udigit@1","udigit@2","udigit@3","vcuspi@1","vcuspi@2","zTIL01","zTIL11","zTIL18","zTIL25","zdiplg","zdiplm","zhuehu","zluxur","znicar"]

# Load reference TFBSs per gene upstream region
refmotifs = {}
coordname = {}
refgene_coords = {}
with open(sp1_motifs,'r') as motifs1:
    for l in motifs1:
        l = l.strip().split('\t')
        refgene = l[3].split('_')[0]
        if refgene not in refgene_coords:
            refgene_coords[refgene] = l[:3] + [l[4]]
        coords = '_'.join(l[5:8])
        motifname = l[8]
        if coords not in coordname:
            coordname[coords] = []
        if 'MA' in motifname and motifname not in coordname[coords]:
            coordname[coords].append(motifname)
        if refgene not in refmotifs:
            refmotifs[refgene] = []
        if motifname != ".":
            motifname_coords = motifname + '@' + '@'.join(l[5:8]+[l[4]])
            refmotifs[refgene].append(motifname_coords)

# Load all motifs in gene upstream regions for query species
qmotifs = {}
with open(sp2_motifs,'r') as motifs2:
    for l in motifs2:
        l = l.strip().split('\t')
        if "Sobic" in l[3]:
            qgene = '.'.join(l[3].split('.')[0:2])
        else:
            qgene = l[3].split('_')[0]
        if qgene not in qmotifs:
            qmotifs[qgene] = []
        if l[8] != ".":
            motif_coords = '@'.join(l[5:8]) + '@' + l[4]
            motifid = l[8].split(':')[0].split(';')[0]
            motifid_coords = motifid + "@" + motif_coords
            qmotifs[qgene].append(motifid_coords)
# Load all the aligned motifs from Cactus
alnmotifs = {}
alngenes = {}
with open(alntable,'r') as motifs3:
    for l in motifs3:
        l = l.strip().split('\t')
        coords = '_'.join(l[0].split('_')[0:3])
        refgene = l[0].split('_')[3]
        subgenome = l[1]
        mismatches = int(l[2])
        gaps = int(l[3])
        qseq = l[5]
        qcoords = l[7]
        len_motifs = int(coords.split('_')[2]) - int(coords.split('_')[1])
        motifs = coordname[coords]
        if coords not in alnmotifs:
            alnmotifs[coords] = []
        if refgene not in alngenes:
            alngenes[refgene] = {}
        # threshold: 2 mismatches (each gap also counts as mismatch)
        if mismatches < len_motifs:
            alnmotifs[coords].append(subgenome)
            if subgenome not in alngenes[refgene]:
                alngenes[refgene][subgenome] = []
            for m in motifs:
                m_counts = m + ";" + qseq + ";" + str(mismatches) + ";" + str(gaps) + ";" + qcoords
                alngenes[refgene][subgenome].append(m_counts)

# Read all syntenic genes and check TFBSs present
#rescue_dict = {}
gene_pairs = {}
with open(gtable,'r') as g:
    for l in g:
        l = l.strip().split('\t')
        #Zm00001eb000170 Ac00001aa044590 achine@2 achine
        subgenome = l[2]
        queryspecies = l[3]
        refgene = l[0]
        querygene = l[1]
        querygene_list = l[4].split(',')
        if refgene not in gene_pairs:
            gene_pairs[refgene] = {}
        if subgenome not in gene_pairs[refgene]:
            gene_pairs[refgene][subgenome] = {'match':querygene,'all_matches':querygene_list}

# Loop over all ref genes and show cactus and denovo matches
for gene,motifs in refmotifs.items():
    for q in query_species:
        syntenic = 'NA'
        denovo_motifs = ['NA']
        alt_syntenic = ['NA']
        # must be paired orthologous collinear gene
        if gene in gene_pairs:
            if q in gene_pairs[gene]:
                syntenic = gene_pairs[gene][q]['match']
                # consider there may be a one-to-many mapping of the ref gene
                alt_syntenic = gene_pairs[gene][q]['all_matches']
                if syntenic in qmotifs:
                    denovo_motifs = qmotifs[syntenic]
        alns = ['NA']
        if gene in alngenes:
            if q in alngenes[gene]:
                alns = alngenes[gene][q]
                
        outline = refgene_coords[gene] + [gene,','.join(motifs),q,','.join(alns),syntenic,','.join(alt_syntenic),','.join(denovo_motifs)]
        print(*outline,sep='\t')


