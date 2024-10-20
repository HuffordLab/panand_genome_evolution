import sys
from collections import defaultdict
import itertools
from Bio import AlignIO
# Usage
bases = ["a","c","t","g","A","C","T","G"]

for multiple_alignment in AlignIO.parse(sys.argv[1], "maf"):
    #aln_dict = {"aburma":0,"achine":0,"agerar":0,"avirgi":0,"blagur":0,"ccitra":0,"crefra":0,"etrips":0,"hconto":0,"hcompr":0,"irugos":0,"ppanic":0,"rrottb":0,"rtuber":0,"sbicol":0,"smicro":0,"snutan":0,"sscopa":0,"tdacts":0,"tdactn":0,"telega":0,"ttrian":0,"tzopol":0,"udigit":0,"vcuspi":0,"zdiplg":0,"zdiplm":0,"zhuehu":0,"zluxur":0,"znicar":0,"zTIL01":0,"zTIL11":0,"zTIL18":0,"zTIL25":0,"zB73v5":0}
    # MAF format requires all sequences to be same length
    ref_record = multiple_alignment[0]
    aln_len = len(ref_record.seq)
    chrom = "".join(str(ref_record.id).split(".")[1:])
    start = int(ref_record.annotations["start"])
    size = int(ref_record.annotations["size"])
    end = start + size
    tripcnt = 0
    nontripcnt = 0
    tripsacumcnt = 0
    present = []
    seqdict = {}
    for record in multiple_alignment:
        sp = str(record.id).split(".")[0]
        seq = str(record.seq)
        percent_gap = seq.count("-") / len(seq)
        # do not count highly gapped alignments
        if sp not in present and percent_gap <=0.5:
            present.append(sp)
            seqdict[sp] = str(record.seq)
        #refseq = str(record.seq)
        #gapcnt = refseq.count("-")
        #ncnt = refseq.count("N") + refseq.count("n")
        #alncnt = len(refseq) - gapcnt - ncnt
        #aln_dict[sp] = alncnt
    outseq = []
    species_present = set()
    for p in present:
        outseq.append(seqdict[p])
        species = p.split('@')[0]
        species_present.add(species)
    
    outline = [chrom,start,end] + [len(present),len(species_present)] + [','.join(present)] + [','.join(outseq)]
    print(*outline,sep="\t")






