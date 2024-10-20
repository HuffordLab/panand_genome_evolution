import sys
from collections import defaultdict
import itertools
from Bio import AlignIO
# Usage
bases = ["a","c","t","g","A","C","T","G"]

for multiple_alignment in AlignIO.parse(sys.argv[1], "maf"):
    ref_record = multiple_alignment[0]
    aln_len = len(ref_record.seq)
    chrom = "".join(str(ref_record.id).split(".")[1:])
    start = int(ref_record.annotations["start"])
    size = int(ref_record.annotations["size"])
    end = start + size
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
    outseq = []
    species_present = set()
    for p in present:
        outseq.append(seqdict[p])
        species = p.split('@')[0]
        species_present.add(species)
    outline = [chrom,start,end] + [len(present),len(species_present)] + [','.join(present)] + [','.join(outseq)]
    print(*outline,sep="\t")






