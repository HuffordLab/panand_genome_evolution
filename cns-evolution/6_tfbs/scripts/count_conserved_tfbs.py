import sys
from Bio import SeqIO
from difflib import SequenceMatcher

def get_number_of_character_differences(string1, string2):
    mismatches = 0
    non_nuc_chars = 0
    for i,s in enumerate(string1):
        if s != string2[i] or s == "-":
            mismatches += 1
            if string2[i] not in ["A","C","T","G"]:
                non_nuc_chars += 1
    return mismatches,non_nuc_chars

with open(sys.argv[1],'r') as infile:
    for l in infile:
        l = l.strip().split('\t')
        seqdict = {}
        coord_dict = {}
        for idx,subgenome in enumerate(l[11].split(',')):
            subgenome_seq = l[12].split(',')[idx]
            seqdict[subgenome] = subgenome_seq
            coord_dict[subgenome] = l[13].split(',')[idx]
        for k,v in seqdict.items():
            if k != 'consensus':
                mismatches,gaps = get_number_of_character_differences(seqdict['zB73v5'].upper(),v.upper())
                outline = l[:10] + [k,seqdict['zB73v5'].upper(),v.upper(),mismatches,gaps,coord_dict[k]]
                print(*outline,sep='\t')
