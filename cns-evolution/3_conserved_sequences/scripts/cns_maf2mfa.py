import sys
import itertools
from Bio import AlignIO

def maf2fa(inmaf,species):
    maf = AlignIO.parse(inmaf, "maf")
    # get list of species
    alndict = {}
    species = species.split(",")
    for line in species:
        #line = line.strip()
        #line = line.split(",")
        if line not in alndict:
            #spdict[line] = []
            alndict[line] = []
    num_sp = len(alndict)
    #print(spdict,csdict)
    for multiple_alignment in maf:
        spcheck = []
        aln_list = []
        # MAF format requires all sequences to be same length
        ref_record = multiple_alignment[0]
        gap_indices = []
        for i,n in enumerate(ref_record.seq):
            if n == "-":
                gap_indices.append(i)

        aln_len = len(ref_record.seq) - len(gap_indices)
        chrom = "".join(str(ref_record.id).split(".")[1:])
        start = int(ref_record.annotations["start"])
        size = int(ref_record.annotations["size"])
        srcsize = int(ref_record.annotations["srcSize"]) 
        end = start + size
        for record in multiple_alignment:
            sp = str(record.id).split(".")[0]
            #assume first record is ref
            chrom = str(record.id)
            species = chrom.split(".")[0]
            #print(species,alndict)
            spcheck.append(species)
            sequence = []
            for i,n in enumerate(record.seq.upper()):
                if i not in gap_indices:
                    sequence.append(n)
            sequence = ''.join(sequence)
            alndict[species].append(sequence)
            #alndict[species].append(str(record.seq.upper()))
            #else:
                # Replace duplicated entries with gaps
            #    alndict[species][-1] = "-" * aln_len
        for sp,seqs in alndict.items():
            if sp not in spcheck:
                alndict[sp].append("-" * aln_len)
    for sp,seq in alndict.items():
        #outfile.write(">" + sp + "\n")
        seq = "".join(seq)
        if (seq.count("-") + seq.count("N") + seq.count("n")) <= len(seq)/2:
            print(sys.argv[1],sp,seq)
        #outfile.write(seq + "\n")

maf2fa(sys.argv[1],sys.argv[2])
