import sys
import random
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord

def trim_fa(infile,max_len,outfa):
    infa = SeqIO.parse(infile, "fasta")
    # if iterating infa object, record will be skipped, thus use new FastaIterator
    seq_len = len(next(SeqIO.parse(infile, "fasta")))
    seq_cut = []
    print(infa)
    print("seqlen",seq_len)
    if seq_len > max_len:
        # sample random indices
        random_indices = random.sample(range(0, seq_len + 1), max_len)
        for alignment in infa:
            # Get max_len seq from start of alignment
            #subseq = alignment.seq[0:max_len]
            subseq = Seq.Seq(''.join(list(map(alignment.seq.__getitem__, random_indices))))
            # Sampling random sites or blocks substantially slows down process
            record = SeqRecord(subseq,alignment.id,"","")
            seq_cut.append(record)
    else:
        seq_cut = infa
    SeqIO.write(seq_cut, outfa, "fasta")
infasta = sys.argv[1]
outfasta = sys.argv[2]
sites = int(sys.argv[3])
trim_fa(infasta,sites,outfasta)
