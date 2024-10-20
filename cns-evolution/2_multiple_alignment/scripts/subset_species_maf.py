import sys
from Bio import AlignIO
# Usage
species = sys.argv[2]
species = species.split(",")


for multiple_alignment in AlignIO.parse(sys.argv[1], "maf"):
    #print("printing a new multiple alignment")
    #print("Alignment length %i" % multiple_alignment.get_alignment_length())
    #print("Alignment depth",multiple_alignment.__len__())
    # if record longer 20bp and has over 4 species alignedi
    rec_count = 0
    for record in multiple_alignment:
        sp = str(record.id).split(".")[0]
        if '@' in sp:
            sp = str(record.id).split("@")[0]
        chrom = "".join(str(record.id).split(".")[1:])
        # start pos
        start = int(record.annotations["start"])
        # ungapped length
        size = int(record.annotations["size"])
        # MAF is 0-indexed and end is non-inclusive like python
        end = start + size
        # get non "N" positions
        refseq = str(record.seq)
        # strand
        if int(record.annotations["strand"]) ==1:
            strand = "+"
        else:
            strand = "-"
        #strand = int(record.annotations["strand"])
        src_size = int(record.annotations["srcSize"])
        outline = ["s",record.id,start,size,strand,src_size,refseq]

        if rec_count == 0:
            #assume first record is ref
            print("")
            print("a")
            # uncomment below to replace refseq with N
            #refseq = "N" * len(refseq)
            outline = ["s",record.id,start,size,strand,src_size,refseq]
            print(*outline,sep='\t')
            rec_count += 1
            # Get closest range based on 100k ranges
            #chrom_size = Zmv5_chrom[chrom]
        else:
            if sp in species:
                outline = ["s",record.id,start,size,strand,src_size,refseq]
                rec_count += 1
                print(*outline,sep='\t')

