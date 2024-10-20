import sys
import copy
from collections import defaultdict
import itertools
from Bio import AlignIO
# Usage
# python filter_maf_by_synteny.py zmays_filt.collinearity.blockregions.100kb.isect.best_score.bed panand_10_lt2.maf

zea = ["tdactn","tdacts","zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm"]
panand = ["aburma","achine","agerar","avirgi","blagur","ccitra","crefra","etrips","hcompr","hconto","irugos","ppanic","rrottb","rtuber","sbicol","smicro","snutan","sscopa","tdactn","tdacts","telega","ttrian","tzopol","udigit","vcuspi","zB73v5","zTIL01","zTIL11","zTIL18","zTIL25","zdiplg","zdiplm","zhuehu","zluxur","znicar"]

spdict = {"Ab":"aburma","Ac":"achine","Ag":"agerar","Av":"avirgi","Bl":"blagur","Cc":"ccitra","Cr":"crefra","Et":"etrips","Hc":"hconto","Hp":"hcompr","Ir":"irugos","Pi":"ppanic","Rr":"rrottb","Rt":"rtuber","Sb":"sbicol","Sm":"smicro","Sn":"snutan","Ss":"sscopa","TdFL":"tdacts","TdKS":"tdactn","Te":"telega","Tt":"ttrian","Tz":"tzopol","Ud":"udigit","Vc":"vcuspi","ZdGi":"zdiplg","ZdMo":"zdiplm","Zh":"zhuehu","Zl":"zluxur","Zn":"znicar","Zv01":"zTIL01","Zv11":"zTIL11","Zx18":"zTIL18","Zx25":"zTIL25","Zm":"zB73v5"}

#zea = ["zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm"]

Zmv5_chrom = {"1":308452471,"2":243675191,"3":238017767,"4":250330460,"5":226353449,"6":181357234,"7":185808916,"8":182411202,"9":163004744,"10":152435371}

bases = ["a","c","t","g","A","C","T","G"]
for multiple_alignment in AlignIO.parse(sys.argv[1], "maf"):
    #print("printing a new multiple alignment")
    #print("Alignment length %i" % multiple_alignment.get_alignment_length())
    #print("Alignment depth",multiple_alignment.__len__())
    # if record longer 20bp and has over 4 species alignedi
    if multiple_alignment.__len__() > 0:
        rec_count = 0
        splist = []
        outlist = []
        bad_block = False
        for record in multiple_alignment:
            sp = str(record.id).split(".")[0]
            chrom = "".join(str(record.id).split(".")[1:])
            # start pos
            start = int(record.annotations["start"])
            # ungapped length
            size = int(record.annotations["size"])
            # MAF is 0-indexed and end is non-inclusive like python
            end = start + size
            # get non "N" positions
            refseq = str(record.seq)
            for character in refseq:
                if character not in bases and character not in ["N","n","-"]:
                    bad_block = True
            # strand
            if int(record.annotations["strand"]) ==1:
                strand = "+"
            else:
                strand = "-"
            #strand = int(record.annotations["strand"])
            src_size = int(record.annotations["srcSize"])
            outline = ["s",record.id,start,size,strand,src_size,refseq]

            if strand == "-":
                start_fwd = src_size - end
                end_fwd = src_size - start
                #print("Strand was negative, coords corrected to ",start_fwd, end_fwd, "for ", outline)
            else:
                start_fwd = start
                end_fwd = end
            rec_count += 1
            if sp not in splist:
                splist.append(sp)
            else:
                bad_block = True
            outlist.append(outline)
        if bad_block:
            print("Redundant species in block!",outlist)
