import sys
import re

def tag2list(aln_tag):
    #Alignment4:score=963.0e_value=2e-64N=21sb1&zm1plus
    blockname = aln_tag.split(':')[0].replace("##","")
    aln_score = aln_tag.split('=')[1].replace("e_value","")
    aln_evalue = aln_tag.split('=')[2].replace("N","")
    aln_gene_string = aln_tag.split('=')[3]
    aln_gene_count = re.findall("[\d]*", aln_gene_string)[0]
    return blockname,aln_score,aln_evalue,aln_gene_count

spdict = {"Ab":"aburma","Ac":"achine","Ag":"agerar","Av":"avirgi","Bl":"blagur","Cc":"ccitra","Cr":"crefra","Et":"etrips","Hc":"hconto","Hp":"hcompr","Ir":"irugos","Pi":"ppanic","Rr":"rrottb","Rt":"rtuber","Sb":"sbicol","Sm":"smicro","Sn":"snutan","Ss":"sscopa","TdFL":"tdacts","TdKS":"tdactn","Te":"telega","Tt":"ttrian","Tz":"tzopol","Ud":"udigit","Vc":"vcuspi","ZdGi":"zdiplg","ZdMo":"zdiplm","Zh":"zhuehu","Zl":"zluxur","Zn":"znicar","Zv01":"zTIL01","Zv11":"zTIL11","Zx18":"zTIL18","Zx25":"zTIL25","Zm":"zB73v5"}

gffdict = {}
with open(sys.argv[1],'r') as gff:
    for l in gff:
        l = l.strip()
        l = l.split('\t')
        gffdict[l[1]] = {'chr':l[0],'start':int(l[2]),'end':int(l[3])}


with open(sys.argv[2],'r') as blocks:
    for l in blocks:
        l = l.strip()
        l = l.split('\t')
        blockname,aln_score,aln_evalue,aln_gene_count = tag2list(l[0])
        chr1_start = gffdict[l[2]]['chr']
        chr2_start =  gffdict[l[3]]['chr']
        start1 = gffdict[l[2]]['start']
        start2 = gffdict[l[3]]['start']
        end1 = gffdict[l[6]]['end']
        end2 = gffdict[l[7]]['end']
        chr1_end = gffdict[l[6]]['chr']
        chr2_end = gffdict[l[7]]['chr']
        chr1_strand = "+"
        chr2_strand = "+"
        if chr1_start == chr1_end and chr2_start == chr2_end:
            # report only on forward strand
            if int(start1) > int(end1):
                start1,end1 = end1,start1
                chr1_strand = "-"
            if int(start2) > int(end2):
                start2,end2 = end2,start2
                chr2_strand = "-"
            # ensure that reference genome comes first
            if "Zm_" in chr1_start:
                chr1_start,start1,end1,chr1_strand,chr2_start,start2,end2,chr2_strand = chr2_start,start2,end2,chr2_strand,chr1_start,start1,end1,chr1_strand
            chr1_start = chr1_start.replace("_",".",1)
            chr2_start = chr2_start.replace("_",".",1)
            chr1_start = spdict[chr1_start.split(".")[0]] + "." + chr1_start.split(".")[1]
            chr2_start = spdict[chr2_start.split(".")[0]] + "." + chr2_start.split(".")[1]
            print(*[chr1_start,start1,end1,chr2_start,start2,end2,blockname,aln_score,aln_evalue,aln_gene_count,chr1_strand,chr2_strand],sep='\t')
        else:
            print(f'Chr not matching for block {blockname}')

