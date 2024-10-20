import sys
import copy
from collections import defaultdict
import itertools
from Bio import AlignIO
# Usage
# python filter_maf_by_synteny.py zmays_filt.collinearity.blockregions.100kb.isect.best_score.bed panand_10_lt2.maf
# zluxur is missing as it has not been annotated due to lack of mRNA
spdict = {"Ab":"aburma","Ac":"achine","Ag":"agerar","Av":"avirgi","Bl":"blagur","Cc":"ccitra","Cr":"crefra","Et":"etrips","Hc":"hconto","Hp":"hcompr","Ir":"irugos","Pi":"ppanic","Rr":"rrottb","Rt":"rtuber","Sb":"sbicol","Sm":"smicro","Sn":"snutan","Ss":"sscopa","TdFL":"tdacts","TdKS":"tdactn","Te":"telega","Tt":"ttrian","Tz":"tzopol","Ud":"udigit","Vc":"vcuspi","ZdGi":"zdiplg","ZdMo":"zdiplm","Zh":"zhuehu","Zl":"zluxur","Zn":"znicar","Zv01":"zTIL01","Zv11":"zTIL11","Zx18":"zTIL18","Zx25":"zTIL25","Zm":"zB73v5"}

tripsacineae = ["tdactn","tdacts","zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm"]


## PREPARE WINDOWED REGIONS DICT
# Windows allow faster lookup

windows = {}
window_size = 100000
max_chrom_size = 310000000
chrom_list = list(range(0,max_chrom_size+window_size,window_size))
for i in range(10):
    chrom = "chr" + str(i+1)
    for idx,m in enumerate(chrom_list[:-1]):
        win_name = chrom + "@" + str(m) + "@" + str(chrom_list[idx+1])
        windows[win_name] = {}


## PREPARE COLLINEAR REGIONS DICT

syndict = {}
with open(sys.argv[1],'r') as window:
    for l in window:
        l = l.strip().split('\t')
        # Add alignment ID to region to ensure uniqueness
        region = '@'.join([l[3].split('.')[1]] + l[4:7])
        region_chrom = l[3].split('.')[1]
        if not 'scaf' in region_chrom:
            region_start = int(l[4])
            region_end = int(l[5]) 
            nearest_start = max(chrom_list[min(range(len(chrom_list)), key = lambda i: abs(chrom_list[i]-region_start))] - window_size,0)
            nearest_end = chrom_list[min(range(len(chrom_list)), key = lambda i: abs(chrom_list[i]-region_end))] + (window_size*2)
            region_chrom_list = list(range(nearest_start,nearest_end,window_size))
            # Now add syn dict info for each
            lspecies = l[0].split('.')[0]
            #lspecies = spdict[species]
            qchrom = '_'.join(l[0].split('.')[1:])
            qregion = '@'.join([qchrom] + l[1:3] + [l[9]])
            syn_info_line = region + "@" + qregion
            #print("INFO:",lspecies,region_chrom,syn_info_line,l,region_chrom_list)
            for idx,m in enumerate(region_chrom_list[:-1]):
                win_name = region_chrom + "@" + str(m) + "@" + str(region_chrom_list[idx+1])
                if lspecies not in windows[win_name]:
                    windows[win_name][lspecies] = []
                windows[win_name][lspecies].append(syn_info_line)
                #print("Added:",win_name,windows[win_name][lspecies])
# window file example rows
#Zm_chr10        94075310        98308866        Ab_ctg_1        1784911 3599972 Alignment0      1223.0  7.3e-87 26
#Zm_chr2 117129112       118210875       Ab_ctg_1        3205240 3466389 Alignment1      581.0   0  12

## PREPARE SUBGENOME DICT

subdict = {}
with open(sys.argv[2],'r') as window:
    for l in window:
        l = l.strip().split('\t')
        species = l[0].split('.')[0]
        qchrom = '_'.join(l[0].split('.')[1:])
        # Add alignment ID to region to ensure uniqueness
        region = '@'.join([qchrom,l[1],l[2],l[6],l[9]])
        subgenome = l[12]
        if species not in subdict:
            subdict[species] = {}
        if qchrom not in subdict[species]:
            subdict[species][qchrom] = {}
        subdict[species][qchrom][region] = subgenome
        #if species == "zB73v5":
        #    print("Added B73 :", species, region,subgenome)
# BED is 0-indexed 
# MAF is 0-indexed

#Zmv5_chrom = {"1":308452471,"2":243675191,"3":238017767,"4":250330460,"5":226353449,"6":181357234,"7":185808916,"8":182411202,"9":163004744,"10":152435371}
bases = ["a","c","t","g","A","C","T","G"]
window_size = 100000
for multiple_alignment in AlignIO.parse(sys.argv[3], "maf"):
    #print("printing a new multiple alignment")
    #print("Alignment length %i" % multiple_alignment.get_alignment_length())
    #print("Alignment depth",multiple_alignment.__len__())
    if multiple_alignment.__len__() > 0:
        rec_count = 0
        splist = []
        covdic = {}
        synteny_info = None
        outdict = defaultdict(list)
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

            if rec_count == 0:
                #assume first record is ref
                rec_count += 1
                ref_chrom = "chr" + "".join(str(record.id).split(".")[1:])
                # Get closest range based on 100k ranges
                #chrom_size = Zmv5_chrom[chrom]
                chrom_size = src_size
                #chrom_list = list(range(0,chrom_size,window_size))
                nearest = chrom_list[min(range(len(chrom_list)), key = lambda i: abs(chrom_list[i]-start_fwd))]
                lookup = ref_chrom + "@" + str(nearest) + "@" + str(nearest+window_size)
                syn_result = windows[lookup]
                # WARNING : chrom prefix hardcoded below!
                overlap_range = "chr" + chrom + "@" + str(start_fwd) + "@" + str(end_fwd)
                sub_path = "NA"
                if sp in subdict:
                    if ref_chrom in subdict[sp]:
                        for sub_region,sub_genome in subdict[sp][ref_chrom].items():
                            #print("Checking aln", query_range, "in",sub_region,sub_genome)
                            sub_region_split = sub_region.split('@')
                            sub_chrom = sub_region_split[0]
                            sub_start = int(sub_region_split[1])
                            sub_end = int(sub_region_split[2])
                            sub_overlap = range(max(start_fwd, sub_start), min(end_fwd, sub_end)+1)
                            if len(sub_overlap)>0 and sub_genome != "NA":
                                sub_path = subdict[sp][ref_chrom][sub_region]
                    sub_sp = outline[1].split('.')[0] + '@' + str(sub_path) + '.' + outline[1].split('.')[1]
                    sub_outline = ['s',sub_sp] + outline[2:]
                    outdict[sub_sp].append(sub_outline)



                #outdict[sp].append(outline)
                
                
                # COLLECT INFO 
                #if overlap_range in syndict:
                #    synteny_info = syndict[overlap_range]
                #    outdict[sp].append(outline)
            else:
                rec_count += 1
                query_range = chrom + "@" + str(start_fwd) + "@" + str(end_fwd)
                is_syntenic = False
                # 1) IS QUERY SYNTENIC TO REFERENCE?
                # ignore self-alignments
                if sp != "zB73v5":
                    # would make sense to have them grouped in subdict by chrom to speed up
                    if sp in syn_result:
                        for a in syn_result[sp]:
                            # does query_range intersect
                            a = a.split('@')
                            query_chrom = a[4]
                            query_start = int(a[5])
                            query_end = int(a[6])
                            query_syn_genes = int(a[7])
                            #if sp == "zluxur" or sp == "znicar":
                            #    query_syn_gene_cutoff = 5
                            #else:
                            #    query_syn_gene_cutoff = 10
                            query_syn_gene_cutoff = 5
                            query_overlap = range(max(start_fwd, query_start), min(end_fwd, query_end)+1)
                            if query_chrom == chrom and len(query_overlap)>0 and query_syn_genes > query_syn_gene_cutoff:
                                is_syntenic = True
                        if is_syntenic:
                            sub_path = "NA"
                            if sp in subdict:
                                #print("Check if ", chrom,"is in: ",subdict[sp])
                                if chrom in subdict[sp]:
                                    for sub_region,sub_genome in subdict[sp][chrom].items():
                                        #print("Checking aln", query_range, "in",sub_region,sub_genome)
                                        sub_region_split = sub_region.split('@')
                                        sub_chrom = sub_region_split[0]
                                        sub_start = int(sub_region_split[1])
                                        sub_end = int(sub_region_split[2])
                                        # if not Zea then assign subgenome from scaffold alone
                                        if sp not in tripsacineae:
                                            #print("Nonzea check",sp,query_range,sub_region,sub_genome)
                                            if sub_chrom == chrom and sub_genome != "NA":
                                                sub_path = subdict[sp][chrom][sub_region]
                                        else:
                                            sub_overlap = range(max(start_fwd, sub_start), min(end_fwd, sub_end)+1)
                                            if sub_chrom == chrom and len(sub_overlap)>0 and sub_genome != "NA":
                                                sub_path = subdict[sp][chrom][sub_region]
                                sub_sp = outline[1].split('.')[0] + '@' + str(sub_path) + '.' + outline[1].split('.')[1]
                                sub_outline = ['s',sub_sp] + outline[2:]
                                outdict[sub_sp].append(sub_outline)
        if len(outdict) > 1:
            print("")
            print("a")
            for k,v in outdict.items():
                #if len(v) == 1:
                print(*v[0],sep='\t')

                            #print("Lookup:",chrom,overlap_range,sp,syn_result[sp])
#                    for ref_reg,q_reg in syndict[sp].items():
#                        ref_chrom = overlap_range.split('@')[0]
#                        q_ref_chrom = ref_reg.split('@')[0]
#                        if ref_chrom == q_ref_chrom:
#                            print("Same chrom",sp, overlap_range,ref_reg,q_reg)
#                            ref_start = int(overlap_range.split('@')[1])
#                            ref_end = int(overlap_range.split('@')[2])
#                            qref_start = int(ref_reg.split('@')[1])
#                            qref_end = int(ref_reg.split('@')[2])
#                            ref_qref_olap = range(max(ref_start, qref_start), min(ref_end, qref_end)+1)
#                            if len(ref_qref_olap) > 0:
#                                print("Found overlap overlap_range,ref_reg,q_reg:",sp,overlap_range,ref_reg,q_reg)
#                            else:
#                                print("No overlap:",len(ref_qref_olap),ref_start,ref_end,qref_start,qref_end)
                    # 2) WHICH QUERY SUBGENOME IS THIS?
                    

                    #if synteny_info and sp in synteny_info:
                    #    sp_syn = synteny_info[sp].split("@")
                    #    sp_syn_chr = sp_syn[0]
                    #    sp_syn_start = int(sp_syn[1])
                    #    sp_syn_end = int(sp_syn[2])
                    #    if chrom == sp_syn_chr and sp_syn_start <= start_fwd and sp_syn_end >= end_fwd:
                    #        outdict[sp].append(outline)
        #if synteny_info:
        #    #for species in outdict:
        #    if len(outdict) > 3:
        #        print("")
        #        print("a")
        #        for k,v in outdict.items():
        #            if len(v) == 1:
        #                print(*v[0],sep='\t')


