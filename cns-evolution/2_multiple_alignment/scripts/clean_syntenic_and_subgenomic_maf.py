import sys
import copy
from collections import defaultdict
import itertools
from Bio import AlignIO
# Usage
# python filter_maf_by_synteny.py zmays_filt.collinearity.blockregions.100kb.isect.best_score.bed panand_10_lt2.maf

def filter_zea(outdict):
    '''
    Only allow maize1 to align to maize1 and maize2 to maize2 in Zea and Tripsacum.
    '''
    onetoone = defaultdict(list)
    # reference genome
    ref_sub = None
    for sp,rows in outdict.items():
        if "zB73v5@" in sp:
            ref_sub = sp.split('.')[0].split('@')[1]
    if ref_sub:
        for sp,rows in outdict.items():
            spname = sp.split('.')[0].split('@')[0]
            spsub = sp.split('.')[0].split('@')[1]
            if spname in zea:
                if spsub == ref_sub:
                    for r in rows:
                        onetoone[sp].append(r)
                elif spname == "tdactn":
                    subsub = ref_sub + "a"
                    newname = spname + "@" + ref_sub
                    # tdactn can have maize1/2 or maize1a/2a but maize1b/2b are ignored
                    if spsub == subsub or spsub == ref_sub:
                        for r in rows:
                            chromscaf = r[1].split('.')[1]
                            newr = ['s'] + [newname+'.'+chromscaf] + r[2:]
                            onetoone[newname].append(newr)
            else:
                for r in rows:
                    onetoone[sp].append(r)

        return(onetoone)


zea = ["tdactn","tdacts","zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm"]
panand = ["aburma","achine","agerar","avirgi","blagur","ccitra","crefra","etrips","hcompr","hconto","irugos","ppanic","rrottb","rtuber","sbicol","smicro","snutan","sscopa","tdactn","tdacts","telega","ttrian","tzopol","udigit","vcuspi","zB73v5","zTIL01","zTIL11","zTIL18","zTIL25","zdiplg","zdiplm","zhuehu","zluxur","znicar"]

spdict = {"Ab":"aburma","Ac":"achine","Ag":"agerar","Av":"avirgi","Bl":"blagur","Cc":"ccitra","Cr":"crefra","Et":"etrips","Hc":"hconto","Hp":"hcompr","Ir":"irugos","Pi":"ppanic","Rr":"rrottb","Rt":"rtuber","Sb":"sbicol","Sm":"smicro","Sn":"snutan","Ss":"sscopa","TdFL":"tdacts","TdKS":"tdactn","Te":"telega","Tt":"ttrian","Tz":"tzopol","Ud":"udigit","Vc":"vcuspi","ZdGi":"zdiplg","ZdMo":"zdiplm","Zh":"zhuehu","Zl":"zluxur","Zn":"znicar","Zv01":"zTIL01","Zv11":"zTIL11","Zx18":"zTIL18","Zx25":"zTIL25","Zm":"zB73v5"}

#zea = ["zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm"]

windows = {}
window_size = 100000
max_chrom_size = 310000000
chrom_list = list(range(0,max_chrom_size+window_size,window_size))
for i in range(10):
    chrom = "chr" + str(i+1)
    for idx,m in enumerate(chrom_list[:-1]):
        win_name = chrom + "@" + str(m) + "@" + str(chrom_list[idx+1])
        windows[win_name] = {}

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


# BED is 0-indexed 
# MAF is 0-indexed

#Zmv5_chrom = {"1":308452471,"2":243675191,"3":238017767,"4":250330460,"5":226353449,"6":181357234,"7":185808916,"8":182411202,"9":163004744,"10":152435371}
bases = ["a","c","t","g","A","C","T","G"]
window_size = 100000
for multiple_alignment in AlignIO.parse(sys.argv[2], "maf"):
    #print("printing a new multiple alignment")
    #print("Alignment length %i" % multiple_alignment.get_alignment_length())
    #print("Alignment depth",multiple_alignment.__len__())
    # if record longer 20bp and has over 4 species alignedi
    if multiple_alignment.__len__() > 0:
        rec_count = 0
        splist = []
        covdic = {}
        outdict = defaultdict(list)
        nadict = defaultdict(list)
        skip_aln = False
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


                nearest = chrom_list[min(range(len(chrom_list)), key = lambda i: abs(chrom_list[i]-start_fwd))]
                lookup = "chr" + chrom + "@" + str(nearest) + "@" + str(nearest+window_size)
                syn_result = windows[lookup]
                if "@NA" not in sp:
                    outdict[sp].append(outline)
                else:
                    skip_aln = True
            else:
                rec_count += 1
                if "@NA" not in sp:
                    if sp in outdict:
                        if '@' in sp:
                            spname = sp.split('@')[0]
                            syntenic_regions = syn_result[spname]
                        else:
                            syntenic_regions = syn_result[sp]
                        #print("Multiple entries for",sp,syntenic_regions)
                        #print("Current entry:", outdict[sp])
                        #print("New entry:", record)
                        cur_syntenic_genes_best = 0
                        new_syntenic_genes_best = 0
                        cur_chrom = "".join(outdict[sp][0][1].split('.')[1:])
                        if outdict[sp][0][4] == "+":
                            cur_start = outdict[sp][0][2]
                            cur_end = outdict[sp][0][2] + outdict[sp][0][3]
                        else:
                            cur_start = outdict[sp][0][5] - (outdict[sp][0][2] + outdict[sp][0][3])
                            cur_end = outdict[sp][0][5] - outdict[sp][0][2]
                        for region in syntenic_regions:
                            #'chr10@129963@126029522@Alignment64477@chr10@137403@119686399@1131'
                            region = region.split('@')
                            if cur_chrom == region[4] and cur_start >= int(region[5]) and cur_end <= int(region[6]):
                                if int(region[7]) > cur_syntenic_genes_best:
                                    cur_syntenic_genes_best = int(region[7])
                            if chrom == region[4] and start_fwd >= int(region[5]) and end_fwd <= int(region[6]):
                                if int(region[7]) > new_syntenic_genes_best:
                                    new_syntenic_genes_best = int(region[7])
                        if new_syntenic_genes_best > cur_syntenic_genes_best:
                            outdict[sp] = [outline]
                            #print("New record is better!",new_syntenic_genes_best,cur_syntenic_genes_best)
                        #else:
                        #    print("Current record is better",new_syntenic_genes_best,cur_syntenic_genes_best)

                        # The below check has been replaced with a better check that ensures the selected alignment has the most syntenic genes supporting it
                        # check if this alignment is on a longer scaf
                        #current_size = outline[5]
                        #if current_size < src_size:
                        #    outdict[sp].append(outline)
                    else:
                        outdict[sp].append(outline)
                else:
                    nadict[sp].append(outline)


                #outdict[sp].append(outline)
        # Check for entirely missing species
        for species in panand:
            species_absent = True
            for sp in outdict:
                spname = sp.split('@')[0]
                if spname == species:
                    species_absent = False
            if species_absent:
                # check the NA dict
                species_na_name = species + "@NA"
                if species_na_name in nadict:
                    # find the na alignment from longest scaf
                    longest_na = -1
                    longest_line = None
                    for i in nadict[species_na_name]:
                        i_len = i[3]
                        if i_len > longest_na:
                            longest_line = i
                            longest_na = i[3]
                    species_new_name = species + "@1"
                    longest_line = ['s'] + [species_new_name + "."+ longest_line[1].split('.')[1]] + longest_line[2:]
                    outdict[species_new_name] = [longest_line]
        #print("outdict",outdict)
        #print("nadict",nadict)
        if len(outdict) > 1 and not skip_aln:
            outdict = filter_zea(outdict)
            print("")
            print("a")
            for k,v in outdict.items():
                if len(v) == 1:
                    print(*v[0],sep='\t')
