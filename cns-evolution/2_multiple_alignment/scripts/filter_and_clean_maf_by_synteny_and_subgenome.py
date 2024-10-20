import sys
import copy
from collections import defaultdict
import itertools
import gzip
from Bio import AlignIO
# Usage
# python filter_and_clean_maf_by_synteny_and_subgenome.py zmays_filt.collinearity.blockregions.100kb.isect.best_score.bed panand_10_lt2.maf remove_non_syntenic single|all

def remove_non_syntenic(synteny_list,ref_chrom_pos):
    '''
    Remove syntenic alignments that dont overlap the reference region.
    '''
    overlapping = []
    block_ref_chrom,block_ref_start,block_ref_end = ref_chrom_pos
    for aln in synteny_list:
        #print(aln)
        aln_ref_chrom = aln[0]
        aln_ref_start = int(aln[1])
        aln_ref_end = int(aln[2])
        if block_ref_chrom == aln_ref_chrom:
            overlap_range = range(max(block_ref_start, aln_ref_start), min(block_ref_end, aln_ref_end)+1)
            if len(overlap_range)>0:
                overlapping.append(aln)
        # keep alignments with no synteny information
        #elif aln_ref_chrom == 'NA':
        #    overlapping.append(aln)
    return overlapping

def filter_aln(annotated_aln):
    '''
    Read in a species alignment dict with subgenomes as keys. Each subgenome value is a nested dict containing one to many alignments with their synteny information. Return a nested list of one highest scoring alignment per subgenome.
    '''
    filt_dict = {}
    #1 if there are multiple aln for a single subgenome, sort by most syntenic genes and pick top
    for subgenome,aln in annotated_aln.items():
        filt_aln = sorted(aln, key=lambda x: int(x[9]),reverse=True)[0]
        filt_dict[subgenome] = filt_aln
    return filt_dict

zea = ["zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm","tdactn","tdacts"]



## Read in collinearity file ##

# Prepare an indexed dict to store collinear regions
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

## Read in subgenome assignments ##

subdict = {}
with open(sys.argv[2],'r') as window:
    for l in window:
        l = l.strip().split('\t')
        species = l[0].split('.')[0]
        qchrom = '_'.join(l[0].split('.')[1:])
        # Add alignment ID to region to ensure uniqueness
        region = '@'.join([qchrom,l[1],l[2],l[6]])
        subgenome = l[12]
        if species not in subdict:
            subdict[species] = {}
        if qchrom not in subdict[species]:
            subdict[species][qchrom] = {}
        subdict[species][qchrom][region] = subgenome

## Parse MAF and print only syntenic alignments ##

# 1) Alignments are ranked by synteny score and the top scoring alignment is selected after subgenome filtering
# 2) If no alignment is in a matched subgenome (Tripsacineae) or has any subgenome information, then the highest ranked by synteny is selected
# 3) Non-syntenic alignments are discarded
# 4) For vcuscpi (and possibly any other species), if the only available alignments are from subgenome "NA" then the top ranked based on synteny can be selected

if sys.argv[3].endswith('.gz'):
    maf_parser = AlignIO.parse(gzip.open(sys.argv[3],"rt"), "maf")
else:
    maf_parser = AlignIO.parse(sys.argv[3],"maf")

print("##maf version=1")
for multiple_alignment in maf_parser:
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
            #print("RECORD:", sp, chrom,start,size,end,refseq,strand)
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

            # For reference sequence
            if rec_count == 0:
                #assume first record is ref
                rec_count += 1
                ref_chrom = "chr" + "".join(str(record.id).split(".")[1:])
                ref_chrom_pos = [ref_chrom,start_fwd,end_fwd]
                # Get closest range based on 100k ranges
                nearest = chrom_list[min(range(len(chrom_list)), key = lambda i: abs(chrom_list[i]-start_fwd))]
                lookup = ref_chrom + "@" + str(nearest) + "@" + str(nearest+window_size)
                # Look up all collinear regions
                syn_result = windows[lookup]
                # WARNING : chrom prefix 'chr' hardcoded below!
                overlap_range = "chr" + chrom + "@" + str(start_fwd) + "@" + str(end_fwd)
                # Assign subgenome to reference
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
                            # assumes subgenome compartments are non-overlapping
                            # regions spanning a boundary will be assigned to the last matching subgenome in the list
                            if len(sub_overlap)>0 and sub_genome != "NA":
                                sub_path = subdict[sp][ref_chrom][sub_region]
                    sub_sp = outline[1].split('.')[0] + '@' + str(sub_path) + '.' + outline[1].split('.')[1]
                    sub_outline = ['s',sub_sp] + outline[2:]
                    outdict[sub_sp].append(sub_outline)
                #print("Analysing alignment at", sp,chrom,start_fwd,end_fwd,sub_path)
            else:
                rec_count += 1
                query_range = chrom + "@" + str(start_fwd) + "@" + str(end_fwd)
                is_syntenic = False
                synteny_hits = []
                # ignore self-alignments
                if sp != "zB73v5":
                    if sp in syn_result:
                        for a in syn_result[sp]:
                            # does query_range intersect
                            a = a.split('@')
                            query_chrom = a[4]
                            query_start = int(a[5])
                            query_end = int(a[6])
                            query_syn_genes = int(a[7])
                            query_overlap = range(max(start_fwd, query_start), min(end_fwd, query_end)+1)
                            #print("Checking overlap",query_chrom,query_start,query_end)
                            #print("Alignment region",chrom, start_fwd, end_fwd)
                            if query_chrom == chrom and len(query_overlap)>0:
                                is_syntenic = True
                                synteny_hits.append(a)
                        if is_syntenic:
                            sub_path = "NA"
                            if sp in subdict:
                                if chrom in subdict[sp]:
                                    for sub_region,sub_genome in subdict[sp][chrom].items():
                                        sub_region_split = sub_region.split('@')
                                        sub_chrom = sub_region_split[0]
                                        sub_start = int(sub_region_split[1])
                                        sub_end = int(sub_region_split[2])
                                        # if not Zea then assign subgenome from scaffold alone
                                        if sp not in zea:
                                            if sub_chrom == chrom and sub_genome != "NA":
                                                sub_path = subdict[sp][chrom][sub_region]
                                        else:
                                            sub_overlap = range(max(start_fwd, sub_start), min(end_fwd, sub_end)+1)
                                            if len(sub_overlap)>0 and sub_genome != "NA":
                                                sub_path = subdict[sp][chrom][sub_region]
                                outline = outline + [synteny_hits] + [sub_path]
                                outdict[sp].append(outline)
                        else:
                            # Allow nonsyntenic alignments to be retained
                            nonsyntenic = [['NA', '0', '0', 'NA', 'NA', '0', '0', '0']]
                            outline = outline + [nonsyntenic] + ['NA']
                            outdict[sp].append(outline)
        if len(outdict) > 1:
            #print("")
            #print("a")
            mysubdict = defaultdict(lambda: defaultdict(list))
            for k,v in outdict.items():
                #if len(v) == 1:
                #print(*v[0],sep='\t')
                #print(k,v)
                #syn_aln = filter_aln(v)
                if not "zB73v5" in k:
                    v_ranked = []
                    #print(k,v)
                    for i in v:
                        #print(i)
                        synteny_list = i[7]
                        # Remove alignments not syntenic with the reference region
                        if sys.argv[4] == "remove_non_syntenic":
                            #print("Prefilter", synteny_list,ref_chrom_pos)
                            synteny_list = remove_non_syntenic(synteny_list,ref_chrom_pos)
                            #print("Postfilter", synteny_list)
                        if any(synteny_list):
                            # Append the best syntenic score from 1 to many scores
                            best_synteny_score = sorted(synteny_list, key=lambda x: int(x[7]),reverse=True)[0][7]
                            i_ranked = i[:7] + [synteny_list,i[8],best_synteny_score]
                            v_ranked.append(i_ranked)
                    # for each subgenome store the best match
                    for ir in v_ranked:
                        subname = k + "@" + ir[8]
                        mysubdict[k][subname].append(ir)
                    #for i in v:
                    #    print(i)
                    #for suba in v:
                    #    if len(suba[7]) >1:
                    #        print(v)
                    #        print("Multiple syntenic:", suba)
                else:
                    ref_subgenome = v[0][1].split('@')[1].split('.')[0]
                    print("")
                    print("a")
                    print(*v[0],sep='\t')
            for k,v in mysubdict.items():
                # print the best subgenome alignment for each species
                #print(k,v)
                best_hits = filter_aln(v)
                #print(best_hits)
                # Drop unassigned subgenome unless no other option
                #if list(best_hits.keys()).count(k+"@NA") < len(list(best_hits.keys())):
                #    best_hits.pop(k+"@NA", None)
                if list(v.keys()).count(k+"@NA") < len(list(v.keys())):
                    v.pop(k+"@NA", None)
                best_hits = filter_aln(v)
                # If the reference is unassigned to a subgenome
                if ref_subgenome == "NA":
                    # for trips/zea just grab the most syntenic subgenome
                    #sort highest synteny score
                    if k in zea:
                        best_synteny_count = 0
                        best_subgenome = list(best_hits.keys())[0]
                        for subgenome,aln in best_hits.items():
                            #synteny_score = aln[7]
                            synteny_count = int(aln[7][0][7])
                            if synteny_count > best_synteny_count:
                                best_synteny_count = synteny_count
                                best_subgenome = subgenome
                        if sys.argv[5] == "single":
                            aln = best_hits[best_subgenome]
                            chrom = aln[1].split('.')[1]
                            subchrom = best_subgenome+'.'+chrom
                            out_aln = [aln[0],subchrom] + aln[2:7]
                            print(*out_aln,sep='\t')
                        elif sys.argv[5] == "all":
                            all_aln = v[best_subgenome]
                            for a in all_aln:
                                chrom = a[1].split('.')[1]
                                subchrom = best_subgenome+'.'+chrom
                                out_aln = [a[0],subchrom] + a[2:7]
                                print(*out_aln,sep='\t')
                            # print every single syntenic alignment
                    #for others just print all subgenomes
                    else:
                        if sys.argv[5] == "single":
                            for subgenome,aln in best_hits.items():
                                query_subgenome = subgenome.split('@')[1]
                                chrom = aln[1].split('.')[1]
                                subchrom = subgenome+'.'+chrom
                                out_aln = [aln[0],subchrom] + aln[2:7]
                                print(*out_aln,sep='\t')
                        elif sys.argv[5] == "all":
                        # To do: add code here to print all aln per subgenome
                            for subgenome,aln in v.items():
                                query_subgenome = subgenome.split('@')[1]
                                for a in aln:
                                    chrom = a[1].split('.')[1]
                                    subchrom = subgenome+'.'+chrom
                                    out_aln = [a[0],subchrom] + a[2:7]
                                    print(*out_aln,sep='\t')
                else:
                    if sys.argv[5] == "single":
                        for subgenome,aln in best_hits.items():
                            query_subgenome = subgenome.split('@')[1]
                            chrom = aln[1].split('.')[1]
                            subchrom = subgenome+'.'+chrom
                            out_aln = [aln[0],subchrom] + aln[2:7]
                            # if Trips then only print matching subgenome
                            if k in zea:
                                if query_subgenome == ref_subgenome or query_subgenome == "NA":
                                    print(*out_aln,sep='\t')
                            else:
                                print(*out_aln,sep='\t')
                    elif sys.argv[5] == "all":
                    # To do: add code here to print all aln per subgenome
                        for subgenome,aln in v.items():
                            query_subgenome = subgenome.split('@')[1]
                            for a in aln:
                                chrom = a[1].split('.')[1]
                                subchrom = subgenome+'.'+chrom
                                out_aln = [a[0],subchrom] + a[2:7]
                                if k in zea:
                                    if query_subgenome == ref_subgenome or query_subgenome == "NA":
                                        print(*out_aln,sep='\t')
                                else:
                                    print(*out_aln,sep='\t')
                # if "NA" subgenome, only print if there is no other subgenome
