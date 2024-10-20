import sys
import random
#
def interval_dist(r1, r2):
     x, y = sorted((r1, r2))
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0


## INPUT FILES ##

# collinear genes and blocks based on zmv5 public
collinear = "../data/collinear_genes_canonical.tsv"
collinear_block = "../data/alignment_information.tsv"
# subgenome information
subgenomes = "../../1_subgenomes/data/panand_subgenomes.txt"
# old get only tdactn and tdacts from here
ortho_old = "../data/Orthogroups_v7.txt"
# new but without tdactn and tdacts
ortho_new = "../data/Orthogroups.tsv"
#-collinear blocks based on panand annotation (not public zmv5)
col_block_old = "../data/zmays_filt.collinearity.blockregions.bed"
# All syntenic blocks for each zmv5 gene in each panand genome
syn_gene = "../data/all_gene_synteny_merged.bed"

gene_locations = '../data/panand_gene_locations.bed'
tandem_dups = '../data/zB73v5_tandem_duplicates.txt'
 

windows = {}
window_size = 100000
max_chrom_size = 310000000

chrom_list = list(range(0,max_chrom_size+window_size,window_size))
for i in range(10):
    chrom = "chr" + str(i+1)
    for idx,m in enumerate(chrom_list[:-1]):
        win_name = chrom + "@" + str(m) + "@" + str(chrom_list[idx+1])
        windows[win_name] = {}

spdict = {"Ab":"aburma","Ac":"achine","Ag":"agerar","Av":"avirgi","Bl":"blagur","Cc":"ccitra","Cr":"crefra","Et":"etrips","Hc":"hconto","Hp":"hcompr","Ir":"irugos","Pi":"ppanic","Rr":"rrottb","Rt":"rtuber","Sb":"sbicol","Sm":"smicro","Sn":"snutan","Ss":"sscopa","TdFL":"tdacts","TdKS":"tdactn","Te":"telega","Tt":"ttrian","Tz":"tzopol","Ud":"udigit","Vc":"vcuspi","ZdGi":"zdiplg","ZdMo":"zdiplm","Zh":"zhuehu","Zl":"zluxur","Zn":"znicar","Zv01":"zTIL01","Zv11":"zTIL11","Zx18":"zTIL18","Zx25":"zTIL25","Zm":"zB73v5"}

gffdict = {"Zm-B73-REFERENCE-NAM-5.0":"zB73v5","Ab-Traiperm_572-DRAFT":"aburma","Ac-Pasquet1232-DRAFT":"achine","Ag-CAM1351-DRAFT":"agerar","Av-Kellogg1287_8-REFERENCE":"avirgi","Bl-K1279B-DRAFT":"blagur","Cc-PI314907-DRAFT":"ccitra","Cr-AUB069-DRAFT":"crefra","Et-Layton_Zhong168-DRAFT":"etrips","Hc-AUB53_1-DRAFT":"hconto","Hp-KelloggPI404118-DRAFT":"hcompr","Ir-Pasquet1136-DRAFT":"irugos","Pi-Clark-DRAFT":"ppanic","Rr-Malcomber3106-DRAFT":"rrottb","Rt-Layton_Zhong169-DRAFT":"rtuber","Sbicolor_454_v3.1.1":"sbicol","Sm-PI203595-DRAFT":"smicro","Sn-CAM1369-DRAFT":"snutan","Ss-CAM1384-DRAFT":"sscopa","Td-FL_9056069_6-DRAFT":"tdacts","Td-KS_B6_1-DRAFT":"tdactn","Te-Pasquet1246-DRAFT":"telega","Tt-AUB21_1-DRAFT":"ttrian","Tz-DC_05_58_3A-DRAFT":"tzopol","Ud-Pasquet1171-DRAFT":"udigit","Vc-Pasquet1098-DRAFT":"vcuspi","Zd-Gigi-REFERENCE":"zdiplg","Zd-Momo-REFERENCE":"zdiplm","Zh-RIMHU001-REFERENCE":"zhuehu","Zn-PI615697-REFERENCE":"znicar","Zv-TIL01-REFERENCE":"zTIL01","Zv-TIL11-REFERENCE":"zTIL11","Zx-TIL18-REFERENCE":"zTIL18","Zx-TIL25-REFERENCE":"zTIL25"}

trips = ["tdactn","tdacts","zdiplg","zluxur","znicar","zTIL01","zTIL11","zTIL18","zTIL25","zB73v5","zhuehu","zdiplm"]

panand = {"aburma":["aburma@1","aburma@2"],"achine":["achine@1","achine@2"],"agerar":["agerar@1","agerar@2","agerar@3"],"avirgi":["avirgi@1"],"blagur":["blagur@1","blagur@2","blagur@3"],"ccitra":["ccitra@1","ccitra@2"],"crefra":["crefra@1"],"etrips":["etrips@1","etrips@2"],"hcompr":["hcompr@1","hcompr@2","hcompr@3","hcompr@4","hcompr@5","hcompr@6"],"hconto":["hconto@1","hconto@2","hconto@3","hconto@4"],"irugos":["irugos@1"],"ppanic":["ppanic@1"],"rrottb":["rrottb@1","rrottb@2"],"rtuber":["rtuber@1"],"sbicol":["sbicol@1"],"smicro":["smicro@1"],"snutan":["snutan@1","snutan@2"],"sscopa":["sscopa@1","sscopa@2"],"telega":["telega@1","telega@2"],"ttrian":["ttrian@1"],"udigit":["udigit@1","udigit@2","udigit@3"],"vcuscpi":["vcuspi@1","vcuspi@2"]}

### 0 Store subgenomes

subdict = {}
with open(subgenomes,'r') as window:
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


### 1 Store orthogroups in dict
odict = {}
with open(ortho_new,'r') as orthogroups:
    header = next(orthogroups)
    header = header.strip().split('\t')
    for l in orthogroups:
        l = l.strip().split('\t')
        ogid = l[0]
        odict[ogid] = {}
        refgenes = []
        for i,clust in enumerate(l[1:],1):
            sp = header[i]
            if 'tdac' not in sp:
                clust = clust.split(',')
                genes_clean = []
                for c in clust:
                     if "Sobic" in c:
                        genes_clean.append('.'.join(c.split('.')[0:2]))
                     else:
                         genes_clean.append(c.split('_')[0].strip())
                odict[ogid][sp] = genes_clean
# for old tdactn and tdacts assemblies
odict_old = {}
with open(ortho_old,'r') as orthogroups_old:
    header = next(orthogroups_old)
    header = header.strip().split('\t')
    for l in orthogroups_old:
        l = l.strip().split('\t')
        ogid = l[0]
        odict_old[ogid] = {}
        refgenes = []
        for i,clust in enumerate(l[1:],1):
            sp = header[i]
            if 'tdac' or 'zmays' in sp:
                clust = clust.split(',')
                genes_clean = []
                for c in clust:
                     if "Sobic" in c:
                        genes_clean.append('.'.join(c.split('.')[0:2]))
                     else:
                         genes_clean.append(c.split('_')[0].strip())
                odict_old[ogid][sp] = genes_clean

### 2 Store list collinear genes per zmv5 gene for each each panand species
# get block sizes
cb_dict = {}
with open(collinear_block,'r') as cb:
    for l in cb:
        l = l.strip().split('\t')
        cb_dict[l[0]] = [int(l[3]),l[4],l[5],l[6]]
# if block is large enough add gene 
c_dict = {}
with open(collinear,'r') as c:
    for l in c:
        l = l.strip().split('\t')
        aln_num = l[0]
        block_genes = cb_dict[aln_num][0]
        if block_genes >= 10:
            sp1 = l[2]
            sp2 = l[3]
            if "Zm_" == sp1[:3]:
                ref = sp1
                query = sp2
            else:
                ref = sp2
                query = sp1
            query_species = spdict[query.split('_')[0]]
            ref_gene = ref.split('_')[1]
            if "Sobic" in query:
                query_gene = '.'.join(query.split('_')[1].split('.')[0:2])
            else:
                query_gene = query.split('_')[1]
            if ref_gene not in c_dict:
                c_dict[ref_gene] = {}
            if query_species not in c_dict[ref_gene]:
                c_dict[ref_gene][query_species] = []
            c_dict[ref_gene][query_species].append(query_gene)

### 3 Store collinearity blocks from old mcscan run with pandand zmays annotation
with open(col_block_old,'r') as window:
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
            qchrom = '_'.join(l[0].split('.')[1:])
            qregion = '@'.join([qchrom] + l[1:3] + [l[9]])
            syn_info_line = region + "@" + qregion
            for idx,m in enumerate(region_chrom_list[:-1]):
                win_name = region_chrom + "@" + str(m) + "@" + str(region_chrom_list[idx+1])
                if lspecies not in windows[win_name]:
                    windows[win_name][lspecies] = []
                windows[win_name][lspecies].append(syn_info_line)

### 4 Store cactus based syntenic blocks per gene
zm = {}
with open(syn_gene,'r') as syn:
    for l in syn:
        l = l.strip().split('\t')
        gname = l[4].split('_')[-1]
        species = l[0]
        coords = [l[1],int(l[2]),int(l[3])]
        if gname not in zm:
            zm[gname] = {}
        if species not in zm[gname]:
            zm[gname][species] = [coords]
        coords_list = []
        for c in zm[gname][species]:
            c_chrom = c[0]
            c_start = c[1]
            c_end = c[2]
            # distance between current coords and existing coords
            if coords[0] == c_chrom:
                interval = interval_dist([coords[1],coords[2]],[c_start,c_end])
                if interval <= 25000:
                    new_start = min(coords[1],c_start)
                    new_end = max(coords[2],c_end)
                    coords_list.append([c_chrom,new_start,new_end])
                else:
                    coords_list.append(c)
                    if coords not in coords_list:
                        coords_list.append(coords)
            else:
                coords_list.append(c)
                if coords not in coords_list:
                    coords_list.append(coords)
        zm[gname][species] = coords_list

### 5 Store gene locations

gdict = {}
with open(gene_locations,'r') as genes:
    for l in genes:
        l = l.strip().split('\t')
        # ignore zlux because it has no annotation
        if l[5] in gffdict:
            species = gffdict[l[5]]
            geneid = l[4]
            if species not in gdict:
                gdict[species] = {}
            gdict[species][geneid] = [l[0],int(l[1]),int(l[2]),l[3]]


zmays_tandem_dup = []
with open(tandem_dups,'r') as tdup:
    for l in tdup:
        gene = l.strip().split('_')[0]
        zmays_tandem_dup.append(gene)

### 6 Loop over each zmv5 gene and show all synteny and orthogroup information for each species
rspecies = 'zB73v5'
counter = 0
for gene,v in gdict[rspecies].items():
    if gene not in zmays_tandem_dup:
        qspecies_list = list(spdict.values())
        qspecies_list.remove('zB73v5')
        chrom = v[0]
        start = v[1]
        end = v[2]
        strand = v[3]
        ## First retrieve all information for the ref gene
        counter += 1
        # orthogroup information
        orthodict = {}
        orthodict_old = {}
        for oid in odict:
            if rspecies in odict[oid]:
                if gene in odict[oid][rspecies]:
                    orthodict = odict[oid]
        for oid in odict_old:
            if rspecies in odict_old[oid]:
                if gene in odict_old[oid][rspecies]:
                    orthodict_old = odict_old[oid]

        # subgenome information per gene
        ref_subgenome = "NA"
        if rspecies in subdict:
            if chrom in subdict[rspecies]:
                ref_subgenome = subdict[rspecies][chrom]

        ## Loop over query species and retrieve relevant info for each
        for qspecies in qspecies_list:
            rorthologs,qorthologs = [],[]
            # Orthologs
            if qspecies not in ["tdactn","tdacts"]:
                if qspecies in orthodict:
                    qorthologs = orthodict[qspecies]
                    rorthologs = orthodict[rspecies]
            else:
                if qspecies in orthodict_old:
                    qorthologs = orthodict_old[qspecies]
                    rorthologs = orthodict_old[rspecies]


            # cactus synteny blocks per gene
            query_cac_syn = []
            if gene in zm:
                for species in zm[gene]:
                    if qspecies in species:
                        for item in zm[gene][species]:
                            item = item + [species]
                            query_cac_syn.append(item)

            # mcscan collinearity blocks based on private zm annotation
            nearest = chrom_list[min(range(len(chrom_list)), key = lambda i: abs(chrom_list[i]-start))]
            lookup = chrom + "@" + str(nearest) + "@" + str(nearest+window_size)
            query_mc_syn = "NA"
            # collinear genes using public zm annotation
            query_collinear_genes = []
            if gene in c_dict:
                if qspecies in c_dict[gene]:
                    query_collinear_genes = c_dict[gene][qspecies]
            #print("Mcscan collinearity:",query_collinear_genes)

            ## Finding comparable genes
            #1 Select the orthologs that are also collinear
            # These are now the candidate comparable genes
            co_ortho = []
            for qo in qorthologs:
                if qo in query_collinear_genes:
                    co_ortho.append(qo)
            #2 For each candidate, get location and check
            # if they are also aligned via cactus
            co_ortho_cac = {}
            for co in co_ortho:
                 # get location of gene
                 co_coords = gdict[qspecies][co]
                 for cac_coords in query_cac_syn:
                     subgenome = cac_coords[3]
                     chrom1 = cac_coords[0]
                     chrom2 = co_coords[0]
                     start1 = cac_coords[1]
                     end1 = cac_coords[2]
                     start2 = co_coords[1]
                     end2 = co_coords[2]
                     if chrom1 == chrom2:
                         distance = interval_dist([start1,end1],[start2,end2])
                         if distance == 0:
                             if subgenome not in co_ortho_cac:
                                 co_ortho_cac[subgenome] = []
                             co_ortho_cac[subgenome].append(co)
            #3 If there a multiple candidates within same
            # subgenome then pick one
            for sub,candidate in co_ortho_cac.items():
                if len(candidate) > 1:
                    sub_ortholog = random.choice(candidate)
                elif len(candidate) == 1:
                    sub_ortholog = candidate[0]
                else:
                    sub_ortholog = "NA"
                full_list = ','.join(candidate)
                print(gene,sub_ortholog,sub,qspecies,full_list)

