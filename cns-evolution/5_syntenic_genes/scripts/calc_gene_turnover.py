import sys
import random
#
def interval_dist(r1, r2):
     x, y = sorted((r1, r2))
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0

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

# collinear genes and blocks based on zmv5 public
subgenomes = "../../1_subgenomes/data/panand_subgenomes.txt" 
# old get only tdactn and tdacts from here
ortho_old = "../data/Orthogroups_v7.txt"
# new but without tdactn and tdacts
ortho_new = "../data/Orthogroups.tsv"
#-collinear blocks based on panand annotation (not public zmv5)
zm_highqual_singletons = "../data/core_expr_chip_syn_genes.txt"
gene_locations = '../data/panand_gene_locations.bed'
tandem_dups = '../data/zB73v5_tandem_duplicates.txt'
 
### Target gene list ###
singletons = []
with open(zm_highqual_singletons,'r') as single:
    for l in single:
        l = l.strip()
        singletons.append(l)

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
    #print(k,v)
    if gene not in zmays_tandem_dup and gene in singletons:
        qspecies_list = list(spdict.values())
        qspecies_list.remove('zB73v5')
        chrom = v[0]
        start = v[1]
        end = v[2]
        strand = v[3]
        ## First retrieve all information for the ref gene
        counter += 1
        #if counter > 750:
        #    sys.exit()
        # orthogroup information
        orthodict = {}
        orthodict_old = {}
        for oid in odict:
            if rspecies in odict[oid]:
                if gene in odict[oid][rspecies]:
                    orthodict = odict[oid]
                    #if qspecies in odict[oid]:
                    #    print("Orthogroup:",odict[oid][qspecies])
        for oid in odict_old:
            if rspecies in odict_old[oid]:
                if gene in odict_old[oid][rspecies]:
                    orthodict_old = odict_old[oid]
                    #print("Orthogroup old:",odict_old[oid][qspecies])

        # subgenome information per gene
        ref_subgenome = "NA"
        if rspecies in subdict:
            if chrom in subdict[rspecies]:
                ref_subgenome = subdict[rspecies][chrom]
        #print("Subgenome: ",ref_subgenome)

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

            #3 If there a multiple candidates within same
            # subgenome then pick one
            if len(rorthologs) == 1:
                for q in qorthologs:
                    if q != '':
                        if q in gdict[qspecies]:
                            location = gdict[qspecies][q]
                            qchrom = location[0]
                            if qspecies in subdict:
                                if qchrom in subdict[qspecies]:
                                    subgenomes = list(set(subdict[qspecies][qchrom].values()))
                                    if 'NA' in subgenomes:
                                        subgenomes.remove('NA')
                                    if len(subgenomes)>0:
                                        if len(subgenomes)==1 or 'maize' in subgenomes[0]:
                                            subgenomes_out = qspecies + '@' + ','.join(subgenomes)
                                            print(gene,qspecies,rorthologs[0],q,subgenomes_out)
                                        else:
                                            # loop through blocks and find match
                                            matched_sub = []
                                            for block,sub in subdict[qspecies][qchrom].items():
                                                if sub != 'NA':
                                                    qstart = location[1]
                                                    qend = location[2]
                                                    bstart = block.split('@')[1]
                                                    bend = block.split('@')[2]
                                                    x = [int(qstart),int(qend)]
                                                    y = [int(bstart),int(bend)]
                                                    overlap = len(range(max(x[0], y[0]), min(x[-1], y[-1])+1))
                                                    if overlap > 0:
                                                        if sub not in matched_sub:
                                                            matched_sub.append(sub)
                                            if len(matched_sub) == 1:
                                                subgenomes_out = qspecies + '@' + matched_sub[0]
                                                print(gene,qspecies,rorthologs[0],q,subgenomes_out)
                                            else:
                                                subgenomes_out = qspecies + '@' + 'unknown'
                                                print(gene,qspecies,rorthologs[0],q,subgenomes_out)

                                    
