import sys
sdict = {}
with open(sys.argv[1],'r') as counts:
    for l in counts:
        l = l.strip().split(' ')
        gene = l[0]
        species = l[1]
        subgenome = l[4]
        if species not in sdict:
            sdict[species] = {}
        if subgenome not in sdict[species]:
            sdict[species][subgenome] = set()
        sdict[species][subgenome].add(gene)

for species,subgenome in sdict.items():
    for sub,genes1 in subgenome.items():
        if 'maize' in sub or 'unknown' in sub:
            for sub,genes2 in subgenome.items():
                # add to all other sets
                sdict[species][sub] = genes1.union(genes2)

for species,subgenome in sdict.items():
    for sub,genes in subgenome.items():
        print(species,sub,len(genes))
