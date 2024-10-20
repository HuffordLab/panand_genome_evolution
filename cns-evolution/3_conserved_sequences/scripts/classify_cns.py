import sys

def cns_class(inlist):
    classes = inlist[3].split(',')
    # UTR supercede introns
    if "five_prime_UTR" in classes:
        location = "five_prime_UTR"
    elif "three_prime_UTR" in classes:
        location = "three_prime_UTR"
    elif "intron" in classes:
        location = "intron"
    else:
        # its a distal CNS
        dist_to_gene = int(inlist[10])
        if dist_to_gene >= 0:
            # downstream
            if abs(dist_to_gene)>1000:
                #distal
                location = "downstream_distal"
            else:
                location = "downstream"
        else:
            #upstream
            if abs(dist_to_gene)>1000:
                #distal
                location = "upstream_distal"
            else:
                location = "upstream"
    return location


# Inline:
# chr1    35788   35795   chromosome,gene,intron,mRNA     chr1    34616   40204   Zm00001eb000010 . +       0

cdict = {}
with open(sys.argv[1],'r') as cns:
    for l in cns:
        l = l.strip().split('\t')
        cns_id = '@'.join(l[:3])
        if cns_id not in cdict:
            cdict[cns_id] = l
        #if l not in cdict[cns_id]:
        #    cdict[cns_id].append(l)

for k,v in cdict.items():
    location = cns_class(v)
    outline = k.split('@') + [location]
    print(*outline,sep='\t')



#1	554856	554931	panand_1.3	1	intron	554542	555072	-	75	1	550833	564766	Zm00001eb000170	.	-
#1	554856	554931	panand_1.3	1	intron	554628	555072	-	75	1	550833	564766	Zm00001eb000170	.	-
#1	554856	554931	panand_1.3	1	intron	554628	555072	-	75	1	550833	564766	Zm00001eb000170	.	-
#1	554856	554931	panand_1.3	1	intron	554628	555072	-	75	1	550833	564766	Zm00001eb000170	.	-
#1	555606	555627	panand_1.6	1	intron	555584	555676	-	21	1	550833	564766	Zm00001eb000170	.	-
#1	555606	555627	panand_1.6	1	intron	555584	555676	-	21	1	550833	564766	Zm00001eb000170	.	-
#1	555606	555627	panand_1.6	1	intron	555584	555676	-	21	1	550833	564766	Zm00001eb000170	.	-
#1	555606	555627	panand_1.6	1	intron	555584	555676	-	21	1	550833	564766	Zm00001eb000170	.	-
#1	679835	679882	panand_1.8	.	.	-1	-1	.	0	1	679956	680801	Zm00001eb000200	.	+
#1	684556	684635	panand_1.27	1	intron	684530	684853	-	79	1	679956	708016	Zm00001eb000190	.	-
