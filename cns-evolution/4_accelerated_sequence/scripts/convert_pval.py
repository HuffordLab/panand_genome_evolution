import sys


with open(sys.argv[1],'r') as results:
    for l in results:
        l = l.strip().split('\t')
        nlog10_pval = abs(float(l[3]))
        pval = 10**-nlog10_pval
        outline = l + [pval]
        print(*outline,sep='\t')
