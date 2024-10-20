import sys
import copy
from collections import defaultdict
import itertools
import gzip
from Bio import AlignIO
from dijkstar import Graph, find_path

def range_dist(r1,r2):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((r1, r2))
     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the differnce of y[0] and x[1]
     #otherwise return 0
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0


if sys.argv[1].endswith('.gz'):
    maf_parser = AlignIO.parse(gzip.open(sys.argv[1],"rt"), "maf")
else:
    maf_parser = AlignIO.parse(sys.argv[1],"maf")

qspecies = sys.argv[2].split(',')
graph_dict = {}
for qsp in qspecies:
    graph_dict[qsp] = Graph()
# start with a dummy start point as the "previous_node_list"

# connect all nodes in first alignment to this with a distance of 1

# then replace previous_node_list with all current nodes


# ID: refchr@refstart@refend@qchr@qstart@qend

previous_alns = defaultdict(list)
init_nodes = {}
for multiple_alignment in maf_parser:
    if multiple_alignment.__len__() > 0:
        rec_count = 0
        q_alns = defaultdict(list)
        for record in multiple_alignment:
            sp = str(record.id).split(".")[0].split('@')[0]
            chrom = "".join(str(record.id).split(".")[1:])
            start = int(record.annotations["start"])
            size = int(record.annotations["size"])
            end = start + size
            refseq = str(record.seq)
            if int(record.annotations["strand"]) ==1:
                strand = "+"
            else:
                strand = "-"
            src_size = int(record.annotations["srcSize"])
            if strand == "-":
                start_fwd = src_size - end
                end_fwd = src_size - start
            else:
                start_fwd = start
                end_fwd = end
            if rec_count == 0:
                #assume first record is ref
                rec_count += 1
                ref_chrom = "chr" + "".join(str(record.id).split(".")[1:])
                ref_chrom_pos = [ref_chrom,start_fwd,end_fwd]
            else:
                alnid = ref_chrom_pos + [chrom,start_fwd,end_fwd]
                q_alns[sp].append(alnid)

        # Append to the graph
        for s,sp_q_alns in q_alns.items():
            if s not in previous_alns:
                init_node = '@'.join(map(str,ref_chrom_pos))
                if s not in init_nodes:
                    init_nodes[s] = init_node
                for qa in sp_q_alns:
                    qid = '@'.join((map(str, qa)))
                    # initial distance is equal to start pos
                    # this heuristic tries to push the path to begin at the start of the chr
                    init_dist = qa[1]
                    if s in graph_dict:
                        graph_dict[s].add_edge(init_node, qid, init_dist)
                previous_alns[s] = sp_q_alns
            else:
                for pa in previous_alns[s]:
                    prev_id = '@'.join(map(str,pa))
                    for qa in sp_q_alns:
                        qid = '@'.join(map(str,qa))
                        prev_chrom,prev_start,prev_end = pa[3:]
                        cur_chrom,cur_start,cur_end = qa[3:]
                        if prev_chrom == cur_chrom:
                            prev_dist = range_dist([prev_start,prev_end],[cur_start,cur_end])
                        else:
                            # between chrom dist arbitrarily set to 100m
                            prev_dist = 100000000
                        if s in graph_dict:
                            graph_dict[s].add_edge(prev_id, qid, prev_dist)
                # update the previous with the current
                previous_alns[s] = sp_q_alns
# Add a final path
# add arbitrary suffix to avoid overwriting previous node
last_node = '@'.join(map(str,ref_chrom_pos)) + "last"

for s,sp_q_alns in previous_alns.items():
    for pa in sp_q_alns:
        prev_id = '@'.join(map(str,pa))
        last_dist = 1
        if s in graph_dict:
            graph_dict[s].add_edge(prev_id, last_node, last_dist)

golden_dict = {}
for s,qgraph in graph_dict.items():
    golden_path = find_path(qgraph, init_nodes[s],last_node)
    golden_dict[s] = golden_path
#golden_nodes = golden_path[0]
#for i in golden_nodes:
    #print(i)
golden_lookup = defaultdict(dict)
for k,v in golden_dict.items():
    for node in v[0]:
        ref_id = '@'.join(node.split('@')[:3])
        q_id = '@'.join(node.split('@')[3:])
        golden_lookup[k][ref_id] = q_id 

if sys.argv[1].endswith('.gz'):
    maf_parser = AlignIO.parse(gzip.open(sys.argv[1],"rt"), "maf")
else:
    maf_parser = AlignIO.parse(sys.argv[1],"maf")

for multiple_alignment in maf_parser:
    if multiple_alignment.__len__() > 0:
        rec_count = 0
        q_alns = defaultdict(list)
        sp_list = []
        for record in multiple_alignment:
            sp = str(record.id).split(".")[0].split('@')[0]
            chrom = "".join(str(record.id).split(".")[1:])
            start = int(record.annotations["start"])
            size = int(record.annotations["size"])
            end = start + size
            refseq = str(record.seq)
            if int(record.annotations["strand"]) ==1:
                strand = "+"
            else:
                strand = "-"
            src_size = int(record.annotations["srcSize"])
            outline = ["s",record.id,start,size,strand,src_size,refseq]
            if strand == "-":
                start_fwd = src_size - end
                end_fwd = src_size - start
            else:
                start_fwd = start
                end_fwd = end
            if rec_count == 0:
                #assume first record is ref
                rec_count += 1
                ref_chrom = "chr" + "".join(str(record.id).split(".")[1:])
                ref_chrom_pos = [ref_chrom,start_fwd,end_fwd]
                print("")
                print("a")
                print(*outline,sep='\t')
            else:
                if sp not in sp_list:
                    q_id = '@'.join(map(str,[chrom,start_fwd,end_fwd]))
                    ref_id = '@'.join(map(str,ref_chrom_pos))
                    golden_aln = golden_lookup[sp][ref_id]
                    if q_id == golden_aln:
                        print(*outline,sep='\t')
                        sp_list.append(sp)
