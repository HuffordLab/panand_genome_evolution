import sys
#import networkx as nx

# To do 09062023
# the one scaf per path restriction can cause issues if 10% of a scaf can be assigned to one path but 90% to another path and the 10% path comes first
# A solution may be to do two loops, freely assigning scafs first and then restricting to one best matching scaf per path






#def compute_shortest_path(start_positions, end_positions, graph):
#    shortest_paths = []
#    for start in start_positions:
#        for end in end_positions:
#            if start != end:
#                # Compute the shortest path with minimum gaps
#                try:
#                    shortest_path = nx.bidirectional_dijkstra(graph, start, end)
#                    shortest_path = shortest_path[1]  # Extract the path from the result
#                    # Append the shortest path to the result
#                    shortest_paths.append(shortest_path)
#                except:
#                    shortest_paths = [0]
#                    return shortest_paths
#    return shortest_paths

def range_dist(r1,r2):
    x, y = sorted((r1, r2),key=lambda r: r.start)
    # check they dont overlap
    #if x[0] <= x[-1] < y[0] and all( y[0] <= y[-1] for y in (r1,r2)):
    # return negative values if overlap
    return y[0] - x[-1]

def best_scoring_range(ranges,range_names,cur_query_pos):
    longest = None
    longest_contiguous = None
    longest_name = None
    longest_contiguous_name = None
    max_length = 0
    max_length_contiguous = 0
    max_score = 0
    max_score_contiguous = 0
    score = 0
    for idx,l in enumerate(ranges):
        # length of range
        range_length = len(l)
        if l != range(0,0):
            name = range_names[idx]
            score = int(align_dict[name][9])
        if cur_query_pos and l != range(0,0):
            pos_key = chrom + '_' + str(l[0]) + '_' + str(l[-1]+1) + '_' + name
            range_query_pos = query_dict[pos_key]
            query_chrom = range_query_pos.split('@')[0]
            cur_query_chrom = cur_query_pos.split('@')[0]
            distance = range_dist(range(int( range_query_pos.split('@')[1]),int(range_query_pos.split('@')[2])),range(int(cur_query_pos.split('@')[1]),int(cur_query_pos.split('@')[2])))
            # avoid large overlaps as they may be poor alignments
            if query_chrom == cur_query_chrom:
                if score > max_score_contiguous:
                    max_score_contiguous = score
                    longest_contiguous = l
                    longest_contiguous_name = name
        if score > max_score:
            max_score = score
            longest = l
            longest_name = name
    if max_score_contiguous > 10:
        return longest_contiguous,longest_contiguous_name
    return longest,longest_name

def longest_range(ranges,range_names,cur_query_pos):
    longest = None
    longest_contiguous = None
    longest_name = None
    longest_contiguous_name = None
    max_length = 0
    max_length_contiguous = 0
    for idx,l in enumerate(ranges):
        # length of range
        range_length = len(l)
        if l != range(0,0):
            name = range_names[idx]
        if cur_query_pos and l != range(0,0):
            pos_key = chrom + '_' + str(l[0]) + '_' + str(l[-1]+1) + '_' + name
            range_query_pos = query_dict[pos_key]
            query_chrom = range_query_pos.split('@')[0]
            cur_query_chrom = cur_query_pos.split('@')[0]
            distance = range_dist(range(int( range_query_pos.split('@')[1]),int(range_query_pos.split('@')[2])),range(int(cur_query_pos.split('@')[1]),int(cur_query_pos.split('@')[2])))
            # avoid large overlaps as they may be poor alignments
            if query_chrom == cur_query_chrom:
                if range_length > max_length_contiguous:
                    max_length_contiguous = range_length
                    longest_contiguous = l
                    longest_contiguous_name = name
        if range_length > max_length:
            max_length = range_length
            longest = l
            longest_name = name
    if max_length_contiguous > 500000:
        return longest_contiguous,longest_contiguous_name
    return longest,longest_name




def select_ranges(ranges,range_names,chrom,path):
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    sorted_names = [x for _, x in sorted(zip(ranges, range_names), key=lambda pair: pair[0][0])]

    clusters = []  # List to store the final cluster
    # 1) in first 1Mb of sequence find biggest block
    # 2) Where this block ends, look in 1Mb up/down region for next biggest block
    # 3) Keep stringing together the biggest blocks until n sets exist, where n=ploidy
    start_pos = sorted_ranges[0][0]
    end_pos = sorted_ranges[-1][-1]
    end_reached = False

    # start at the beginning
    cur_pos = None
    #print("Starting greedy search")
    #print("Input", sorted_ranges)
    syn_path = []
    syn_path_names = []
    query_path = []
    cur_query_pos = None
    while not end_reached:
        # look 1Mb left and right of start_pos and get longest block
        local_ranges = []
        local_ranges_names = []
        if not cur_pos:
            cur_pos = sorted_ranges[0][0]
        else:
            #print(syn_path)
            cur_pos = syn_path[-1][-1]
        #print("Cur_pos: ", cur_pos)
        extend_bp = 1000000
        while len(local_ranges) == 0:
            #print("Search range is :", extend_bp /1000000, "Mb")
            #print("Local ranges",len(local_ranges),local_ranges)
            #print("Local range names",len(local_ranges_names),local_ranges_names)
            for idx,s in enumerate(sorted_ranges):
                range_start = s[0]
                range_end = s[-1]

                # Check if a query scaf has already been assigned to an alt path
                # each scaf can only belong to one subgenome and one path
                on_alt_path = False
                s_key = chrom + '_' + str(range_start) + '_' + str(range_end+1) + '_' + sorted_names[idx]
                query_scaf = query_dict[s_key].split('@')[0]
                #for k,v in path_scafs.items():
                #    if k!=path:
                #        if query_scaf in v:
                #            on_alt_path = True

                # does new range start within search range of cur_pos
                if not on_alt_path and range_start >= (cur_pos - 1000000) and range_start <= (cur_pos + extend_bp) and range_end > cur_pos:
                #if range_start >= (cur_pos - 1000000) and range_start <= (cur_pos + extend_bp) and range_end > cur_pos:
                    if len(syn_path)>0:
                        cur_range = syn_path[-1]
                        overlap_frac = len(range(max(cur_range[0], s[0]), min(cur_range[-1], s[-1])+1)) / len(cur_range)
                        if overlap_frac < 0.1:
                            local_ranges.append(s)
                            local_ranges_names.append(sorted_names[idx])
                    else:
                        local_ranges.append(s)
                        local_ranges_names.append(sorted_names[idx])
            extend_bp += 500000
            if cur_pos + extend_bp > end_pos:
                end_reached = True
                # append dummy to break while loop
                local_ranges.append(range(0,0))
        if len(query_path) > 0:
            cur_query_pos = query_path[-1]
        #longest_local,ll_name = longest_range(local_ranges,local_ranges_names,cur_query_pos)
        longest_local,ll_name = best_scoring_range(local_ranges,local_ranges_names,cur_query_pos)

        # while be false if range(0,0) dummy was used
        if longest_local:
            syn_path.append(longest_local)
            pos_key = chrom + '_' + str(longest_local[0]) + '_' + str(longest_local[-1]+1) + '_' + ll_name
            query_scaf = query_dict[pos_key].split('@')[0]
            # increment how many ref bp assigned to query scaf per path
            if query_scaf in path_scafs[path]:
                path_scafs[path][query_scaf] += ((longest_local[-1]+1) - longest_local[0])
            else:
                path_scafs[path][query_scaf] = ((longest_local[-1]+1) - longest_local[0])
            #path_scafs[path].append(query_scaf)
            query_path.append(query_dict[pos_key])
            syn_path_names.append(ll_name)
    return syn_path,query_path,syn_path_names

def select_ranges_round2(ranges,range_names,chrom,path):
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    sorted_names = [x for _, x in sorted(zip(ranges, range_names), key=lambda pair: pair[0][0])]

    clusters = []  # List to store the final cluster
    # 1) in first 1Mb of sequence find biggest block
    # 2) Where this block ends, look in 1Mb up/down region for next biggest block
    # 3) Keep stringing together the biggest blocks until n sets exist, where n=ploidy
    start_pos = sorted_ranges[0][0]
    end_pos = sorted_ranges[-1][-1]
    end_reached = False

    # start at the beginning
    cur_pos = None
    #print("Starting greedy search")
    #print("Input", sorted_ranges)
    syn_path = []
    syn_path_names = []
    query_path = []
    cur_query_pos = None
    while not end_reached:
        # look 1Mb left and right of start_pos and get longest block
        local_ranges = []
        local_ranges_names = []
        if not cur_pos:
            cur_pos = sorted_ranges[0][0]
        else:
            #print(syn_path)
            cur_pos = syn_path[-1][-1]
        #print("Cur_pos: ", cur_pos)
        extend_bp = 1000000
        while len(local_ranges) == 0:
            #print("Search range is :", extend_bp /1000000, "Mb")
            #print("Local ranges",len(local_ranges),local_ranges)
            #print("Local range names",len(local_ranges_names),local_ranges_names)
            for idx,s in enumerate(sorted_ranges):
                range_start = s[0]
                range_end = s[-1]

                # Check if a query scaf has already been assigned to an alt path
                # each scaf can only belong to one subgenome and one path
                on_alt_path = False
                s_key = chrom + '_' + str(range_start) + '_' + str(range_end+1) + '_' + sorted_names[idx]
                query_scaf = query_dict[s_key].split('@')[0]
                max_path_bp = 0
                for k,v in path_scafs.items():
                    if k!=path:
                        for a in v:
                            if query_scaf == a:
                                on_alt_path = True


                # does new range start within search range of cur_pos
                if not on_alt_path and range_start >= (cur_pos - 1000000) and range_start <= (cur_pos + extend_bp) and range_end > cur_pos:
                    if len(syn_path)>0:
                        cur_range = syn_path[-1]
                        overlap_frac = len(range(max(cur_range[0], s[0]), min(cur_range[-1], s[-1])+1)) / len(cur_range)
                        if overlap_frac < 0.1:
                            local_ranges.append(s)
                            local_ranges_names.append(sorted_names[idx])
                    else:
                        local_ranges.append(s)
                        local_ranges_names.append(sorted_names[idx])
            extend_bp += 500000
            if cur_pos + extend_bp > end_pos:
                end_reached = True
                # append dummy to break while loop
                local_ranges.append(range(0,0))
        if len(query_path) > 0:
            cur_query_pos = query_path[-1]
        #longest_local,ll_name = longest_range(local_ranges,local_ranges_names,cur_query_pos)
        longest_local,ll_name = best_scoring_range(local_ranges,local_ranges_names,cur_query_pos)

        # while be false if range(0,0) dummy was used
        if longest_local:
            syn_path.append(longest_local)
            pos_key = chrom + '_' + str(longest_local[0]) + '_' + str(longest_local[-1]+1) + '_' + ll_name
            query_path.append(query_dict[pos_key])
            syn_path_names.append(ll_name)
    return syn_path,query_path,syn_path_names

def unique_length(chr_ranges):
    '''
    Take a list of chr@start@end strings representing alignments,
    and return the total length of non overlapping alignments.
    '''
    chr_dict = {}
    chr_cov = {}
    for a in chr_ranges:
        chrom = a.split('@')[0]
        #a_r = range(int(a.split('@')[1]), int(a.split('@')[2])+1)
        a_r = [int(a.split('@')[1]), int(a.split('@')[2])+1]
        if chrom in chr_dict:
            chr_dict[chrom].append(a_r)
        else:
            chr_dict[chrom] = [a_r]
    genome_len = 0
    for k,v in chr_dict.items():
        merged = unique_range(v)
        chrom_len = 0
        for i in merged:
            chrom_len += len(range(i[0],i[1]+1))
        chr_cov[k] = chrom_len
        genome_len += chrom_len
    return genome_len,chr_cov

def unique_range(intervals):
    '''
    Take a list of ranges and output their
    non overlapping total length
    '''
    # Sort the array on the basis of start values of intervals.
    intervals.sort()
    stack = []
    # insert first interval into stack
    stack.append(intervals[0])
    for i in intervals[1:]:
        # Check for overlapping interval,
        # if interval overlap
        if stack[-1][0] <= i[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], i[-1])
        else:
            stack.append(i)
    #total_len = 0
    #for i in stack:
    #    total_len += len(range(i[0],i[1]+1))
    return stack

def greedy_split(scaf_dict,ploidy):
    '''
    Split scaffolds into n groups where n is the
    expected ploidy. A greedy algorithm can
    maximize the total unique ref genome 
    coverage across the n groups.
    '''

    # Ideas:
    # Initialize n groups
    # Iterate over scafs
    # Just always dump each scaf into the group with least overlap
    # if multiple then choose randomly

    # The these groups can be passed to the other greedy algo that finds a path for each group

    subgenomes = {}
    for i in range(ploidy):
        ploid = i+1
        subgenomes[ploid] = []

    for scaf,aln in scaf_dict.items():
        min_overlap = 1
        min_sub = None
        #print("Assigning",scaf)
        for sub,sub_scaf in subgenomes.items():
            overlap_frac = interval_overlap(scaf,sub_scaf,scaf_dict)
            #print("Subgenome, overlap_fraction",sub,overlap_frac)
            if overlap_frac < min_overlap:
                min_overlap = overlap_frac
                min_sub = sub
                #print("Overlap is less than minimum")
        if min_sub:
            subgenomes[min_sub].append(scaf)
    return(subgenomes)
    
def interval_overlap(scaf,sub_scaf,scaf_dict):
    # get both lists of chr_pos
    # how many bases of the total bases in scaf alignments is covered in sub_scaf alignments
    scaf_aln = scaf_dict[scaf]
    scaf_aln_bp = 0
    sa_overlap_bp = 0
    sub_scaf_aln = []
    for s in sub_scaf:
        for a in scaf_dict[s]:
            sub_scaf_aln.append(a)
    for sa in scaf_aln:
        chrom = sa.split('@')[0]
        start = int(sa.split('@')[1])
        end = int(sa.split('@')[2])+1
        a_r = [start, end]
        scaf_aln_bp += (end - start)
        partial_overlaps = []
        for ssa in sub_scaf_aln:
            ssa_chrom = ssa.split('@')[0]
            ssa_start = int(ssa.split('@')[1])
            ssa_end = int(ssa.split('@')[2])+1
            ssa_r = [ssa_start,ssa_end]
            if chrom == ssa_chrom:
                overlap = range(max(start, ssa_start), min(end, ssa_end))
                if len(overlap)>0:
                    partial_overlaps.append([overlap[0],overlap[-1]+1])
        if len(partial_overlaps)>0:
            merged = unique_range(partial_overlaps)
            for i in merged:
                sa_overlap_bp += len(range(i[0],i[1]+1))
    overlap_fraction = 0
    if sa_overlap_bp>0:
        overlap_fraction = sa_overlap_bp / scaf_aln_bp
    return overlap_fraction

# Example usage
if __name__ == "__main__":
    #in_ranges = "test_input_Chr01_etrip.txt"
    in_ranges = sys.argv[1]
    pvagin_chromsize = {"pvagin.Chr01":54347098,"pvagin.Chr02":51907944,"pvagin.Chr03":47877141,"pvagin.Chr04":45857669,"pvagin.Chr05":57965693,"pvagin.Chr06":42929029,"pvagin.Chr07":42926309,"pvagin.Chr08":30709806,"pvagin.Chr09":44791682,"pvagin.Chr10":44532209} 
    #path_scafs = {}
    #for p in range(int(sys.argv[2])):
    #    path_i = p + 1
    #    path_scafs[path_i] = {}
    with open(in_ranges) as f:
        for line in f:
            pass
        last_line = line.strip().split('\t')
    #etrips.scaf_11  101018164       101204671       pvagin.Chr01    1868095 2157007 Alignment13464  374.08.2e-18  8       +       +
    
    #############
    ## ROUND 1 ##
    #############

    pvagin_len = sum((pvagin_chromsize.values()))
    path_scafs = {}
    expected_paths = int(sys.argv[2])
    #for p in range(expected_paths):
    #    path_i = p + 1
    #    path_scafs[path_i] = {}
    align_dict = {}
    scaf_dict = {}
    with open(in_ranges,'r') as granges:
        for g in granges:
            g = g.strip().split('\t')
            if int(g[9]) > 10:
                #edges.append((int(g[4]),int(g[5])))
                #ranges.append(range(int(g[4]),int(g[5])))
                #range_names.append(g[6])
                ref_pos = g[3] + "@" + g[4] + "@" + g[5]
                if g[0] in scaf_dict:
                    scaf_dict[g[0]].append(ref_pos)
                else:
                    scaf_dict[g[0]] = [ref_pos]
            pos_str = g[0] + "@" + g[1] + "@" + g[2] + "@" + g[3] + "@" + g[4] + "@" + g[5]
            align_dict[pos_str] = g

    path_scafs = greedy_split(scaf_dict,expected_paths)
    #bag1 = result[1]
    #bag2 = result[2]
    #path_scafs[1] = bag1
    #path_scafs[2] = bag2

    #############
    ## ROUND 2 ##
    #############
    with open(in_ranges,'r') as granges:
        expected_paths = int(sys.argv[2])
        ranges = []
        range_names = []
        # holds target genome information
        align_dict = {}
        query_dict = {}
        edges = []
        chrom = None
        assigned_aln = []
        for g in granges:
            g = g.strip().split('\t')
            if not chrom:
                chrom = g[3]
            # switch to next chrom
            if g[3] != chrom or g == last_line:
                # find the paths
                while expected_paths > 0:
                    #print("Path scafs:",path_scafs)
                    selected_ranges,query_path,path_names = select_ranges_round2(ranges,range_names,chrom,expected_paths)
                    #selected_names = []
                    sum_ranges = 0
                    for s in selected_ranges:
                        sum_ranges += len(s)
                    chrom_size = pvagin_chromsize[chrom]
                    frac_cov = sum_ranges / chrom_size
                    #print(chrom,expected_paths,frac_cov,len(ranges),selected_ranges,query_path,path_names)
                    #print("Path names:", path_names)
                    #print("Range names:", range_names)
                    pop_indices = []
                    for idx,n in enumerate(range_names):
                        #print("Checking if ",n," in pathnames")
                        if n in path_names:
                            #print("Selected:",n)
                            #print("Removing", ranges[idx],range_names[idx])
                            pop_indices.append(idx)
                    for idx in sorted(pop_indices,reverse=True):
                        del ranges[idx]
                        del range_names[idx]

                    for idx,alignment in enumerate(query_path):
                        outline = align_dict[path_names[idx]] + [expected_paths, frac_cov, len(ranges)]
                        print(*outline,sep='\t')
                        assigned_aln.append(path_names[idx])
                    expected_paths -= 1
                # update chrom
                chrom = g[3]
                # reset ranges
                ranges = []
                range_names = []
                expected_paths = int(sys.argv[2])
            # require 15-gene block minimum
            if int(g[9]) > 10:
                edges.append((int(g[4]),int(g[5])))
                ranges.append(range(int(g[4]),int(g[5])))
                range_names.append(g[6])
                pos_str = g[3] + "_" + g[4] + "_" + g[5] + "_" +g[6]
                query_chrpos = g[0] + "@" + g[1] + "@" + g[2]
                query_dict[pos_str] = query_chrpos
            align_dict[g[6]] = g





for k,v in align_dict.items():
    if k not in assigned_aln:
        outline = v + ["NA","NA","NA"]
        print(*outline,sep='\t')



    #paths = [0]
    #while paths == [0]:
        # djikstra might work, if you just extend each edge by a fixed amount until you can find a path trhough the graph
        # Define start and end positions
        #sorted_ranges = sorted(ranges, key=lambda x: x[0])
        #start_pos = [sorted_ranges[0][0]]
        #end_pos = [sorted_ranges[-1][-1]]
        # Compute the shortest paths
        #graph = nx.Graph()
        #graph.add_edges_from(edges)
        #paths = compute_shortest_path(start_pos, end_pos, graph)
        #print(paths)
        #print(edges)
        #for idx,i in enumerate(edges):
        #    i = list(i)
        #    i[0] = i[0] - 100000
        #    i[1] = i[1] + 100000
        #    if i[0] <0:
        #        i[0] = 0
        #    if i[1] > end_pos[0]:
        #        i[1] = end_pos[0]
        #    edges[idx] = tuple(i)
    #print(paths)
