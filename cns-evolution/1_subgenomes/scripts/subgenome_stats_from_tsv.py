import sys

def mergeIntervals(arr):

    # Sorting based on the increasing order
    # of the start intervals
    arr.sort(key=lambda x: x[0])

    # Stores index of last element
    # in output array (modified arr[])
    index = 0

    # Traverse all input Intervals starting from
    # second interval
    for i in range(1, len(arr)):

        # If this is not first Interval and overlaps
        # with the previous one, Merge previous and
        # current Intervals
        if (arr[index][1] >= arr[i][0]):
            arr[index][1] = max(arr[index][1], arr[i][1])
        else:
            index = index + 1
            arr[index] = arr[i]
    merged_arr = []
    for i in range(index+1):
        merged_arr.append(arr[i])
    return merged_arr

def subgenome_score(scaf_dict,path_scafs,genome_size):
    sub_bp = {}
    for sub,scafs in path_scafs.items():
        sub_regions_dict = {}
        unique_bp = 0
        for s in scafs:
            sub_regions = scaf_dict[s]
            for sr in sub_regions:
                refchrom = sr.split('@')[0]
                start = int(sr.split('@')[1])
                end = int(sr.split('@')[2])
                if refchrom not in sub_regions_dict:
                    sub_regions_dict[refchrom] = []
                else:
                    sub_regions_dict[refchrom] = sub_regions_dict[refchrom] + [[start,end]]
        for c,regions in sub_regions_dict.items():
            # merge overlapping regions
            unique_regions = mergeIntervals(regions)
            # sum ranges
            for merged_region in unique_regions:
                unique_bp += merged_region[1] - merged_region[0]


        sub_bp[sub] = unique_bp

    total_score = 0
    sub_scores = []
    #print("Path merged bp:", sub_bp)
    for s,bp in sub_bp.items():
        score = bp / genome_size
        sub_scores.append(score)
        total_score += score
    outscores = [total_score] + sub_scores
    return(outscores)






# Example usage
if __name__ == "__main__":
    #in_ranges = "test_input_Chr01_etrip.txt"
    in_ranges = sys.argv[1]
    #pvagin_chromsize = {"pvagin.Chr01":54347098,"pvagin.Chr02":51907944,"pvagin.Chr03":47877141,"pvagin.Chr04":45857669,"pvagin.Chr05":57965693,"pvagin.Chr06":42929029,"pvagin.Chr07":42926309,"pvagin.Chr08":30709806,"pvagin.Chr09":44791682,"pvagin.Chr10":44532209} 
    gsizes = {"pvagin":651047655,"aburma":1946409527,"achine":1335640055,"agerar":4084992212,"avirgi":905553168,"blagur":3243580366,"ccitra":2391180077,"crefra":807966905,"etrips":4764026201,"hcompr":6136823671,"hconto":2089701397,"irugos":704846194,"ppanic":799622294,"rrottb":1446746720,"rtuber":1572892435,"sbicol":708863705,"smicro":921361800,"snutan":4374693094,"sscopa":2543548284,"tdactn":3694024913,"tdacts":2976690102,"telega":3027872127,"ttrian":1034932660,"udigit":4703656903,"vcuspi":1191887642,"zB73v5":2182075994,"zTIL01":2598056910,"zTIL11":2380575982,"zTIL18":2470011168,"zTIL25":2100365973,"zdiplg":2690685548,"zdiplm":2153424507,"zhuehu":2244829945,"zluxur":3064462296,"znicar":2490707474}
    # 1) for each q species collect unique bp per scaffold with any synteny info
    # 2) for each q species subgenome collect unique pvag bp per chromosome
    # 3) for each q species collect unique bp per scaffold assigned to any subgenome
    # 4) Calculate % of q species unassigned, % of q species syntenic unassigned, average % pvag covered per subgenome, and list of % covered per subgenome

    # how much of query aligned
    self_dict = {}
    # how much of ref aligned
    ref_dict = {}

    with open(in_ranges,'r') as infile:
        for l in infile:
            l = l.split('\t')
            species = l[0].split('.')[0]
            qchrom = l[0].split('.')[1]
            qstart = int(l[1])
            qend = int(l[2])
            subgenome = l[12]
            
            rchrom = l[3]
            rstart = int(l[4])
            rend = int(l[5])

            if species not in self_dict:
                self_dict[species] = {}
            if species not in ref_dict:
                ref_dict[species] = {}
            if subgenome not in self_dict[species]:
                self_dict[species][subgenome] = {}
            if qchrom not in self_dict[species][subgenome]:
                self_dict[species][subgenome][qchrom] = []
            if subgenome not in ref_dict[species]:
                ref_dict[species][subgenome] = {}
            if rchrom not in ref_dict[species][subgenome]:
                ref_dict[species][subgenome][rchrom] = []

            if species == "tdactn":
                for m in ["maize1a","maize1b","maize2a","maize2b"]:
                    if m not in ref_dict[species]:
                        ref_dict[species][m] = {}
                        ref_dict[species][m][rchrom] = []
                    elif rchrom not in ref_dict[species][m]:
                        ref_dict[species][m][rchrom] = []
                    if m not in self_dict[species]:
                        self_dict[species][m] = {}
                        self_dict[species][m][qchrom] = []
                    elif qchrom not in self_dict[species][m]:
                        self_dict[species][m][qchrom] = []

                if subgenome == "maize1":
                    ref_dict[species]["maize1a"][rchrom].append([rstart,rend])
                    self_dict[species]["maize1b"][qchrom].append([qstart,qend])
                elif subgenome == "maize2":
                    ref_dict[species]["maize2a"][rchrom].append([rstart,rend])
                    self_dict[species]["maize2b"][qchrom].append([qstart,qend])
                else:
                    ref_dict[species][subgenome][rchrom].append([rstart,rend])
            else:
                ref_dict[species][subgenome][rchrom].append([rstart,rend])
                self_dict[species][subgenome][qchrom].append([qstart,qend])

    for sp,sub in ref_dict.items():
        #sp_size = gsizes[sp]
        ref_size = gsizes["pvagin"]
        for scaf,regions in sub.items():
            unique_bp = 0
            for k,v in regions.items(): 
                if len(v)>0:
                    unique_regions = mergeIntervals(v)
                    # sum ranges
                    for merged_region in unique_regions:
                        unique_bp += merged_region[1] - merged_region[0]
            print("ref",sp,scaf,unique_bp/ref_size,unique_bp)

    for sp,sub in self_dict.items():
        sp_size = gsizes[sp]
        sp_unique_bp = 0
        sp_size_syntenic = 0
        for scaf,regions in sub.items():
            unique_bp = 0
            for k,v in regions.items(): 
                if len(v)>0:
                    unique_regions = mergeIntervals(v)
                    # sum ranges
                    for merged_region in unique_regions:
                        unique_bp += merged_region[1] - merged_region[0]
                        if scaf != "NA":
                            sp_unique_bp += merged_region[1] - merged_region[0]
            sp_size_syntenic += unique_bp
            print("self",sp,scaf,unique_bp/sp_size,unique_bp)
        print("self",sp,"all",sp_unique_bp/sp_size,sp_unique_bp)
        print("self",sp,"all_syn",sp_unique_bp/sp_size_syntenic,sp_unique_bp)






