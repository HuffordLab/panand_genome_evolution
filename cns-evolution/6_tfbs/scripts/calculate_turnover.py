import sys
from collections import Counter
# Note that there will be aligned motifs without a denovo match
# Reason1: Even one mutation can prevent a denovo prediction
# Reason2: Perfect matching TFBS could be outside of 1000bp region

def pair_aln_denovo(qmotifs_aln,qmotifs_denovo):
    # denovo motifs with alignment
    denovo_match = []
    perfect_nomatch = []
    for qa in qmotifs_aln:
        qa_list = qa.split(';')
        qa_name = qa_list[0]
        mismatches = int(qa_list[2])
        gaps = int(qa_list[3])
        sum_mg = mismatches + gaps
        qa_positions = qa_list[4:]
        motif_has_match = False
        motif_within_100kb = False
        # some alignments are split and have multiple non-adjacent positions
        for qa_pos in qa_positions:
            qa_chrom = qa_pos.split('@')[0]
            qa_start = qa_pos.split('@')[1]
            qa_end = qa_pos.split('@')[2]
            for qd in qmotifs_denovo:
                qd_list = qd.split('@')
                qd_name = qd_list[0]
                qd_chrom = qd_list[1]
                qd_start = qd_list[2]
                qd_end = qd_list[3]
                

                # Check for match
                if qa_chrom == qd_chrom:
                    #print("Check potential match:",qa_name,qa_pos,qd)
                    # if coords overlap
                    x = [int(qa_start),int(qa_end)]
                    y = [int(qd_start),int(qd_end)]
                    overlap = len(range(max(x[0], y[0]), min(x[-1], y[-1])+1))
                    if overlap > 0 and qa_name == qd_name:
                        #print("Found good match:", qa_pos,qd)
                        if qd not in denovo_match:
                            denovo_match.append(qd)
                            motif_has_match = True
                    elif abs(max(x[0], y[0]) - min(x[-1], y[-1])) < 100000:
                        motif_within_100kb = True
        if not motif_has_match and sum_mg == 0 and motif_within_100kb:
            if qa not in perfect_nomatch:
                perfect_nomatch.append(qa)
    return denovo_match,perfect_nomatch

with open(sys.argv[1],'r') as tfbs:
    for l in tfbs:
        l = l.strip().split('\t')
        refgene = l[4]
        rmotifs = l[5].split(',')
        subgenome = l[6]
        # require a collinear gene with denovo predictions
        if l[10] != "NA":

            qmotifs_denovo = l[10].split(',')
            # if qgene is duplicated, only a single one is selected
            qgene_all = l[9].split(',')
            qmotifs_denovo = l[10].split(',')
            high_qual_motifs = []
            high_qual_motifs_denovo = []
            # if there are no alignments, then keep high_qual lists empty
            if l[7] != "NA":
                qmotifs_aln = l[7].split(',')
                denovo_matches,perfect_nomatches = pair_aln_denovo(qmotifs_aln,qmotifs_denovo)
                high_qual_motifs = []
                high_qual_motifs_denovo = []
                for dm in denovo_matches:
                    high_qual_motifs.append(dm.split('@')[0])
                    high_qual_motifs_denovo.append(dm.split('@')[0])
                for pn in perfect_nomatches:
                    high_qual_motifs.append(pn.split(';')[0])
                #print(denovo_matches,perfect_nomatches)
            
            rmotif_names = []
            for n in rmotifs:
                name = n.split('@')[0]
                rmotif_names.append(name)
            qmotif_names = []
            for n in qmotifs_denovo:
                name = n.split('@')[0]
                qmotif_names.append(name)
            
            # remove the ref motifs with a high-quality alignment
            r_wo_hq = list((Counter(rmotif_names) - Counter(high_qual_motifs)).elements())
            q_wo_hq = list((Counter(qmotif_names) - Counter(high_qual_motifs_denovo)).elements())
            # For each remaining ref motif match it with a denovo motif (excluding the high qual)
            rq_wo_hq = list((Counter(r_wo_hq) & Counter(q_wo_hq)).elements())
            r_w_hq = list((Counter(rmotif_names) & Counter(high_qual_motifs)).elements())
            r_wo_hq_lq = list((Counter(r_wo_hq) - Counter(q_wo_hq)).elements())
            q_wo_hq_lq = list((Counter(q_wo_hq) - Counter(r_wo_hq)).elements())
            # Number of ref motifs total
            r_tot = len(rmotif_names)
            # Number of ref motifs total with high-qual match (denovo or perfect distal)
            r_hq = len(r_w_hq)
            # Number of ref motifs with low-qual match (denovo only)
            r_lq = len(rq_wo_hq)
            # Number of ref motifs with no match
            r_orphan = len(r_wo_hq_lq)
            # Number of denovo motifs total
            q_tot = len(qmotif_names)
            # Number of denovo motifs with high-qual match
            q_hq = len(high_qual_motifs_denovo)
            # Number of denovo motifs with low-qual match, same as ref motifs with low-qual match
            q_lq = len(rq_wo_hq)
            # Number of denovo motifs with no match
            q_orphan = len(q_wo_hq_lq)

            ## PRINT RESULTS
            out_stats = [r_tot,r_hq,r_lq,r_orphan,q_tot,q_hq,q_lq,q_orphan]

            outline = l + out_stats
            print(*outline,sep='\t')

