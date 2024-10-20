bedtools intersect -wao -a cns_cov0.25_len8_merge_maskref_nonCDS_min5bp_closestGene_blastxfilt.bed  -b ../../../zB73v5_annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.with_intron.gff3 |cut -f1-3,17|bedtools merge -i - -d -1 -o distinct -c 4 > cns_cov0.25_len8_merge_maskref_nonCDS_min5bp_closestGene_blastxfilt_classes.bed
bedtools closest -D b -t first -a cns_cov0.25_len8_merge_maskref_nonCDS_min5bp_closestGene_blastxfilt_classes.bed -b <(sed 's/^/chr/' ../../phastcons/old/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.genes.bed| sed 's/^chrscaf/scaf/' | sort -k1,1 -k2,2n) > cns_cov0.25_len8_merge_maskref_nonCDS_min5bp_closestGene_blastxfilt_classes_genedist.bed
python classify_cns.py cns_cov0.25_len8_merge_maskref_nonCDS_min5bp_closestGene_blastxfilt_classes_genedist.bed > cns_cov0.25_len8_merge_maskref_nonCDS_min5bp_closestGene_blastxfilt_single_class.bed
