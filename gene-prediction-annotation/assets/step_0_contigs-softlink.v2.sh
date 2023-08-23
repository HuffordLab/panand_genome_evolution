rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Andropogon-chinensis/1_hifiasm.v2/Andropogon-chinensis.bp.p_ctg.fasta.gz Andropogon-chinensis_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Andropogon-tenuilfolius/1_hifiasm.v2/Andropogon-tenuilfolius.bp.p_ctg.fasta.gz Andropogon-tenuilfolius_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Cymbopogon-citratus/1_hifiasm.v2/Cymbopogon-citratus.bp.p_ctg.fasta.gz Cymbopogon-citratus_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Ischaemum-rugosum/1_hifiasm.v2/Ischaemum-rugosum.bp.p_ctg.fasta.gz Ischaemum-rugosum_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Pogonatherum-paniceum/1_hifiasm.v2/Pogonatherum-paniceum.bp.p_ctg.fasta.gz Pogonatherum-paniceum_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Rhytacne-rottboelloides/1_hifiasm.v2/Rhytacne-rottboelloides.bp.p_ctg.fasta.gz Rhytacne-rottboelloides_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Rottboellia-tuberculosa/1_hifiasm.v2/Rottboellia-tuberculosa.bp.p_ctg.fasta.gz Rottboellia-tuberculosa_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Schizachyrium-microstachyum/1_hifiasm.v2/Schizachyrium-microstachyum.bp.p_ctg.fasta.gz Schizachyrium-microstachyum_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Thelopogon-elegans/1_hifiasm.v2/Thelopogon-elegans.bp.p_ctg.fasta.gz Thelopogon-elegans_v0.fasta.gz
rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Tripsacum-zopolitense/1_hifiasm.v2/Tripsacum-zopolitense.bp.p_ctg.fasta.gz Tripsacum-zopolitense_v0.fasta.gz
#rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Zea-diploperennis-momo/1_hifiasm.v2/Zea-diploperennis-momo.bp.p_ctg.fasta.gz Zea-diploperennis-momo_v0.fasta.gz
#rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Zea-mays-ssp-huehuetenangensis/1_hifiasm.v2/Zea-mays-ssp-huehuetenangensis.bp.p_ctg.fasta.gz Zea-mays-ssp-huehuetenangensis_v0.fasta.gz
#rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Zea-mays-ssp-maysXmexicana/1_hifiasm.v2/Zea-mays-ssp-maysXmexicana.bp.p_ctg.fasta.gz Zea-mays-ssp-maysXmexicana_v0.fasta.gz
#rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Zea-mays-ssp-maysXparviglumis/1_hifiasm.v2/Zea-mays-ssp-maysXparviglumis.bp.p_ctg.fasta.gz Zea-mays-ssp-maysXparviglumis_v0.fasta.gz
#rsync -avP /lss/research/mhufford-lab/arnstrm/PanAnd/Zea-nicaraguensis/1_hifiasm.v2/Zea-nicaraguensis.bp.p_ctg.fasta.gz Zea-nicaraguensis_v0.fasta.gz
for f in *.fasta.gz; do pigz -d -p8 $f; done
