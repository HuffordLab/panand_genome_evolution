
name="$1"
echo "$name"| sed 's/chr//'| sed 's/_/\t/g'| sed 's/scaf\t/scaf_/' > bed/$name
chrom=`echo "$name"|sed 's/scaf_/scaf@/'| cut -f1 -d'_' | sed 's/scaf@/scaf_/'`
mafsInRegion bed/${name} ${name}.maf maf/${chrom}_phastcons.maf
python cns_maf2mfa.py ${name}.maf aburma,achine,agerar,avirgi,blagur,ccitra,crefra,etrips,hconto,hcompr,irugos,ppanic,rrottb,rtuber,sbicol,smicro,snutan,sscopa,telega,ttrian,udigit,vcuspi,tdactn,tdacts,zB73v5,zdiplg,zdiplm,zhuehu,zluxur,znicar,zTIL01,zTIL11,zTIL18,zTIL25,tzopol,sspont,pvagin >> all_species_cns.txt
rm -f bed/$name ${name}.maf
