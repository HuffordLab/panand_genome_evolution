# Genome Assembly


## Contig assembly

### Species with PacBio CLR data

CLR subreads were concatenated to a single fasta file, and error correction was run using `falcon` (config file `faclon_ec.cfg`). 

```bash
conda activate denovo_asm
cfg=faclon_ec.cfg
fc_run ${cfg}
 ```

All corrected reads were concatenated to a single file `genus-ec-reads.fa.gz` and `canu` was used for contig assembly:

```bash
conda activate denovo_asm
tdate=$(date '+%Y%m%d')
aname="genus"
cfg=schyzosachirium.cfg
ecreads="${genus}-ec-reads.fa.gz"
canu \
   -trim-assemble \
   -p $aname \
   -d canu-${tdate} \
   -s ${cfg} \
   -pacbio-corrected $ecreads
```

### Species with PacBio HiFi data

The species with HiFi data, the bam files were converted to `fasta` format and `hifiasm` was used to generate the assembly.

```bash
./runHiFiasm.v2.sh genus-hifi.fasta.gz Genus-species 
```

### Species with ONT data

The error correction was performed using `necat` and assembly was generated using `canu`

```bash

```

## Scaffolding using BioNano optical map


## Scaffolding/Pseudomolecule construction using HiC data (Tripsacum genomes only)


## Pseudomolecule construction (for _Zea_ species)




