# Gene Prediction

Scripts used below can be found in the `assets` folder.

RNAseq data generated for various tissues were used for evidence 


Data was inspected to ensure the sequence quality was satisfactory

```bash
mkdir Genus-species
cd Genus-species
for fq in *.fastq.gz; do
  fastqc --threads ${SLURM_CPUS_ON_NODE} $fq
done
multiqc .
```


## Scripts used for annotation:



- The `step_1_make-star-slurm.sh` was used to generate slurm script for mapping reads using STAR mapping program. The script was submitted using `sbatch Genus-species_0.sub`. Once the run was complete, the directory was organized using `step_1_clean-dir.sh`. It was run within the `Genus-species` folder, to move the STAR Db, misc files, bam files and slurm job files to respective directory. Mapping stats were collated using `multiqc`

- Genome summary stats, including BUSCO runs, separate primary and alternative haplotypes from the genome `fasta` files based on the scaffold names was carried out using `step_2_genome-summary-stats-full.sh` script.

- The `step_3_prepare-star-primary.sh` was used to generate slurm script for mapping reads using STAR mapping program (for only primary scaffolds). Once the run was complete, the directory was organized using `step_3_clean-star-run.sh`. It was run within the `Genus-species` folder, to move the STAR Db, misc files, bam files and slurm job files to respective directory. Mapping stats were collated using `multiqc`

- The merged BAM file generated in the previous step was then used for genome guided transcript assembly. The `step_4_make-transcript-assemblies.sh` was used to run various assemblers and BRAKER. It was launched as

```bash
step_4_make-transcript-assemblies.sh Genus-species
```

- Once the Trinity run was complete, the `step_4_gmap-trinity-transcripts.sh` was used to map the Trinity genome guided transcriptomes back to the genome assembly.

```bash
step_4_gmap-trinity-transcripts.sh Genus-species
```

- The assembled transcriptomes were then finalized using `step_4_make-mikado-container.sh` that runs the mikado pipeline to pick the best transcript for each locus.

```bash
step_4_make-mikado-container.sh Genus-species
```

- The BRAKER predictions, mikado transcripts and homology predictions (rice) were them combined using the `step_4_run-gemoma.sh` script.

```bash
step_4_run-gemoma.sh Genus-species
```

- The annotations were then finalized using `step_6_finalize-gff.sh` scripts - which renames the gene/mRNA ids, adds the exon features, creates CDS and peptide fasta files for subsequent analyses. It also calculates detailed annotation statistics.

```bash
step_6_finalize-gff.sh Genus-species
```


## Quality control

- Counts for transcript assemblies and for evidence based gene models were computed using `step_4_count-trascripts-and-genes.sh` 

- BUSCO for genome assemblies (`step_2_genome-summary-stats-full.sh`) and predictions (`step_6_run-busco-on-predictions.sh`) were computed and compared (`step_6_get-busco-stats.sh`)

- Phylostratiography (using `phylostratR`) was performed using `step_6_run-phylostratr-on-prot.sh` and `step_6_run-phylostratr-on-predictions.sh` scripts

- TESorter was used to detect any TE containing genes in the predictions (`step_6_screen-te-in-genes.sh` and `step_6_filter-te-from-gff.sh`)

- Gene/mRNA summary stats were computed using `agat_sp_statistics.pl` (included in the `step_6_finalize-gff.sh` script).

- Final naming of gene models was perfomed according to maizeGDB guidelines following [MaizeGDB_gff3_format](https://github.com/HuffordLab/MaizeGDB_gff3_format) guide.


