# Gene molecular evolution analysis

Scripts used below can be found in the `scripts` folder.

The analyses are done for a filtered set orthoFinder-based orthogroups(OGs).
You can find the OG sets in `data` folder

## PHAST-based modeling for molecular evolution associated with perennial-annual transition:
Run script `01A_phyloP_PHASTModeling.sh`
This includes several steps of MSA processing/filtering, neutral evolution model fitting and phyloP calculation.
The output is then summarized and further analyzed in the Rscript `01B_phyloP_resSummary.R`

## Test for non-random occurrence of premature stop codon associated with perennial-annual transition:
The Rscript `02_preMatureStopCodonAnalysis.R` takes in the MSA of each OG and performs Fisher's exact test for the occurrence of premature stop codon occurrence.

## Visualization
Coding for plotting can be found in `scripts`: 
FigA for the phylogenetic tree; 
FigB for the numbers of candidate genes showing the opposite effects.

