#!/bin/bash
# run mikado pick
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus=${cpus:-36}
targets=$(find $(pwd) -name "uniprot-sprot_viridiplantae.fasta")
blastxml=$(find $(pwd) -name "mikado.blast.xml")
orfs=$(find $(pwd) -name "mikado_prepared.fasta.transdecoder.bed")
# run configure
mikado serialise --start-method spawn --procs ${cpus} --blast_targets ${targets} --json-conf config.toml --xml ${blastxml} --orfs ${orfs}
