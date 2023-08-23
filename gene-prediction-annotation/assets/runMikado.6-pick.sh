#!/bin/bash
cpus=$SLURM_JOB_CPUS_PER_NODE
cpus=${cpus:-36}
# run configure
mikado pick --start-method spawn --procs ${cpus} --json-conf config.toml --subloci-out mikado_subloci.gff3


