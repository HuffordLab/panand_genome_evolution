#!/bin/bash
cpus=$SLURM_JOB_CPUS_PER_NODE
# run prepare
mikado prepare --procs ${cpus} --json-conf config.toml

