#!/bin/bash
#MSUB -q normal
#MSUB -n 1
#MSUB -c 8
#MSUB -E --mem-per-cpu=10G
#MSUB -T 7200
#MSUB -o logs/job.%I.out
#MSUB -e logs/job.%I.err

ENV=../benchmark/env/inti_el9
module load micromamba
micromamba activate -p $ENV

Rscript run_full_models.R
