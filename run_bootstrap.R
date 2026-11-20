library(batchtools)
library(benchbatch)
library(data.table)
ncpus = as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
setDTthreads(ncpus)
reg =loadRegistry("bootstrap", writeable=TRUE)
resources = list(
    walltime=60*60*24,
    ntasks=300,
    mem="10G",
    log="logs",
    partition="normal",
    template="../benchmark/hpc/pegasus-simple.tmpl")
benchbatch::run(reg, cluster="pegasus", resources=resources, resubmit_not_finished=T,
    resubmit_errors=F)
