library(purrr)
library(data.table)
library(batchtools)
library(RGCCA)
library(r.jive)
library(tidyverse)
library(dtplyr)

reg =loadRegistry("bootstrap", writeable=FALSE)

ncpus = as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
setDTthreads(ncpus)

bootstrap_results = reduceResultsDataTable(ids=findDone(),
    fun=\(x) x$learner_id)

job_table = getJobTable()[,.(job.id, prob.pars, algorithm, repl)] |> unwrap()

job_table = bootstrap_results %>%
    unwrap() %>%
    left_join(job_table) %>%
    rename(learner_id=result.1) %>%
    tibble()

job_table = job_table %>%
    mutate(
        learner_id=case_when(
            algorithm=="sjive" ~ sub("jive", "sjive", learner_id),
            .default=learner_id
        )
    )

uniq_lrn = job_table %>%
    select(-job.id, -repl) %>%
    #filter(algorithm=="rgcca") %>%
    group_by(learner_id, cancer) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    mutate(uniq_id=1:n())

# sjive has the wrong learner_id, modify it before creating the DAG!!
# (jive_ instead of sjive_)

fwrite(uniq_lrn, "uniq_lrn.csv")
fwrite(job_table, "job_table.csv")

uniq_lrn = fread("uniq_lrn.csv")
job_table = fread("job_table.csv")

resources = list(
    template="./pegasus-simple.tmpl",
    walltime=3600*24*1,
    CpW=1,
    ntasks=73,
    mem="10G",
    logs="logs",
    partition="normal",
    qos="default")

#manually creating pegasus DAG 
dag = vector("character", nrow(uniq_lrn))
for (i in 1:nrow(uniq_lrn)) {
    print(i)
    dag[i] = paste0("TASK ", i, " -c ", resources$CpW,
        " bash -Eeo pipefail -c ",
        "\"Rscript extract_model.R -i ", i, "\"") 
}
fileConn = file(paste0("pegasus_bootstrap.dag"))
writeLines(dag, fileConn)
close(fileConn)

#creating shell script to run pegasus with asked resources
brew_env = new.env()
brew_env$resources = resources
brew_env$job = list(
    registry="bootstrap_results", # only for names
    dag="pegasus_bootstrap.dag"
)
brew::brew(file=resources$template,
    output="pegasus_bootstrap.sh",
    envir=brew_env)
