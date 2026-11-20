library(purrr)
library(data.table)
library(batchtools)
library(RGCCA)
library(r.jive)
library(tidyverse)
library(dtplyr)

ncpus = as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
setDTthreads(ncpus)

cancers = c(
    "LAML",
    "ESCA",
    "PAAD",
    "SARC",
    "LIHC",
    "COAD",
    "KIRP",
    "OV",
    "KIRC",
    "SKCM",
    "STAD",
    "BLCA",
    "UCEC",
    "LGG",
    "HNSC",
    "LUSC",
    "LUAD",
    "BRCA"
)
cancers = sort(cancers)
#dat = map(paste0("../benchmark/data/", cancers, ".csv"), data.table::fread)
dat = paste0("data/herrmann2021/", cancers, ".csv")
names(dat) = cancers

# get weights from previous registry
reg = loadRegistry("bootstrap_orig_models", writeable=FALSE)

weights = reduceResultsDataTable() %>%
    left_join(unwrap(getJobTable()[, .(job.id, problem, algorithm, algo.pars, prob.pars)])) %>%
    rename(weights_orig=result)

rgcca_params = weights %>%
    filter(algorithm=="rgcca") %>%
    select(-j, -a, -eta, -job.id, -problem) %>%
    group_by(cancer, clin, sup) %>%
    group_by(clin, sup) %>%
    arrange(cancer, .by_group=TRUE) %>%
    ungroup()
jive_params = weights %>%
    filter(algorithm=="jive") %>%
    select(-n, -tau, -sup, -eta, -job.id, -problem) %>%
    group_by(cancer, clin) %>%
    group_by(clin) %>%
    arrange(cancer, .by_group=TRUE) %>%
    ungroup()
ipca_params = weights %>%
    filter(algorithm=="ipca") %>%
    select(-j, -n, -tau, -sup, -eta, -job.id, -problem) %>%
    group_by(cancer, clin) %>%
    group_by(clin) %>%
    arrange(cancer, .by_group=TRUE) %>%
    ungroup()
sjive_params = weights %>%
    filter(algorithm=="sjive") %>%
    select(-n, -tau, -sup, -job.id, -problem) %>%
    group_by(cancer, clin) %>%
    group_by(clin) %>%
    arrange(cancer, .by_group=TRUE) %>%
    ungroup()
ades = list(rgcca=rgcca_params %>% select(-algorithm, -cancer),
    jive=jive_params %>% select(-algorithm, -cancer),
    sjive=sjive_params %>% select(-algorithm, -cancer),
    ipca=ipca_params %>% select(-algorithm, -cancer)
)

#reg = makeExperimentRegistry(file.dir="bootstrap", seed=1)
reg =loadRegistry("bootstrap", writeable=TRUE)

clearRegistry()
multiomic_bootstrap = function(data, job, cancer) {
    data = data.table::fread(paste0("data/herrmann2021/", cancer, ".csv"))
    bts_smp = c(
        sample(which(data$status==0), replace=TRUE),
        sample(which(data$status==1), replace=TRUE)
    )
    data = data[bts_smp,]
    time = data$time
    status = data$status
    ysurv = survival::Surv(time=time, event=status)
    blocks = c("mutation", "rna", "mirna", "cnv", "clinical")
    x_list = lapply(blocks,
        \(bl) data[,.SD, .SDcols=grepl(paste0("_", bl), colnames(data))])
    names(x_list) = blocks
    x_list = lapply(x_list, as.matrix)
    clin_data = x_list[["clinical"]]
    x_list = x_list[c("mutation", "rna", "mirna", "cnv")]
    return(list(x=x_list, clin=clin_data, y=ysurv))
}
addProblem(name="bootstrap", data=NULL, fun=multiomic_bootstrap, seed=42)

source("jdr_functions.R")
addAlgorithm(name="rgcca", fun=rgcca_wrapper)
addAlgorithm(name="jive", fun=jive_wrapper)
addAlgorithm(name="sjive", fun=sjive_wrapper)
addAlgorithm(name="ipca", fun=ipca_wrapper)

pdes = list(bootstrap=data.frame(cancer=cancers))
addExperiments(prob.designs=pdes, algo.designs=ades, repls=99, combine="bind")
