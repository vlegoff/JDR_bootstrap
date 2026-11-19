library(purrr)
library(data.table)
library(batchtools)
library(RGCCA)
library(r.jive)
library(tidyverse)
library(dtplyr)

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
dat = map(paste0("../benchmark/data/", cancers, ".csv"), data.table::fread)
names(dat) = cancers

#reg = makeExperimentRegistry(file.dir="bootstrap_orig_models", seed=1)
reg =loadRegistry("bootstrap_orig_models", writeable=TRUE)

#clearRegistry()
multiomic_no_bootstrap = function(data, job, cancer) {
    data = data.table::fread(paste0("../benchmark/data/", cancer, ".csv"))
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
addProblem(name="no_bootstrap", data=NULL, fun=multiomic_no_bootstrap, seed=42)

source("jdr_functions.R")
addAlgorithm(name="rgcca", fun=rgcca_wrapper)
addAlgorithm(name="jive", fun=jive_wrapper)
addAlgorithm(name="ipca", fun=ipca_wrapper)
addAlgorithm(name="sjive", fun=sjive_wrapper)

rgcca_params = fread("rgcca_params.csv") %>%
    filter(task_id %in% cancers) %>%
    mutate(
        clin=ifelse(clin=="with", TRUE, FALSE),
        sup=ifelse(sup=="sup", TRUE, FALSE)
    ) %>%
    group_by(clin, sup) %>%
    arrange(task_id, .by_group=TRUE) %>%
    ungroup()
jive_params = fread("jive_params.csv") %>%
    filter(task_id %in% cancers) %>%
    mutate(
        clin=ifelse(clin=="with", TRUE, FALSE)
    ) %>%
    group_by(clin) %>%
    arrange(task_id, .by_group=TRUE) %>%
    ungroup()
sjive_params = fread("sjive_params.csv") %>%
    filter(task_id %in% cancers) %>%
    mutate(
        clin=ifelse(clin=="with", TRUE, FALSE)
    ) %>%
    group_by(clin) %>%
    arrange(task_id, .by_group=TRUE) %>%
    ungroup()
ipca_params = fread("ipca_params.csv") %>%
    filter(task_id %in% cancers) %>%
    mutate(
        clin=ifelse(clin=="with", TRUE, FALSE)
    ) %>%
    group_by(clin) %>%
    arrange(task_id, .by_group=TRUE) %>%
    ungroup()

ades = list(rgcca=select(rgcca_params, -task_id),
            jive=select(jive_params, -task_id),
            sjive=select(sjive_params, -task_id),
            ipca=select(ipca_params, -task_id))
pdes = list(no_bootstrap=data.frame(cancer=cancers))
addExperiments(prob.designs=pdes, algo.designs=ades, repls=1, combine="bind")
