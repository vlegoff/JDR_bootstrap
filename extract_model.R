library(purrr)
library(data.table)
library(batchtools)
library(RGCCA)
library(r.jive)
library(tidyverse)
library(dtplyr)
library(argparse)

parser = ArgumentParser()
parser$add_argument("-i", type="integer", default=1)
args = parser$parse_args()
i = args$i

reg =loadRegistry("bootstrap_rgcca_only", writeable=FALSE)
uniq_lrn = fread("uniq_lrn.csv")
job_table = fread("job_table.csv")

cancer = uniq_lrn[i,]$cancer
learner_id = uniq_lrn[i,]$learner_id
ids = job_table %>%
    filter(
        learner_id==uniq_lrn[i,]$learner_id,
        cancer==uniq_lrn[i,]$cancer
    ) %>%
    pull(job.id)
ids = intersect(ids, findDone(ids)$job.id)
tmp_results = list()
tmp_results_cox = list()
for (j in seq_along(ids)) {
    message(paste0("model: ", i, "/", nrow(uniq_lrn), "\tbootstrap: ", j, "/", length(ids)))
    id = ids[[j]]
    weights = reduceResultsList(id)[[1]]
    cox_dt = data.table(weights$cox)
    colnames(cox_dt) = c("weight")
    cox_dt$type = "cox"
    tmp = strsplit(names(weights$cox), "\\.")
    cox_dt$comp = names(weights$cox)
    cox_dt$cancer = cancer
    cox_dt$learner_id = learner_id
    jdr_dt = list()
    if (weights$model=="rgcca" | weights$model=="ipca") {
        for (b in c("cnv", "rna", "mirna", "mutation")) {
            colnames(weights$jdr[[b]]) = paste0(b, ".", 1:ncol(weights$jdr[[b]]))
            jdr_dt[[b]] = weights$jdr[[b]] %>%
                data.table() %>%
                mutate(var=rownames(weights$jdr[[b]])) %>%
                pivot_longer(!var, names_to="comp", values_to="weight") %>%
                data.table()
            jdr_dt[[b]]$block = b
            jdr_dt[[b]]$learner_id = learner_id
            jdr_dt[[b]]$cancer = learner_id
        }
    } else if (weights$model=="jive" | weights$model=="sjive") { # check things
        if (!is.null(weights$jdr$joint)) {
            b = "joint"
            colnames(weights$jdr$joint) = paste0(b, ".", 1:ncol(weights$jdr$joint))
            jdr_dt[[b]] = weights$jdr$joint %>%
                data.table() %>%
                mutate(var=rownames(weights$jdr[[b]])) %>%
                pivot_longer(!var, names_to="comp", values_to="weight") %>%
                data.table()
            jdr_dt[[b]]$block = "joint"
            jdr_dt[[b]]$learner_id = learner_id
            jdr_dt[[b]]$cancer = learner_id
        }
        if (!is.null(weights$jdr$indiv)) {
            for (b in c("cnv", "rna", "mirna", "mutation")) {
                colnames(weights$jdr$indiv[[b]]) = paste0(b, ".", 1:ncol(weights$jdr$indiv[[b]]))
                jdr_dt[[b]] = weights$jdr$indiv[[b]] %>%
                    data.table() %>%
                    mutate(var=rownames(weights$jdr$indiv[[b]])) %>%
                    pivot_longer(!var, names_to="comp", values_to="weight") %>%
                    data.table()
                jdr_dt[[b]]$block = b
                jdr_dt[[b]]$learner_id = learner_id
                jdr_dt[[b]]$cancer = learner_id
            }
        }
    }
    tmp_results = c(tmp_results, list(rbindlist(jdr_dt)))
    tmp_results_cox = c(tmp_results_cox, list(cox_dt))
}
tmp_results = rbindlist(tmp_results)
tmp_results_cox = rbindlist(tmp_results_cox)
tmp_results %>%
    group_by(learner_id, cancer, block, comp, var) %>%
    summarize(
        weight_sd=sd(c(weight, rep(0, 99-n()))),
        weight_mean=mean(c(weight, rep(0, 99-n()))),
        n=n()
    ) %>%
    ungroup() %>%
    mutate(bts_ratio = weight_mean/weight_sd) %>%
    mutate(pval = 2*pnorm(abs(bts_ratio), lower.tail=FALSE)) %>%
    write_csv(file=paste0("results/", cancer, "_", learner_id, ".csv"))
tmp_results_cox %>%
    group_by(learner_id, cancer, comp) %>%
    summarize(weight_sd=sd(weight), weight_mean=mean(weight), n=n()) %>%
    ungroup() %>%
    mutate(bts_ratio = weight_mean/weight_sd) %>%
    mutate(pval = 2*pnorm(abs(bts_ratio), lower.tail=FALSE)) %>%
    write_csv(file=paste0("results/", cancer, "_", learner_id, "_cox.csv"))
