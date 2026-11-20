library(data.table)
library(tidyverse)
library(dtplyr)
library(purrr)

data_dir = "../benchmark/results/"
results_dir = "../benchmark/results/"
registries = c(
    "RGCCA",
    "RGCCA_woClin",
    "JIVE",
    "JIVE_woClin",
    "SJIVE",
    "SJIVE_woClin",
    "PCA",
    "PCA_woClin"
)
cv_files = paste0(data_dir, "/cv_logs/", registries, "/raw_cv_results.csv")

cv_df_raw = rbindlist(map(cv_files, fread), fill=TRUE)
cv_df = fread(paste0(results_dir, "/results_cross_validated.csv"))

cv_df = cv_df %>%
    mutate(clin=ifelse(grepl("clin_fav", learner_id), "with", "without")) %>%
    mutate(learner=sub("clin_fav_", "", learner_id))
cv_df_raw = cv_df_raw %>%
    mutate(clin=ifelse(grepl("clin_fav", learner_id), "with", "without")) %>%
    mutate(learner=sub("clin_fav_", "", learner_id))

cv_jive = cv_df %>%
    filter(grepl("jive", learner_id)) %>%
    filter(!grepl("sjive", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("jive", "_j", j=".*", "_a", a=".*")) %>%
    group_by(clin, task_id, j, a, learner_id) %>%
    summarize(count=n()) %>%
    ungroup() %>%
    slice_max(count, by=c(clin, task_id)) %>%
    left_join(cv_df_raw %>% select(learner_id, task_id, lambda, iteration)) %>%
    group_by(clin, task_id, j, a) %>%
    summarize(lambda=median(lambda))

# separating ties from the rest
cv_jive_ties = cv_jive %>%
    group_by(clin, task_id) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    filter(n>1) %>%
    left_join(cv_jive)
cv_jive_no_ties = cv_jive %>%
    group_by(clin, task_id) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    filter(n==1) %>%
    left_join(cv_jive)

cv_jive = cv_df_raw %>%
    filter(grepl("jive", learner_id)) %>%
    filter(!grepl("sjive", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("jive", "_j", j=".*", "_a", a=".*")) %>%
    group_by(clin, task_id, j, a, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    select(-learner_id, -lambda) %>%
    inner_join(cv_jive_ties) %>%
    slice_min(cv_metric, by=c(clin, task_id)) %>%
    select(-cv_metric) %>%
    rbind(cv_jive_no_ties) %>%
    arrange(clin, task_id) %>%
    select(-n)


cv_jive2 = cv_df_raw %>%
    filter(grepl("jive", learner_id)) %>%
    filter(!grepl("sjive", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("jive", "_j", j=".*", "_a", a=".*")) %>%
    group_by(clin, task_id, j, a, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    slice_min(cv_metric, by=c(clin, task_id)) %>%
    select(-learner_id, -cv_metric)

cv_df_raw %>% filter(learner_id=="jive_j5_a0")

cv_sjive = cv_df %>%
    filter(grepl("sjive", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("sjive", "_j", j=".*", "_a", a=".*", "_eta", eta=".*")) %>%
    group_by(clin, task_id, j, a, eta, learner_id) %>%
    summarize(count=n()) %>%
    ungroup() %>%
    slice_max(count, by=c(clin, task_id)) %>%
    left_join(cv_df_raw %>% select(learner_id, task_id, lambda, iteration)) %>%
    group_by(clin, task_id, j, a, eta) %>%
    summarize(lambda=median(lambda))

# separating ties from the rest
cv_sjive_ties = cv_sjive %>%
    group_by(clin, task_id) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    filter(n>1) %>%
    left_join(cv_sjive)
cv_sjive_no_ties = cv_sjive %>%
    group_by(clin, task_id) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    filter(n==1) %>%
    left_join(cv_sjive)

cv_sjive = cv_df_raw %>%
    filter(grepl("sjive", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("sjive", "_j", j=".*", "_a", a=".*", "_eta", eta=".*")) %>%
    group_by(clin, task_id, j, a, eta, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    select(-learner_id, -lambda) %>%
    inner_join(cv_sjive_ties) %>%
    slice_min(cv_metric, by=c(clin, task_id)) %>%
    select(-cv_metric) %>%
    rbind(cv_sjive_no_ties) %>%
    arrange(clin, task_id) %>%
    select(-n)

cv_rgcca = cv_df %>%
    filter(grepl("rgcca", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("rgcca_", sup=".*", "_n", n=".*",
                   "_tau", tau=".*")) %>%
    group_by(clin, task_id, sup, n, tau, learner_id) %>%
    summarize(count=n()) %>%
    ungroup() %>%
    slice_max(count, by=c(clin, task_id, sup)) %>%
    left_join(cv_df_raw %>% select(learner_id, task_id, lambda, iteration)) %>%
    group_by(clin, task_id, sup, n, tau) %>%
    summarize(lambda=median(lambda))
cv_rgcca_ties = cv_rgcca %>%
    group_by(clin, sup, task_id) %>%
    summarize(ties=n()) %>%
    ungroup() %>%
    filter(ties>1) %>%
    left_join(cv_rgcca)
cv_rgcca_no_ties = cv_rgcca %>%
    group_by(clin, sup, task_id) %>%
    summarize(ties=n()) %>%
    ungroup() %>%
    filter(ties==1) %>%
    left_join(cv_rgcca)
cv_rgcca = cv_df_raw %>%
    filter(grepl("rgcca", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("rgcca_", sup=".*", "_n", n=".*",
                   "_tau", tau=".*")) %>%
    group_by(clin, task_id, sup, n, tau, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    select(-learner_id, -lambda) %>%
    inner_join(cv_rgcca_ties) %>%
    slice_min(cv_metric, by=c(clin, task_id)) %>%
    select(-cv_metric) %>%
    rbind(cv_rgcca_no_ties) %>%
    arrange(clin, task_id, sup) %>%
    select(-ties)

cv_rgcca2 = cv_df_raw %>%
    filter(grepl("rgcca", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("rgcca_", sup=".*", "_n", n=".*",
                   "_tau", tau=".*")) %>%
    group_by(clin, task_id, sup, n, tau, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    slice_min(cv_metric, by=c(clin, task_id, sup)) %>%
    select(-learner_id, -cv_metric)

cv_ipca = cv_df %>%
    filter(grepl("ipca", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("ipca_a", a=".*")) %>%
    group_by(clin, task_id, a, learner_id) %>%
    summarize(count=n()) %>%
    ungroup() %>%
    slice_max(count, by=c(clin, task_id)) %>%
    left_join(cv_df_raw %>% select(learner_id, task_id, lambda, iteration)) %>%
    group_by(clin, task_id, a) %>%
    summarize(lambda=median(lambda))

cv_ipca_ties = cv_jive %>%
    group_by(clin, task_id) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    filter(n>1) %>%
    left_join(cv_ipca)
cv_ipca_no_ties = cv_jive %>%
    group_by(clin, task_id) %>%
    summarize(n=n()) %>%
    ungroup() %>%
    filter(n==1) %>%
    left_join(cv_ipca)

cv_ipca = cv_df_raw %>%
    filter(grepl("ipca", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("ipca_a", a=".*")) %>%
    group_by(clin, task_id, a, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    select(-learner_id, -lambda) %>%
    inner_join(cv_ipca_ties) %>%
    slice_min(cv_metric, by=c(clin, task_id)) %>%
    select(-cv_metric) %>%
    rbind(cv_ipca_no_ties) %>%
    arrange(clin, task_id) %>%
    select(-n)

cv_ipca2 = cv_df_raw %>%
    filter(grepl("ipca", learner_id)) %>%
    separate_wider_regex(learner,
        patterns=c("ipca_a", a=".*")) %>%
    group_by(clin, task_id, a, learner_id) %>%
    summarize(across(c(cv_metric, lambda), \(x) median(x, na.rm=TRUE))) %>%
    ungroup() %>%
    slice_min(cv_metric, by=c(clin, task_id)) %>%
    select(-learner_id, -cv_metric)

fwrite(cv_jive, "data/jive_params.csv")
fwrite(cv_sjive, "data/sjive_params.csv")
fwrite(cv_rgcca, "data/rgcca_params.csv")
fwrite(cv_ipca, "data/ipca_params.csv")
