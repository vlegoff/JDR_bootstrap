library(biomaRt)
library(fgsea)
library(rbioapi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(tidyverse)
library(dtplyr)
library(msigdbr)

ncpus = as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
setDTthreads(ncpus)

cancers = c("LAML", "ESCA", "PAAD", "SARC", "LIHC", "KIRP", "COAD", "OV", "KIRC",
    "SKCM", "STAD", "BLCA", "UCEC", "LUSC", "LGG", "LUAD", "HNSC", "BRCA")
methods = c("jive_", "jive_clin_fav", "rgcca_sup_", "rgcca_sup_clin_fav",
    "rgcca_unsup_", "rgcca_unsup_clin_fav", "ipca_", "ipca_clin_fav",
    "sjive_", "sjive_clin_fav")

cox_df = list()
jdr_df = list()
for (c in cancers) {
    for (m in methods) {
        message(paste0(c, "\t", m))
        cox_tmp = fread(paste0("results/", c, "_", m, "_cox.csv"))
        cox_df = c(cox_df, list(cox_tmp))
        jdr_tmp = fread(paste0("results/", c, "_", m, ".csv"))
        jdr_tmp$cancer = c
        jdr_df = c(jdr_df, list(jdr_tmp))
    }
}
jdr_df = rbindlist(jdr_df)
cox_df = rbindlist(cox_df)
jdr_df[, p.adj:=p.adjust(pval, method="BH")][, qval:=qvalue::qvalue(pval)$qvalues]
cox_df[, p.adj:=p.adjust(pval, method="BH")]

fwrite(jdr_df, "results/aggregated_jdr_weights.csv")
fwrite(cox_df, "aggregated_cox_weights.csv")


jdr_df = fread("results/aggregated_jdr_weights.csv")
cox_df = fread("results/aggregated_cox_weights.csv")

jdr_for_barplot = jdr_df[,
    .(n=sum(p.adj<=0.05, na.rm=TRUE)), by=.(cancer, learner_id, comp)
]
fwrite(jdr_for_barplot, "jdr_for_barplot.csv")
cox_for_barplot = cox_df[,
    .(n=sum(p.adj<=0.05, na.rm=TRUE)), by=.(cancer, learner_id, comp)
]
fwrite(cox_for_barplot, "cox_for_barplot.csv")

jdr_sig = jdr_df[p.adj<=0.05, .N, by=.(learner_id, cancer, block, comp)]
sig_comps = cox_df %>%
    subset(select=c(learner_id, cancer, comp, p.adj)) %>%
    filter(p.adj<=0.05) %>%
    left_join(jdr_sig %>%
        mutate(comp=sub("joint", "common", comp))
    ) %>%
    filter(!is.na(N)) %>%
    subset(select=c(cancer, learner_id, comp, p.adj, N)) %>%
    mutate(cancer=forcats::fct_relevel(cancer, cancers)) %>%
    arrange(cancer, learner_id, comp)
sig_comps_plus = sig_comps %>%
    left_join(select(cox_df, cancer, learner_id, comp, weight_mean, weight_sd, n)) %>%
    rename(
        n_genes_sig=N,
        cox_mean=weight_mean,
        cox_sd=weight_sd,
        n_bootstrap_sample=n
    ) %>%
    left_join(
        enrich_df %>%
            filter(gene_set %in% c("KEGG", "TAM")) %>%
            mutate(comp=case_when(
                grepl("common", comp) ~ sub("(rna|mirna|mutation|cnv)\\.", "", comp),
                .default=comp
            )) %>%
            mutate(padj=p.adjust(pval, method="BH")) %>%
            group_by(cancer, learner_id, comp) %>%
            summarize(n_pathway_signif=sum(padj<=0.05), n_pathway_total=n()) %>%
            ungroup() %>%
            mutate(prop_pathway_signif=n_pathway_signif/n_pathway_total)
    ) %>%
    select(cancer, learner_id, comp, p.adj, cox_mean, cox_sd, n_genes_sig,
        n_pathway_signif, prop_pathway_signif, n_pathway_total, n_bootstrap_sample)
fwrite(sig_comps_plus, "cox_comps_sig.csv")

gene_map = fread("../wahdidou/vlegoff_gene_id_map.csv")
dup_ens = gene_map |>
    _[["ensembl_gene_id"]] |>
    table() |>
    data.frame() |>
    subset(Freq>1) |>
    _[["Var1"]] |>
    as.character()
dup_sym = gene_map |>
    _[["orig_symbol"]] |>
    table() |>
    data.frame() |>
    subset(Freq>1) |>
    _[["Var1"]] |>
    as.character()
gene_map_uniq = gene_map |>
    subset(!(ensembl_gene_id %in% dup_ens) & !(orig_symbol %in% dup_sym)) |>
    data.table()
mirna_map = fread("../wahdidou/vlegoff_mirna_gene_map.csv")
uplift_map = fread("../id_conversion/uplift_ensg_map.csv")

#gencode = "/env/cnrgh/proj/math_stats/scratch/gloaguen/These_VLG/enrichment/gencode.v22.annotation.gtf"
url_gencode22 = "ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz"
file_gencode22 = "gencode.v22.annotation.gtf.gz"
download.file(url_gencode22, destfile=file_gencode22)

gtf            = rtracklayer::import(gzfile(file_gencode))
gtf_df         = as.data.frame(gtf)
gtf_df$gene_id = sapply(strsplit(x = gtf_df$gene_id, split = "\\."), function(x) x[[1]])

jdr_df %>% 
  filter(cancer == "LGG") %>%
  filter(learner_id %in% c("jive_", "jive_clin_fav", "rgcca_sup_", "rgcca_unsup_")) %>%
  filter(grepl("_cnv$", var)) %>%
  mutate(var=sub("_cnv$", "", var)) %>%
  left_join(gtf_df %>%
    filter(type=="gene") %>%
    select(gene_id, seqnames, start, end) %>%
    rename(var=gene_id, chr=seqnames)
  ) %>%
  fwrite("LGG_cnv.csv")

jdr_df %>% 
  filter(cancer == "LGG") %>%
  filter(learner_id %in% c("jive_", "jive_clin_fav", "rgcca_sup_", "rgcca_unsup_")) %>%
  filter(grepl("_rna$", var)) %>%
  mutate(var=sub("_rna$", "", var)) %>%
  left_join(gtf_df %>%
    filter(type=="gene") %>%
    select(gene_id, seqnames, start, end) %>%
    rename(var=gene_id, chr=seqnames)
  ) %>%
  fwrite("LGG_rna.csv")

jdr_df %>% 
  filter(cancer == "LGG") %>%
  filter(grepl("_mutation$", var)) %>%
  mutate(var=sub("_mutation$", "", var)) %>%
  left_join(gene_map_uniq %>%
    select(orig_symbol, ensembl_gene_id) %>%
    rename(var=orig_symbol)) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  rename(symbol=var) %>%
  rename(var=ensembl_gene_id) %>%
  left_join(gtf_df %>%
    filter(type=="gene") %>%
    select(gene_id, seqnames, start, end) %>%
    rename(var=gene_id, chr=seqnames)
  ) %>%
  fwrite("LGG_mutation.csv")

genesets_H = msigdbr(species="Homo sapiens", collection="H")
fwrite(genesets_H, "genesets_H.csv")
genesets = genesets_H[, c("gs_name", "ensembl_gene")]
tam_functions_for_profiler = clusterProfiler::read.gmt("TAM_Function_precursor.gmt")

#sig_to_enrich = sig_comps[!grepl("mirna", comp)][!grepl("common", comp)]
sig_to_enrich = sig_comps[!grepl("mirna", comp)][grepl("rna", comp)]
ekegg_df = list()
for (i in 1:nrow(sig_to_enrich)) {
    message(paste0(i, "/", nrow(sig_to_enrich), "\t", sig_to_enrich[i,]$comp))
    cox_weight = cox_df %>%
        filter(cancer==sig_to_enrich[i,]$cancer,
               learner_id==sig_to_enrich[i,]$learner_id,
               comp==sig_to_enrich[i,]$comp) %>%
        pull(weight_mean)
    comp_for_gsea = jdr_df[
        cancer==sig_to_enrich[i,]$cancer & 
            learner_id==sig_to_enrich[i,]$learner_id & 
            comp==sig_to_enrich[i,]$comp,
        #.(var=sub("_(rna|cnv|mutation)$", "", var), stat=-1*sign(cox_weight)*sign(weight_mean)*log(pval))
        .(var=sub("_(rna|cnv|mutation)$", "", var), stat=bts_ratio)
    ][
        order(stat, decreasing=TRUE)
    ]
    if (grepl("^mutation", sig_to_enrich[i,]$comp)) {
        comp_for_gsea = gene_map_uniq[,
            .(var=orig_symbol, ensembl_gene_id)
        ][
            comp_for_gsea, on=.NATURAL
        ][,
            .(var=ensembl_gene_id, stat)
        ][
            !is.na(var)
        ]
    }
    comp_for_gsea = comp_for_gsea[
        order(stat, decreasing=TRUE), .(var, stat)
    ]
    vec_for_gsea = comp_for_gsea$stat
    names(vec_for_gsea) = comp_for_gsea$var
    ekegg = GSEA(
        geneList=vec_for_gsea,
        TERM2GENE=genesets,
        pvalueCutoff=1,
        eps=0
    )
    ekegg_df[[i]] = as.data.frame(ekegg)
    if(nrow(ekegg_df[[i]])>0) {
        ekegg_df[[i]]$cancer = sig_to_enrich[i,]$cancer
        ekegg_df[[i]]$learner_id = sig_to_enrich[i,]$learner_id
        ekegg_df[[i]]$comp = sig_to_enrich[i,]$comp
    }
}

gseaplot(ekegg, geneSetID="HALLMARK_INTERFERON_ALPHA_RESPONSE")

# common components
sig_to_enrich = sig_comps[grepl("common", comp)]
ecommon_df = list()
for (i in 1:nrow(sig_to_enrich)) {
    message(paste0(i, "/", nrow(sig_to_enrich), "\t", sig_to_enrich[i,]$comp))
    #cox_weight = cox_df %>%
    #    filter(cancer==sig_to_enrich[i,]$cancer,
    #           learner_id==sig_to_enrich[i,]$learner_id,
    #           comp==sig_to_enrich[i,]$comp) %>%
    #    pull(weight_mean)
    comp_to_select = sub("common", "joint", sig_to_enrich[i,]$comp)
    full_comp = jdr_df[
        cancer==sig_to_enrich[i,]$cancer & 
            learner_id==sig_to_enrich[i,]$learner_id & 
            comp==comp_to_select
    ]
    ecommon_tmp = list()
    #for (b in c("rna", "mutation", "cnv")) {
    for (b in c("rna")) {
        message(b)
        comp_for_gsea = full_comp[
            grepl(paste0("_", b, "$"), var),
            #.(var=sub("_(rna|cnv|mutation)$", "", var), stat=-1*sign(cox_weight)*sign(weight_mean)*log(pval))
            .(var=sub("_(rna|cnv|mutation)$", "", var), stat=bts_ratio)
        ][
            order(stat, decreasing=TRUE)
        ]
        if (b=="mutation") {
            comp_for_gsea = gene_map_uniq[,
                .(var=orig_symbol, ensembl_gene_id)
            ][
                comp_for_gsea, on=.NATURAL
            ][,
                .(var=ensembl_gene_id, stat)
            ][
                !is.na(var)
            ]
        }
        comp_for_gsea = comp_for_gsea[
            order(stat, decreasing=TRUE), .(var, stat)
        ]
        vec_for_gsea = comp_for_gsea$stat
        names(vec_for_gsea) = comp_for_gsea$var
        # KEGG
        ecommon = GSEA(
            geneList=vec_for_gsea,
            TERM2GENE=genesets,
            pvalueCutoff=1,
            eps=0
        )
        ecommon_tmp[[b]] = as.data.frame(ecommon)
        ecommon_tmp[[b]]$gene_set = "KEGG"
        if(nrow(ecommon_tmp[[b]])>0) {
            ecommon_tmp[[b]]$cancer = sig_to_enrich[i,]$cancer
            ecommon_tmp[[b]]$learner_id = sig_to_enrich[i,]$learner_id
            ecommon_tmp[[b]]$comp = paste0(b, ".", sig_to_enrich[i,]$comp)
        }
    }
    message("mirna")
    comp_for_gsea = full_comp[
        grepl("_mirna$", var),
        #.(var=sub("_mirna$", "", var), stat=-1*sign(cox_weight)*sign(weight_mean)*log(pval))
        .(var=sub("_mirna$", "", var), stat=bts_ratio)
    ][
        order(stat, decreasing=TRUE)
    ][,
        .(var=gsub("\\.", "-", var), stat)  
    ]
    vec_for_gsea = comp_for_gsea$stat
    names(vec_for_gsea) = comp_for_gsea$var
    emirna = GSEA(
        geneList=vec_for_gsea,
        pvalueCutoff=1,
        TERM2GENE=tam_functions_for_profiler,
        minGSSize=10,
        eps=0
    )
    ecommon_tmp[["mirna"]] = as.data.frame(emirna)
    ecommon_tmp[["mirna"]]$gene_set = "TAM"
    if(nrow(ecommon_tmp[["mirna"]])>0) {
        ecommon_tmp[["mirna"]]$cancer = sig_to_enrich[i,]$cancer
        ecommon_tmp[["mirna"]]$learner_id = sig_to_enrich[i,]$learner_id
        ecommon_tmp[["mirna"]]$comp = paste0("mirna.", sig_to_enrich[i,]$comp)
    }
    ecommon_df[[i]] = rbindlist(ecommon_tmp)
}

sig_mirna = sig_comps[grepl("mirna", comp)]
emirna_df = list()
for (i in 1:nrow(sig_mirna)) {
    message(paste0(i, "/", nrow(sig_mirna), "\t", sig_mirna[i,]$comp))
    cox_weight = cox_df %>%
        filter(cancer==sig_to_enrich[i,]$cancer,
               learner_id==sig_to_enrich[i,]$learner_id,
               comp==sig_to_enrich[i,]$comp) %>%
        pull(weight_mean)
    comp_for_gsea = jdr_df[
        cancer==sig_mirna[i,]$cancer & 
            learner_id==sig_mirna[i,]$learner_id & 
            comp==sig_mirna[i,]$comp,
        #.(var=sub("_mirna$", "", var), stat=-1*sign(cox_weight)*sign(weight_mean)*log(pval))
        .(var=sub("_mirna$", "", var), stat=bts_ratio)
    ][
        order(stat, decreasing=TRUE)
    ][,
        .(var=gsub("\\.", "-", var), stat)  
    ]
    vec_for_gsea = comp_for_gsea$stat
    names(vec_for_gsea) = comp_for_gsea$var
    emirna = GSEA(
        geneList=vec_for_gsea,
        pvalueCutoff=1,
        TERM2GENE=tam_functions_for_profiler,
        minGSSize=10,
        eps=0
    )
    emirna_df[[i]] = as.data.frame(emirna)
    emirna_df[[i]]$gene_set = "TAM"
    if(nrow(emirna_df[[i]])>0) {
        emirna_df[[i]]$cancer = sig_mirna[i,]$cancer
        emirna_df[[i]]$learner_id = sig_mirna[i,]$learner_id
        emirna_df[[i]]$comp = sig_mirna[i,]$comp
    }
}

ecommon_df = rbindlist(ecommon_df)
emirna_df = rbindlist(emirna_df)
ekegg_df = rbindlist(ekegg_df)
enrich_df = rbind(mutate(ekegg_df, gene_set="H"), emirna_df, ecommon_df) %>%
    select(cancer, learner_id, comp, Description, enrichmentScore,
           NES, setSize, pvalue, core_enrichment, gene_set) %>%
    mutate(leading_edge=core_enrichment) %>%
    mutate(core_enrichment=sapply(strsplit(core_enrichment, "/"), length)) %>%
    mutate(gene_set=ifelse(gene_set=="KEGG", "H", gene_set)) %>%
    rename(
        pathway=Description,
        ES=enrichmentScore,
        size=setSize,
        pval=pvalue
    )
fwrite(enrich_df, "enrichment_H.csv")
