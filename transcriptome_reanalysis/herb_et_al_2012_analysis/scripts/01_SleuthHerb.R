library(sleuth)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lemon)
library(cowplot)
library(xlsx)

base_dir <- "../herb_et_al_2012_analysis/kallisto/"
sample_id <- dir (file.path(base_dir))
metadata <- data.frame(sample = sample_id,condition = factor(c(rep("RN",6),rep("FO",6))))
paths <- file.path(base_dir, sample_id)
paths
names(paths) <- sample_id
paths
metadata$path <- paths
metadata
metadata$condition <- relevel(metadata$condition, ref = "FO")
design <- ~ condition
so <- sleuth_prep(metadata, full_model = design, num_cores = 4L, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, transform_fun_counts = function (x) log2(x+0.5))
so <- sleuth_fit(so)
models(so)
DE_RN <- sleuth_wt(so, which_beta = 'conditionRN')
sleuth_results_RN <- sleuth_results(DE_RN,test = 'conditionRN', show_all = TRUE)
write.xlsx(sleuth_results_RN, file='../transcriptome_reanalysis/herb_et_al_2012_analysis/results_sleuth/sleuth_results_Herb.xlsx', sheetName = "DE_RN")

