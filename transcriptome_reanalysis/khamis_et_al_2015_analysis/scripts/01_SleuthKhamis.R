library(sleuth)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lemon)
library(cowplot)
library(xlsx)
base_dir <- "../khamis_et_al_2015_analysis/kallisto"
sample_id <- dir (file.path(base_dir))
metadata <- data.frame(sample = sample_id,condition = factor(c(rep("Forrager",8),rep("Nurse",8))))
paths <- file.path(base_dir, sample_id)
paths
names(paths) <- sample_id
paths
metadata$path <- paths
metadata
metadata$condition <- relevel(metadata$condition, ref = "Nurse")
design <- ~ condition
so <- sleuth_prep(metadata, full_model = design, num_cores = 4L, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, transform_fun_counts = function (x) log2(x+0.5))
so <- sleuth_fit(so)
models(so)
DE_Forrager <- sleuth_wt(so, which_beta = 'conditionForrager')
sleuth_results_Forrager <- sleuth_results(DE_Forrager,test = 'conditionForrager', show_all = TRUE)
write.xlsx(sleuth_results_Forrager, file='../transcriptome_reanalysis/khamis_et_al_2015_analysis/results_sleuth/sleuth_results_Khamis.xlsx', sheetName = "DE_Forrager")

