library(sleuth)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(lemon)
library(cowplot)
library(xlsx)
base_dir <- "../warner_et_al_2019_analysis/kallisto"
sample_id <- dir (file.path(base_dir))

####L4 larvae

metadata <- data.frame(sample = sample_id,condition = factor(c(rep("LQ4",3),rep("LW4",3))))
paths <- file.path(base_dir, sample_id)
paths
names(paths) <- sample_id
paths
metadata$path <- paths
metadata
metadata$condition <- relevel(metadata$condition, ref = "LW4")
design <- ~ condition
so <- sleuth_prep(metadata, full_model = design, num_cores = 4L, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, transform_fun_counts = function (x) log2(x+0.5))
so <- sleuth_fit(so)
models(so)
DE_LQ4 <- sleuth_wt(so, which_beta = 'conditionLQ5')
sleuth_results_LQ4 <- sleuth_results(DE_LQ4,test = 'conditionLQ4', show_all = TRUE)
write.xlsx(sleuth_results_LQ4, file='../transcriptome_reanalysis/warner_et_al_2019_analysis/results_sleuth/sleuth_resultsLW4xLQ4.xlsx', sheetName = "DE_LQ5")


####L5 larvae
metadata <- data.frame(sample = sample_id,condition = factor(c(rep("LQ5",3),rep("LW5",3))))
paths <- file.path(base_dir, sample_id)
paths
names(paths) <- sample_id
paths
metadata$path <- paths
metadata
metadata$condition <- relevel(metadata$condition, ref = "LW5")
design <- ~ condition
so <- sleuth_prep(metadata, full_model = design, num_cores = 4L, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, transform_fun_counts = function (x) log2(x+0.5))
so <- sleuth_fit(so)
models(so)
DE_LQ5 <- sleuth_wt(so, which_beta = 'conditionLQ5')
sleuth_results_LQ5 <- sleuth_results(DE_LQ5,test = 'conditionLQ5', show_all = TRUE)
write.xlsx(sleuth_results_LQ5, file='../transcriptome_reanalysis/warner_et_al_2019_analysis/results_sleuth/sleuth_resultsLW5xLQ5.xlsx', sheetName = "DE_LQ5")

