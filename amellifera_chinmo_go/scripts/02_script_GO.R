#Installing necessary libraries
library(AnnotationDbi)
library(ggplot2)

BiocManager::install("AnnotationForge")
library(AnnotationForge)

BiocManager::install("GO.db")
library(GO.db)

BiocManager::install("clusterProfiler")
library(clusterProfiler)

BiocManager::install("DOSE")

#Creating a local honey bee package for gene ontology annotation####


baseorg <- read.csv2("GFA_gene_reports/gene_ontology_results.csv", header = T, sep=",", dec=",")

gene_info <- data.frame(baseorg$gene_id, baseorg$symbol, baseorg$gene_name)
colnames(gene_info) <- c("GID", "SYMBOL", "GENENAME")

chromosome <-data.frame(baseorg$gene_id,baseorg$chromosome)
colnames(chromosome) <- c("GID","chromosome")

GO <- data.frame(baseorg$gene_id,baseorg$GO_ID, baseorg$Evidence_Code)
colnames(GO) <- c("GID","GO","EVIDENCE")

#removing redundant data
gene_info <- unique(gene_info)
chromosome <- unique(chromosome)
GO <- unique(GO)

anyDuplicated(gene_info)
anyDuplicated(chromosome)    
anyDuplicated(GO)           


#Creating the Org package
makeOrgPackage(
  gene_info = gene_info,
  chromosome = chromosome,
  go = GO,
  version = "1.0.0",
  maintainer = "Arthur Roque <a.roque0@usp.br>",
  author = "Arthur Roque",
  outputDir = "../amellifera_chinmo_go",
  tax_id = "7460",
  genus = "Apis",
  species = "mellifera",
  goTable = "go"
)

install.packages("../amellifera_chinmo_go/org.Amellifera.eg.db", repos = NULL, type = "source")
library(org.Amellifera.eg.db)
library(readxl)

#GO analysis for this study's transcriptome####
hj_vs_ac <- read_excel("data/hj_vs_acSig.xlsx",    col_types = c("text", "text", "text", 
"text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

genes_to_test_justino_JH_vs_AC <- as.list(hj_vs_ac$geneName[
  hj_vs_ac$padj <= 0.05 & 
    !is.na(hj_vs_ac$padj) &
    (hj_vs_ac$log2FoldChange <= -0.5 | hj_vs_ac$log2FoldChange >= 0.5)
])

GO_results_justino_BP <- enrichGO(gene = genes_to_test_justino_JH_vs_AC,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_justino_MF <- enrichGO(gene = genes_to_test_justino_JH_vs_AC,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "MF")
GO_results_justino_CC <- enrichGO(gene = genes_to_test_justino_JH_vs_AC,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "CC")

as.data.frame(GO_results_justino_BP)
as.data.frame(GO_results_justino_MF)
as.data.frame(GO_results_justino_CC)

# Saving results as a csv file.
write.csv(as.data.frame(GO_results_justino_CC), file = "tables/GO_results_justino_CC.csv", row.names = FALSE)
write.csv(as.data.frame(GO_results_justino_MF), file = "tables/GO_results_justino_MF.csv", row.names = FALSE)
write.csv(as.data.frame(GO_results_justino_BP), file = "tables/GO_results_justino_BP.csv", row.names = FALSE)

top_BP <- as.data.frame(GO_results_justino_BP) %>% head(5)
top_MF <- as.data.frame(GO_results_justino_MF) %>% head(5)
top_CC <- as.data.frame(GO_results_justino_CC) %>% head(5)

# Add an ONTOLOGY column
top_MF$ONTOLOGY <- "MF"
top_CC$ONTOLOGY <- "CC"
top_BP$ONTOLOGY <- "BP"
combined_results <- rbind(top_MF, top_CC, top_BP)

fit_justino <- ggplot(combined_results, aes(x = RichFactor, y = Description)) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(size = Count, colour = p.adjust)) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free") +
  scale_colour_gradient(low = "blue", high = "#FFCC00", guide = "colorbar") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),   
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),              
    strip.text = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Rich Factor",
    y = "GO Term",
    colour = "Adjusted p-value",
    size = "Gene Count"
  )


fit_justino




##GO for all genes with Chinmo's motif####
all_promoters_motif <- read.csv2("data/all_promoters_with_motifs_chinmo.csv", header = T, sep=",", dec=",")
all_promoters_motif_list <- as.list(unique(all_promoters_motif$gene_names))

GO_results_all_genes_BP <- enrichGO(gene = all_promoters_motif_list,pvalueCutoff = 0.5, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_all_genes_MF <- enrichGO(gene = all_promoters_motif_list,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "MF")
GO_results_all_genes_CC <- enrichGO(gene = all_promoters_motif_list,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "CC")

as.data.frame(GO_results_all_genes_BP)
#'None categories found enriched. Not writing results.
as.data.frame(GO_results_all_genes_MF)
#'None categories found enriched. Not writing results.
as.data.frame(GO_results_all_genes_CC)


#DEGs of this study with Chinmo's motif####
justino_motif <- read.csv2("data/genes_with_motif_Justino_JH_vs_AC.csv", header = T, sep=",", dec=",")
genes_to_test_justino_motif <- as.list(unique(justino_motif$gene_names))

GO_results_justino_motif_BP <- enrichGO(gene = genes_to_test_justino_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_justino_motif_MF <- enrichGO(gene = genes_to_test_justino_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "MF")
GO_results_justino_motif_CC <- enrichGO(gene = genes_to_test_justino_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "CC")

as.data.frame(GO_results_justino_motif_BP)
as.data.frame(GO_results_justino_motif_MF)
as.data.frame(GO_results_justino_motif_CC)
bp<-as.data.frame(GO_results_justino_motif_BP)
bp$GeneRatio
# Saving results as a csv file
write.csv(as.data.frame(GO_results_justino_motif_BP), file = "tables/GO_results_justino_motif_BP.csv", row.names = FALSE)
write.csv(as.data.frame(GO_results_justino_motif_MF), file = "tables/GO_results_justino_motif_MF.csv", row.names = FALSE)
write.csv(as.data.frame(GO_results_justino_motif_CC), file = "tables/GO_results_justino_motif_CC.csv", row.names = FALSE)

#
top_MF_motif <- as.data.frame(GO_results_justino_motif_MF) %>% head(5)
top_BP_motif <- as.data.frame(GO_results_justino_motif_BP) %>% head(5)
top_CC_motif <- as.data.frame(GO_results_justino_motif_CC) %>% head(5)

# Add an ONTOLOGY column
top_MF_motif$ONTOLOGY <- "MF"
top_BP_motif$ONTOLOGY <- "BP"
top_CC_motif$ONTOLOGY <- "CC"

combined_results_motif <- rbind(top_MF_motif, top_BP_motif, top_CC_motif)

fit_justino_motif <- ggplot(combined_results_motif, aes(x = RichFactor, y = Description)) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(size = Count, colour = p.adjust)) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free") +
  scale_colour_gradient(low = "blue", high = "#FFCC00", guide = "colorbar") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),   
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),                
    strip.text = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Rich Factor",
    y = "GO Term",
    colour = "Adjusted p-value",
    size = "Gene Count"
  )


fit_justino_motif


library(ggpubr)

combinedplot<- ggarrange(annotate_figure(fit_justino,top = text_grob("GO enriched terms for DEGs", color = "black",
size = 25,face = "bold", hjust = 0.5)), 
          annotate_figure(fit_justino_motif,top = text_grob("GO enriched terms for DEGs with Chinmo motif", color = "black",
size = 25,face = "bold", hjust = 0.5)),
ncol = 1,nrow=2,labels = c("A", "B"),font.label = list(size = 20, face = "bold", family = "Arial"))


ggsave(filename = "plots/GO_plot.png",
       plot = combinedplot,
       width = 4250, 
       height = 2750, 
       units = "px", 
       dpi = 300)


#DEGs of Khamis study with Chinmo's motif####
khamis_motif <- read.csv2("data/genes_with_motif_Khamis.csv", header = T, sep=",", dec=",")
genes_to_test_khamis_motif <- as.list(unique(khamis_motif$gene_names))

GO_results_khamis_motif_BP <- enrichGO(gene = genes_to_test_khamis_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_khamis_motif_MF <- enrichGO(gene = genes_to_test_khamis_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "MF")
GO_results_khamis_motif_CC <- enrichGO(gene = genes_to_test_khamis_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "CC")

as.data.frame(GO_results_khamis_motif_BP)
as.data.frame(GO_results_khamis_motif_MF)
as.data.frame(GO_results_khamis_motif_CC)

#'None categories found enriched. Not writing results.

#DEGs of Warner study with Chinmo's motif####

warner_motif <- read.csv2("data/genes_with_motif_Warner_L5.csv", header = T, sep=",", dec=",")
genes_to_test_warner_motif <- as.list(unique(warner_motif$gene_names))

GO_results_warner_motif_BP <- enrichGO(gene = genes_to_test_warner_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "BP")
GO_results_warner_motif_MF <- enrichGO(gene = genes_to_test_warner_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "MF")
GO_results_warner_motif_CC <- enrichGO(gene = genes_to_test_warner_motif,pvalueCutoff = 0.05, pAdjustMethod = "BH", OrgDb = "org.Amellifera.eg.db", keyType = "SYMBOL", ont = "CC")

as.data.frame(GO_results_khamis_motif_BP)
#'None categories enriched Not writing results.
as.data.frame(GO_results_khamis_motif_MF)
#'None categories enriched Not writing results.
as.data.frame(GO_results_khamis_motif_CC)
#'None categories enriched Not writing results.

#Plot common genes#####
read.csv()
# Loading tables
df_justino <- read_csv("data/genes_with_motif_Justino_JH_vs_AC.csv", locale = locale(encoding = "latin1"))
df_khamis <- read_csv("data/genes_with_motif_Khamis.csv", locale = locale(encoding = "latin1"))
df_warner <- read_csv("data/genes_with_motif_Warner_L5.csv", locale = locale(encoding = "latin1"))

# Extracting gene_names
genes_justino <- unique(df_justino$gene_names[!is.na(df_justino$gene_names)])
genes_khamis <- unique(df_khamis$gene_names[!is.na(df_khamis$gene_names)])
genes_warner <- unique(df_warner$gene_names[!is.na(df_warner$gene_names)])

library('VennDiagram')
library(tidyverse)
library(gridExtra)
venn.plot <- venn.diagram(
  x = list(
    `\nThis study` = genes_justino,
    `\nKhamis et\n al. (2015)` = genes_khamis,
    `Warner et al. (2019)` = genes_warner
  ),
  filename = NULL,
  fill = c("#99CCFF", "#FFFF99", "#6666FF"), 
  alpha = 0.8,
  label.col = "black",
  cex = 2,
  cat.cex = 1.5, 
  cat.fontface = "bold"
)

grid.draw(venn.plot)


venn_recorded <- recordPlot()

png("plots/venn.png", width = 4500, height = 4500, units = "px", res = 600)

replayPlot(venn_recorded)

dev.off()

