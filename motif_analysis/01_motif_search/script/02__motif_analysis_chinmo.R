# alorenzetti 20230408

# description ####
# in this script, we'll
# input the binding motif
# of a transcription factor (Chinmo)
# and search it upstream of
# genes of interest in the 
# honeybee genome

# setup ####
# checking if I am at the right directory
if(!dir.exists("results")){stop("I do not see the results directory necessary to run this analysis. Aborting. Set up the directory containing *results* as the working directory and try again.")}

# loading libs ####
# biocmanager is required for loading and insalling packages
if(!require("BiocManager")){install.packages("BiocManager"); library("BiocManager")}

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# choosing number of threads to use in this analysis
threads=6

# list of CRAN packages
requiredPacks = c("tidyverse",
                  "GenomicRanges",
                  "memes", "BSgenome",
                  "rtracklayer", "Biostrings",
                  "openxlsx")

# loading required cran packages
p_load(char=requiredPacks)

# checking if meme is installed
stopifnot(meme_is_installed())

# getting started ####
# creating a dir to output files
if(!dir.exists("results")){dir.create("results")}

# reading Amellifera genome
ame = readDNAStringSet(filepath = "../misc/Amellifera.fa")
ame_names_bckp = names(ame)
names(ame) = str_replace(names(ame), " .*$", "")
ame_seqinfo = Seqinfo(seqnames = names(ame), seqlengths = width(ame))

# reading annotation
annot = rtracklayer::import(con = "../misc/Amellifera.gff")
genes = annot[annot$type == "gene"]

# reading a dictionary
# of geneIDs and locus tags
# based on the annotation
dict = read_tsv(file = "data/dict.tsv")

# creating background frequency using
# 1000 random promoter sequences of amellifera
promoter_size = 1000
promoters = flank(x = genes, width = promoter_size, both = F)
seqlevels(promoters) = seqlevels(ame_seqinfo)
seqinfo(promoters) = ame_seqinfo
promoters = trim(x = promoters)
promoters_seq = getSeq(ame, promoters)
names(promoters_seq) = promoters$ID

set.seed(667)
to_sample = sample(x = 1:length(promoters), size = 1000, replace = F)

promoters_samp = promoters[to_sample]
promoters_samp_seq = getSeq(ame, promoters_samp)
writeXStringSet(x = promoters_samp_seq, filepath = "results/random_promoters_background.fa", format = "fasta")

# running fasta-get-markov
command = "/Users/alorenzetti/meme/bin/fasta-get-markov"
args = paste("-dna", "results/random_promoters_background.fa", "results/random_promoters_background_freq.txt")
system2(command = command, args = args)

# converting motif to meme
# format
command = "/Users/alorenzetti/meme/bin/chen2meme"
args = paste("-bg", "results/random_promoters_background_freq.txt", "data/TableS3Grmai.txt")
system2(command = command, args = args, stdout = "results/chinmo_motif.txt")

# running fimo for all
fimo = memes::runFimo(sequences = promoters_seq, motifs = "results/chinmo_motif.txt", outdir = "chinmo_all")
fimo = as.data.frame(fimo) %>% 
  as_tibble() %>% 
  rename(seqnames = "gffid")

# getting gene IDS for DEGs with JH vs AC
degs_hj_vs_ac = read.delim("results/hj_vs_acSig_funcat.tsv", sep = "\t") %>% 
  as_tibble()
# degs_hj_vs_ac %>% dim()

fimo_hj_vs_ac = fimo %>% 
  left_join(x = .,
            y = degs_hj_vs_ac, 
            by = "gffid") %>% 
  filter(gffid %in% degs_hj_vs_ac$gffid)
# fimo_hj_vs_ac %>% dim()

# getting gene IDS for DEGs with JH vs NT
degs_hj_vs_nt = read.delim("results/hj_vs_ntSig_funcat.tsv", sep = "\t") %>% 
  as_tibble()
# degs_hj_vs_nt %>% dim()

fimo_hj_vs_nt = fimo %>% 
  left_join(x = .,
            y = degs_hj_vs_nt, 
            by = "gffid") %>% 
  filter(gffid %in% degs_hj_vs_nt$gffid)
# fimo_hj_vs_nt %>% dim()

# writing preliminary results
write_tsv(x = left_join(x = fimo,
                        y = dict,
                        by = c("gffid" = "gff_id")), file = "results/all_promoters_with_motifs_chinmo.tsv")
write.xlsx(x = left_join(x = fimo,
                         y = dict,
                         by = c("gffid" = "gff_id")), file = "results/all_promoters_with_motifs_chinmo.xlsx")

write_tsv(x = fimo_hj_vs_ac,
          file = "results/hj_vs_ac_degs_with_motifs_chinmo.tsv")
write.xlsx(x = fimo_hj_vs_ac,
           file = "results/hj_vs_ac_degs_with_motifs_chinmo.xlsx")

write_tsv(x = fimo_hj_vs_nt,
          file = "results/hj_vs_nt_degs_with_motifs_chinmo.tsv")
write.xlsx(x = fimo_hj_vs_nt,
           file = "results/hj_vs_nt_degs_with_motifs_chinmo.xlsx")

