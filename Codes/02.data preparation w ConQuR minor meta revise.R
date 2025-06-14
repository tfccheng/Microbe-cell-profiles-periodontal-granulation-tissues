library(microbiome)
library(phyloseq)
library(knitr)
library(dplyr)
library(readxl)


setwd("")

set.seed(123)
#### prepare data####
asv_tab <- read.csv("ASVs_counts.txt", sep = "\t", row.names = 1) %>% as.matrix()
tax_tab.silva <- read.csv("ASV_silva_species_combined_genus.txt", sep = "\t", row.names = 1) %>% as.matrix()
colnames(tax_tab.silva)[1]  <- "Domain"
meta_tab <- read.csv("meta.csv")
colnames(meta_tab)[1] <- "seq_sample_name"
meta_tab$Sample <- as.factor(meta_tab$Sample)
meta_tab$Batch <- as.factor(meta_tab$Batch)
clinic_meta <- read_xlsx("(DL) Subject questionnaire + baseline data 31.10.2023.xlsx")[-c(1,50:52),]
#clinic_meta <- read_xlsx("subject clinical data  for R.xlsx")
clinic_meta.require <- clinic_meta %>%
  select(Case, Sex, Age, BMI, `BMI:\r\nNormal = 1\r\nOverweight/Obese = 2\r\nWHO`,
         `Diagnosis code (stage)`, `Full Diagnosis`, Extent, Grade)
colnames(clinic_meta.require)[5:7] <- c("Obese", "Stage", "Full_Diagnosis")

meta_tab.final <- meta_tab %>%
  dplyr::left_join(clinic_meta.require, by = "Case" ) %>%
  tibble::column_to_rownames("seq_sample_name")
meta_tab.final$Sex <- factor(meta_tab.final$Sex)
meta_tab.final$Stage <- factor(meta_tab.final$Stage)
meta_tab.final$Grade <- factor(meta_tab.final$Grade)
meta_tab.final <- meta_tab.final %>% tibble::column_to_rownames("seq_sample_name")
meta_tab.final$BMI <- as.numeric(meta_tab.final$BMI)
write.csv(meta_tab.final, "meta.final.csv")
meta_tab.finalest <- read.csv("meta.finalest.csv", row.names = 1) %>%
  mutate(across(c(where(is.character), -c(Case, True.Name)), as.factor))
meta_tab.finalest$Batch <- as.factor(meta_tab.finalest$Batch)
saveRDS(meta_tab.finalest, "data preparation/meta_tab.finalest.rds")

phy_tre <- read_tree("GT_nhx.tre")
ref_seq <- Biostrings::readDNAStringSet("ASVs.fa", format="fasta")
pseq.silva <- phyloseq(otu_table(asv_tab, taxa_are_rows = T),
                       tax_table(tax_tab.silva),
                       sample_data(meta_tab.finalest),
                       phy_tree(phy_tre),
                       refseq(ref_seq))
saveRDS(pseq.silva, "data preparation/pseq.silva.rds")
pseq.silva.minimal.filter <- filter_taxa(pseq.silva, function (x) {sum(x > 0) > 1}, prune=TRUE) 
saveRDS(pseq.silva.minimal.filter, "data preparation/pseq.silva.minimal.filter.rds")
#pseq.silva.minimal.filter <- readRDS("pseq.silva.minimal.filter.rds")

pseq.silva.minimal.filter.f <- microbiome::add_besthit(pseq.silva.minimal.filter)
taxa_names(pseq.silva.minimal.filter.f) <- 
  gsub(":.*\\.", ":", taxa_names(pseq.silva.minimal.filter.f)) #clean names
saveRDS(pseq.silva.minimal.filter.f, "data preparation/pseq.silva.minimal.filter.f.rds")

#pseq.silva.minimal.filter.f <- readRDS("pseq.silva.minimal.filter.f.rds")

#### remove batch effects####
library(ConQuR)
library(doParallel)
##### convert phyloseq object to  TreeSummarizedExperiment object#####
library(mia)
tse.silva <- makeTreeSummarizedExperimentFromPhyloseq(pseq.silva.minimal.filter.f)
saveRDS(tse.silva, "data preparation/tse.silva.rds")


taxa.filt <- as.matrix(assay(tse.silva)) %>% t()
saveRDS(taxa.filt, "data preparation/taxa.filt.rds")


batchid <- meta_tab.finalest$Batch
summary(batchid)
covar <- meta_tab.finalest[, c('Sample', 'Sex', 'Age', 'BMI', 'Obese', 'Full_Diagnosis', 'Smoker', 'DM')]

summary(covar)


####fine tune ConQuR on HPC
result_tuned.new1 = Tune_ConQuR(tax_tab=taxa.filt, batchid=batchid, covariates=covar,
                           batch_ref_pool=c("1", "2", "3", "4", "5"),
                           logistic_lasso_pool= c(T,F), 
                           quantile_type_pool=c("standard", "lasso", "composite"),
                           simple_match_pool=c(T,F),
                           lambda_quantile_pool=c(NA, "2p/n", "2p/logn"),
                           interplt_pool=c(T,F),
                           frequencyL=0,
                           frequencyU=1,
                           num_core = 64)
result_tuned.new1 <- readRDS("result_tuned.new1.rds")

result_tuned.new1$method_final
taxa_optimal.new1 <- result_tuned.new1$tax_final


#Check how ConQuR fine tune works
par(mfrow=c(1, 2))

Plot_PCoA(TAX=taxa.filt, factor=batchid, main="Before Correction, Bray-Curtis")
Plot_PCoA(TAX=taxa_optimal.new1, factor=batchid, main="Fine-Tuned ConQuR, Bray-Curtis")



Plot_PCoA(TAX=taxa.filt, factor=batchid, dissimilarity="Aitch", main="Before Correction, Aitchison")

Plot_PCoA(TAX=taxa_optimal.new1, factor=batchid, dissimilarity="Aitch", main="Fine-Tuned ConQuR, Aitchison")


tse.silva.subsample <- tse.silva 
for (m in c("Family", "Genus", "Species")){
  altExp(tse.silva.subsample, m) <- agglomerateByRank(tse.silva.subsample, rank = m, onRankOnly = T)
}


tse.silva.ConQuR_BE_remove.new1 <- tse.silva
assay(tse.silva.ConQuR_BE_remove.new1) <- taxa_optimal.new1 %>% t()
for (m in c("Family", "Genus", "Species")){
  altExp(tse.silva.ConQuR_BE_remove.new1, m) <- agglomerateByRank(tse.silva.ConQuR_BE_remove.new1, rank = m, onRankOnly = T)
}


names(rowData(tse.silva.subsample))[1] <- "Kingdom"
saveRDS(tse.silva.subsample, "data preparation/tse.silva.subsample.rds")

names(rowData(tse.silva.ConQuR_BE_remove.new1))[1] <- "Kingdom"
saveRDS(tse.silva.ConQuR_BE_remove.new1, "data preparation/tse.silva.ConQuR_BE_remove.new1.rds")

pseq.subsample <- makePhyloseqFromTreeSummarizedExperiment(tse.silva.subsample)

pseq.ConQuR_BE_remove.new1 <- makePhyloseqFromTreeSummarizedExperiment(tse.silva.ConQuR_BE_remove.new1)

saveRDS(pseq.subsample, "data preparation/pseq.subsample.rds")

saveRDS(pseq.ConQuR_BE_remove.new1, "data preparation/pseq.ConQuR_BE_remove.new1.rds")

