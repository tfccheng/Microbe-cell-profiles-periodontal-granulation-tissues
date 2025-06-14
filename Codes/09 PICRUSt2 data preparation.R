library(microbiome)
library(phyloseq)
library(biomformat)
library(dplyr)

setwd("")

set.seed(123)
#load prefiltered data
pseq.micro.BE.new1 <- readRDS("data preparation/pseq.ConQuR_BE_remove.new1.rds")
a <- read.csv("results/PICRUSt2/asv_tab.txt", sep = "\t" )
# b<-data.frame(taxa_names(pseq.micro.BE.new1),a$X)
# cc<- stringr::str_split_fixed(b$`taxa_names.pseq.micro.BE.new1.`, pattern = ":", n=2)
taxa_names(pseq.micro.BE.new1) <- a$X
pseq.micro.BE.new1 %>%
  refseq() %>%
  Biostrings::writeXStringSet("results/PICRUSt2/micro.BE.new1.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")
asv_tab_for_picrust2.micro.BE.new1 <- as.data.frame(otu_table(pseq.micro.BE.new1))
write.table(asv_tab_for_picrust2.micro.BE.new1, 
            "results/PICRUSt2/asv_tab_for_picrust2.micro.BE.new1.txt", sep = "\t", quote = F, col.names = NA)




