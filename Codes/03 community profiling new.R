library(microbiome)
library(phyloseq)
library(knitr)
library(dplyr)
library(vegan)
library(ggplot2)


set.seed(123)

setwd("")



#load prefiltered data

tse.micro.BE.new1 <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds")

pseq.ConQuR_BE_remove.new1 <-readRDS("data preparation/pseq.ConQuR_BE_remove.new1.rds")

## Use microbiome package
# data transformation
pseq.ConQuR_BE_remove.new1.comptrans <- microbiome::transform(pseq.ConQuR_BE_remove.new1,
                                              "compositional")

total_samples <- phyloseq::nsamples(pseq.ConQuR_BE_remove.new1.comptrans)
#pseq.silva.comptrans.phylum <- aggregate_taxa(pseq.silva.comptrans,level = "Phylum")
pseq.ConQuR_BE_remove.new1.comptrans.phylum <- aggregate_rare(pseq.ConQuR_BE_remove.new1.comptrans,
                                             level = "Phylum",
                                             detection = 1/100,
                                             prevalence = 1/100)
pseq.ConQuR_BE_remove.new1.comptrans.genus <- aggregate_rare(pseq.ConQuR_BE_remove.new1.comptrans,
                                       level = "Genus",
                                       detection = 1/100,
                                       prevalence = 10/100)
pseq.ConQuR_BE_remove.new1.comptrans.family <- aggregate_rare(pseq.ConQuR_BE_remove.new1.comptrans,
                                             level = "Family",
                                             detection = 1/100,
                                             prevalence = 25/total_samples)
# Limit the analysis on core taxa and specific sample group
library(hrbrthemes)
library(RColorBrewer)
#library(pals)
mycolor <- colorRampPalette(brewer.pal(12, "Paired"))(35)
#mycolor <- colorRampPalette(brewer.pal(8, "Set1"))(32)
mycolor1 <- colorRampPalette(brewer.pal(12, "Paired"))(15)
p.rel.abund.all.phylum <- plot_composition(pseq.ConQuR_BE_remove.new1.comptrans.phylum,
                      taxonomic.level = "Phylum",
                      group_by = "Sample", sample.sort = "Firmicutes",
                      otu.sort = "abundance") +
  scale_fill_manual("Phylum", values = mycolor1) +
  labs(y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.position = "bottom")
p.rel.abund.all.phylum

p.rel.abund.all.phylum.avebysource <- plot_composition(pseq.ConQuR_BE_remove.new1.comptrans.phylum,
                                                       otu.sort = "abundance",
                                                       average_by = "Sample") +
  scale_fill_manual("Phylum", values = mycolor1) +
  labs(y = "Relative abundance") +
  theme_ipsum(grid= F,
              axis_title_size = 12) +
  theme(axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5))
p.rel.abund.all.phylum.avebysource
write.csv(p.rel.abund.all.phylum.avebysource$data, "Figures/composition/p.rel.abund.all.phylum.avebysource.csv")

p.rel.abund.all.genus <- plot_composition(pseq.ConQuR_BE_remove.new1.comptrans.genus,
                                          taxonomic.level = "Genus",
                                          group_by = "Sample",
                                          otu.sort = "abundance") +
  scale_fill_manual("Genus", values = mycolor) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(face = "italic"))
p.rel.abund.all.genus

p.rel.abund.all.genus.avebysource <- plot_composition(pseq.ConQuR_BE_remove.new1.comptrans.genus,
                                                      otu.sort = "abundance",
                                                      average_by = "Sample") +
  scale_fill_manual("Genus", values = mycolor) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.text = element_text(face = "italic"))

p.rel.abund.all.genus.avebysource
write.csv(p.rel.abund.all.genus.avebysource$data, "Figures/composition/p.rel.abund.all.genus.avebysource.csv")

p.rel.abund.all.family <- plot_composition(pseq.ConQuR_BE_remove.new1.comptrans.family,
                                          taxonomic.level = "Family",
                                          group_by = "Sample",
                                          otu.sort = "abundance") +
  scale_fill_manual("Family", values = mycolor1) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F,
              axis_title_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5),
        legend.position = "bottom")
p.rel.abund.all.family

p.rel.abund.all.family.avebysource <- plot_composition(pseq.ConQuR_BE_remove.new1.comptrans.family,
                                                      otu.sort = "abundance",
                                                      average_by = "Sample") +
  scale_fill_manual("Family", values = mycolor1) +
  labs(x = "Samples",
       y = "Relative abundance") +
  theme_ipsum(grid=F, 
              axis_title_size = 12) +
  theme(axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0.5))
p.rel.abund.all.family.avebysource
write.csv(p.rel.abund.all.family.avebysource$data, "Figures/composition/p.rel.abund.all.family.avebysource.csv")




