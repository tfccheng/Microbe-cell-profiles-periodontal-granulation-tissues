library(microbiome)
library(phyloseq)
library(ggpubr)
library(dplyr)
library(knitr)
library(eulerr)
library(RColorBrewer)

set.seed(123)

setwd("")



#load prefiltered data

tse.micro.BE.new1 <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds")

pseq.ConQuR_BE_remove.new1 <-readRDS("data preparation/pseq.ConQuR_BE_remove.new1.rds")


table(meta(pseq.ConQuR_BE_remove.new1)$Sample, meta(pseq.ConQuR_BE_remove.new1)$Sample)


#1. core genera heatmaps
#by source
Sources <- c("GT", "PT", "RT", "ST")
list_core.source <- c() # an empty object to store information
list_core.source.g <- c() 
ps.subsource.rel.list <- list()
p.core.gen.list <- list()
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)
for (n in Sources){ # for each variable n in Sample
  print(paste0("Identifying Core Taxa for ", n))
  
  ps.subsource <- subset_samples(pseq.ConQuR_BE_remove.new1, Sample == n) 
  ps.subsource.prune <- prune_taxa(taxa_sums(ps.subsource) > 0, ps.subsource)
  ps.subsource.rel <- microbiome::transform(ps.subsource.prune, "compositional")
  core_m <- core_members(ps.subsource.rel, #
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) 
  ps.subsource.rel.gen <- aggregate_taxa(ps.subsource.rel, "Genus") #plot core genus
  ps.subsource.rel.gen <- subset_taxa(ps.subsource.rel.gen, Genus!="Unknown") 
  core_m.g <- core_members(ps.subsource.rel.gen, #
                        detection = 0.001, # 0.001 in at least 90% samples 
                        prevalence = 0.5)
  p.core.gen <- plot_core(ps.subsource.rel.gen, 
                  plot.type = "heatmap", 
                  colours = rev(brewer.pal(5, "RdBu")),
                  prevalences = prevalences, 
                  detections = detections, min.prevalence = .5) +
    xlab("Detection Threshold\n(Relative Abundance (%))")+
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30, whitespace_only = F)) +
    theme_bw() + ggtitle(n) +
    ylab("Genus")
  list_core.source[[n]] <- core_m # add to a list core taxa for each group.
  list_core.source.g[[n]] <- core_m.g # add to a list core taxa for each group.
  ps.subsource.rel.list[[n]] <- ps.subsource.rel
  p.core.gen.list[[n]] <- p.core.gen #store core plot in list
  print(core_m.g)
}

#plot core genera in different sources i.e. case/control in three samples
p.core.genus.combine <- ggarrange(plotlist = p.core.gen.list, 
          ncol = 4, nrow =1, 
          legend = "right", common.legend = TRUE)
p.core.genus.combine

#Venn diagram of core taxa
ggvenn::ggvenn(list_core.source.g, 
       fill_color = c("#0f77c1", "#00a087", "#e64b35", "#631879"),
       stroke_size = 0.5, set_name_size = 5)

overlap_core.g <- intersect(list_core.source.g[[1]], intersect(list_core.source.g[[2]], 
                                                               intersect(list_core.source.g[[3]], list_core.source.g[[4]])))
setdiff(list_core.source.g[[2]], list_core.source.g[[3]])
