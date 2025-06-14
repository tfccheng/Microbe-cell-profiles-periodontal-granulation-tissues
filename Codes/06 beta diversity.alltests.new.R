library(mia)
library(miaViz)
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(scater)
library(ggsignif)
library(ggpubr)
library(scales)
library(cowplot)
set.seed(123)

setwd("")


#load prefiltered data

tse.micro.BE.new1 <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds")


tse.micro.BE.new1 <- tse.micro.BE.new1 %>%
  subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05) %>% #filter with 0.05% detection & 5% prevalence
  transformAssay(assay.type = "counts", method = "relabundance") %>%
  transformAssay(assay.type = "counts", method = "log10", pseudocount = T) %>%
  transformAssay(assay.type = "counts", method = "z", MARGIN = "features") %>%
  transformAssay(assay.type = "relabundance", method = "clr", pseudocount = T) %>%
  transformAssay(assay.type = "relabundance", method = "rclr")  #transform counts to robust clr

  


#### Perform PCoA with Bray Curtis dissimilarity####

tse.micro.BE.new1.bray <- runMDS(tse.micro.BE.new1, FUN = vegan::vegdist, 
                        method = "bray", name = "PCoA_BC", exprs_values = "relabundance")

p.pcoa.bray <- plotReducedDim(tse.micro.BE.new1.bray, "PCoA_BC", 
                              colour_by = "Sample", 
                              point_size =5)
e <- attr(reducedDim(tse.micro.BE.new1.bray, "PCoA_BC"), "eig")
rel_eig <- e/sum(e[e>0])  
permanova.pcoa.bray <- adonis2(t(assay(tse.micro.BE.new1.bray,"relabundance")) ~ Sample,
                               by = "margin", 
                               data = colData(tse.micro.BE.new1.bray),
                               method = "bray",
                               permutations = 999)
write.csv(permanova.pcoa.bray, "Figures/beta/permanova.pcoa.bray.csv")
permanova.pcoa.bray.pairwise <- pairwiseAdonis::pairwise.adonis2(t(assay(tse.micro.BE.new1.bray,"relabundance")) ~ Sample, 
                                      data = colData(tse.micro.BE.new1.bray))
writexl::write_xlsx(lapply(permanova.pcoa.bray.pairwise, as.data.frame),"Figures/beta/permanova.pcoa.bray.pairwise.xlsx")

p.pcoa.bray <- p.pcoa.bray + 
  labs(subtitle = paste0("Bray-Curtis; P: ", permanova.pcoa.bray$`Pr(>F)`),
       x = paste("PCoA1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
       y = paste("PCoA2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))+
  stat_ellipse(aes(color = colour_by)) +
  scale_color_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879"))+
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw()+
  theme(aspect.ratio=1) 
p.pcoa.bray





tse.micro.BE.new1.bray.nmds <- runNMDS(tse.micro.BE.new1, FUN = vegan::vegdist, 
                                       method = "bray", name = "NMDS_BC", 
                                       exprs_values = "relabundance")
p.pcoa.bray.nmds <- plotReducedDim(tse.micro.BE.new1.bray.nmds, "NMDS_BC", 
                                   colour_by = "Sample", 
                                   point_size =5)+ 
  labs(subtitle = paste0("Bray-Curtis; P: ", permanova.pcoa.bray$`Pr(>F)`),
       x = "NMDS1", y = "NMDS2")+
  stat_ellipse(aes(color = colour_by)) +
  scale_color_manual(name = "Sample",values = c("#0f77c1", "#00a087", "#e64b35", "#631879"))+
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw()+
  theme(aspect.ratio=1) 
p.pcoa.bray.nmds

ggarrange(p.pcoa.bray, p.pcoa.bray.nmds, ncol = 2)

####Unifrac####
tse.micro.BE.new1.unifrac_w <- runMDS(tse.micro.BE.new1, FUN = mia::calculateUnifrac, 
                        name = "PCoA_Unifrac", weighted = TRUE,
                        tree = rowTree(tse.micro.BE.new1),
                        ntop = nrow(tse.micro.BE.new1),
                        exprs_values = "counts")
p.pcoa.unifrac_w <- plotReducedDim(tse.micro.BE.new1.unifrac_w, "PCoA_Unifrac", 
                                       colour_by = "Sample", shape_by = "Sample",
                                       point_size =5)

e <- attr(reducedDim(tse.micro.BE.new1.unifrac_w, "PCoA_Unifrac"), "eig")
rel_eig <- e/sum(e[e>0])  
permanova.pcoa.unifrac <- adonis2(t(assay(tse.micro.BE.new1.unifrac_w,"counts")) ~ Sample,
                                        by = "margin", 
                                        data = colData(tse.micro.BE.new1.unifrac_w),
                                        permutations = 999)
pairwiseAdonis::pairwise.adonis2(t(assay(tse.micro.BE.new1.unifrac_w,"counts")) ~ Sample, 
                 data = colData(tse.micro.BE.new1.unifrac_w))

p.pcoa.unifrac_w <- p.pcoa.unifrac_w + 
  labs(subtitle = paste0("Weighted UniFrac; P: ", permanova.pcoa.unifrac$`Pr(>F)`),
       x = paste("PCoA1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
       y = paste("PCoA2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))+
  stat_ellipse(aes(color = colour_by)) +
  scale_color_manual(values = c("#0f77c1", "#00a087", "#e64b35", "#631879"))+
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw()+
  theme(aspect.ratio=1) 
p.pcoa.unifrac_w

tse.micro.BE.new1.unifrac_w <- runNMDS(tse.micro.BE.new1, FUN = mia::calculateUnifrac, 
                                       name = "NMDS_BC", weighted = TRUE,
                                       tree = rowTree(tse.micro.BE.new1),
                                       ntop = nrow(tse.micro.BE.new1),
                                       exprs_values = "counts")
p.pcoa.unifrac_w <- plotReducedDim(tse.micro.BE.new1.unifrac_w, "NMDS_BC", 
                                   colour_by = "Sample", shape_by = "Sample",
                                   point_size =5)





####Perform PCoA with robust Aitchison distance####


tse.micro.BE.new1 <- runMDS(tse.micro.BE.new1, FUN = vegan::vegdist, 
              method = "euclidean", name = "PCoA_Aitchison", assay.type = "rclr")
p.pcoa.aitchison.bysample <- plotReducedDim(tse.micro.BE.new1, "PCoA_Aitchison", 
                                            colour_by = "Sample", shape_by = "Sample",
                                            point_size =5)
e <- attr(reducedDim(tse.micro.BE.new1, "PCoA_Aitchison"), "eig")
rel_eig <- e/sum(e[e>0])  
permanova.pcoa.aitchison.bysample <- adonis2(t(assay(tse.micro.BE.new1,"rclr")) ~ Sample,
                                             by = "margin", 
                                             data = colData(tse.micro.BE.new1),
                                             method = "robust.aitchison",
                                             permutations = 999)
p.pcoa.aitchison.bysample <- p.pcoa.aitchison.bysample + 
  labs(subtitle = paste0("Aitchison; P: ", permanova.pcoa.aitchison.bysample$`Pr(>F)`),
       x = paste("PCoA1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
       y = paste("PCoA2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))+
  scale_color_manual(values = c("#0f77c1", "#00a087", "#e64b35"))+
  scale_x_continuous(labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_bw()+
  theme(aspect.ratio=1) 
p.pcoa.aitchison.bysample


tse.micro.BE.new1 <- runTSNE(tse.micro.BE.new1, name = "tSNE", assay.type = "log10")
plotReducedDim(tse.micro.BE.new1, "tSNE", 
               colour_by = "Sample", shape_by = "Sample",
               point_size =5)
tse.micro.BE.new1 <- runUMAP(tse.micro.BE.new1, name = "UMAP", assay.type = "log10")
plotReducedDim(tse.micro.BE.new1, "UMAP", 
               colour_by = "Sample", shape_by = "Sample",
               point_size =5)
