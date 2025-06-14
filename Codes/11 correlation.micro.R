library(dplyr)
library(reshape2)
library(mia)
library(psych)
library(corrplot)
library(ggpubr)
library(verification)
library(MLeval)
library(pROC)
library(patchwork)
library(ComplexHeatmap)
library(mixOmics)
library(readxl)
library(writexl)
setwd("D:/Research Data/Group Members/TFC/GT 16S/personal analysis/R/analysis")
set.seed(123)
tse.micro <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds") %>%
  transformAssay(assay.type = "counts", method = "relabundance", 
                 MARGIN = "samples") %>% 
  transformAssay(assay.type = "relabundance", method = "clr", MARGIN = "samples", pseudocount = TRUE)

taxatable <- rowData(tse.micro)[,taxonomyRanks(tse.micro)] %>% as.data.frame() %>%
  tibble::rownames_to_column("ASV")
sample_meta <- data.frame(colData(tse.micro))



####prepare data####
tse.micro.subsamples <- splitOn(tse.micro, f = "Sample", use_names = TRUE, MARGIN = 2)



#### microb correlation with picrust2 global ANCOMBC2 sig res####
PW_global <- read_xlsx("results/PICRUSt2/PW_global_desc.xlsx") %>%
  filter(q_val < 0.01) %>%
  filter(passed_ss == "TRUE") 
PW_adj_global <- read_xlsx("results/PICRUSt2/PW_adj_global_desc.xlsx") %>%
  filter(q_val < 0.01) %>%
  filter(passed_ss == "TRUE") 
KO_global <- read_xlsx("results/PICRUSt2/KO_global_desc.xlsx") %>%
  filter(q_val < 0.05) %>%
  filter(passed_ss == "TRUE")
KO_adj_global <- read_xlsx("results/PICRUSt2/KO_adj_global_desc.xlsx") %>%
  filter(q_val < 0.05) %>%
  filter(passed_ss == "TRUE")
PW.sig <- unique(c(PW_global$pathway, PW_adj_global$pathway))
KO.sig <- unique(c(KO_global$`function`, KO_adj_global$`function`))
p2PW <- read.csv("results/PICRUSt2/p2PW.csv") %>%
  tibble::column_to_rownames("X")
p2KO <- read.csv("results/PICRUSt2/p2KO.csv") %>%
  tibble::column_to_rownames("X")

p2PW.sig <- p2PW[PW.sig,] %>% t()
p2KO.sig <- p2KO[KO.sig,] %>% t()



df.family.tax <- read_xlsx("Figures/MaAsLin2/df.family.tax.xlsx") %>%
  filter(N.not.zero > 30)
df.family.adj.tax <- read_xlsx("Figures/MaAsLin2/df.family.adj.tax.xlsx") %>%
  filter(N.not.zero > 30) %>%
  filter(Family != "Mitochondria") # Remove mitochondria
df.genus.tax <- read_xlsx("Figures/MaAsLin2/df.genus.tax.xlsx") %>%
  filter(N.not.zero > 30)
df.genus.adj.tax <- read_xlsx("Figures/MaAsLin2/df.genus.adj.tax.xlsx") %>%
  filter(N.not.zero > 30)
df.sp.tax <- read_xlsx("Figures/MaAsLin2/df.sp.tax.xlsx") %>%
  filter(N.not.zero > 30)
df.sp.adj.tax <- read_xlsx("Figures/MaAsLin2/df.sp.adj.tax.xlsx") %>%
  filter(N.not.zero > 30)

family.sig <- unique(c(df.family.tax$Family, df.family.adj.tax$Family))
genus.sig <- unique(c(df.genus.tax$Genus, df.genus.adj.tax$Genus))
sp.sig <- unique(c(df.sp.tax$Species, df.sp.adj.tax$Species))

####PW correlation to family####
p2PW_desc <- rbind(PW_global, PW_adj_global) %>% 
  dplyr::select(`pathway`, `description`) %>%
  unique()

df.sig.f.clr <- t(as.data.frame(assay(altExp(tse.micro, "Family"), "clr")))[,family.sig]

corr.family <- corr.test(x = df.sig.f.clr, 
                        y = p2PW.sig,
                        use = "pairwise",method="spearman",adjust="BH",
                        alpha=.05,ci=TRUE,minlength=5,normal=TRUE)


significant_cor <- corr.family$r 
significant_cor[corr.family$p.adj> 0.05] <- NA # Remove rows and columns with less than 5 (inclusive) significant
significant_cor <- significant_cor[rowSums(!is.na(significant_cor)) > 5, colSums(!is.na(significant_cor)) > 5]
corr.family.filt <- list()
corr.family.filt$r <- corr.family$r[rownames(significant_cor), colnames(significant_cor)]
corr.family.filt$padj <- corr.family$p.adj[rownames(significant_cor), colnames(significant_cor)]

r<- corr.family.filt$r
p <- corr.family.filt$padj
cl.all1 <- rbind(df.family.tax,
                df.family.adj.tax) %>%
  distinct(Family, .keep_all =TRUE)
cl.all2 <- rbind(df.genus.tax,
                df.genus.adj.tax) %>%
  distinct(Genus, .keep_all =TRUE)
cl.all <- unique(c(cl.all1$Phylum, cl.all2$Phylum))
col2 <- setNames(c(paletteer::paletteer_d("ggsci::category20_d3", 
                                          length(cl.all))), c(cl.all))
cl <- rbind(df.family.tax,
            df.family.adj.tax) %>%
  distinct(Family, .keep_all =TRUE)
cl <- cl %>% filter(Family %in% rownames(corr.family.filt$r))

require(RColorBrewer)
col_fun = colorRampPalette(colors = rev(brewer.pal(9,"RdBu")))
# level_group <- factor(df.clinical.num.group$Group, levels = c("PD", "NPD"))
cell_fun = function(j, i, x, y, w, h, fill) {
  if(p[i, j] < 0.001 & !is.na(p[i, j])) {
    grid.text("***", x, y)
  } else if(p[i, j] < 0.01& !is.na(p[i, j])) {
    grid.text("**", x, y)
  }
  else if(p[i, j] < 0.05& !is.na(p[i, j])) {
    grid.text("*", x, y)
  }
}
ha2 = rowAnnotation(Phylum = cl$Phylum, show_annotation_name = F,
                    simple_anno_size = unit(0.3,"cm"), 
                    col = list(Phylum = col2[unique(cl$Phylum)]))
p.family_PW <- Heatmap(r, name = "Spearman's \u03C1",
                      left_annotation = ha2, 
                      cluster_rows = T, cluster_columns = T, row_km = 3, column_km = 3,
                      row_names_side = "left", column_names_side = "top",
                      show_column_dend = F, show_row_dend = F,
                      row_names_gp = gpar(fontsize = 16),
                      column_names_gp = gpar(fontsize = 14),
                      row_title_rot = 90,
                      heatmap_legend_param = list(legend_gp = gpar(fontsize = 14)),
                      cluster_column_slices =F, column_title=NULL,
                      cell_fun = cell_fun, 
                      col = col_fun(9), column_names_rot = 45, 
                      width = ncol(r)*unit(6, "mm"), 
                      height = nrow(r)*unit(6, "mm")) +
  rowAnnotation(block = anno_block(gp = gpar(fill = 2:4, col = NA),
                                   labels = c("1", "2", "3")),
                show_annotation_name = F, width = unit(0.5, "cm"))
print(p.family_PW)


corr.family.PW.r <- corr.family[["r"]] %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.family.PW.r) <- c("Family", "p2PW", "Spearman_rho")
corr.family.PW.p <- corr.family[["p"]] %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.family.PW.p) <- c("Family", "p2PW", "p_val")
corr.family.PW.p.adj <- corr.family[["p.adj"]] %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.family.PW.p.adj) <- c("Family", "p2PW", "p_adj")
corr.family.PW.summary<-cbind(corr.family.PW.r, corr.family.PW.p, corr.family.PW.p.adj) %>%
  dplyr::select(c(1:3,6,9)) %>%
  dplyr::left_join(p2PW_desc, by = c("p2PW" = "pathway")) %>%
  relocate(description, .after = "p2PW") %>%
  arrange(p_adj) %>%
  filter(p_adj < 0.05)
file.sig <- paste0("results/PICRUSt2/", 
                   paste("Family","p2PW", "corr.microb.p2.summary.csv", sep = "_"))
write.csv(corr.family.PW.summary, file.sig, row.names = F)

####PW correlation to genus####
df.sig.g.clr <- t(as.data.frame(assay(altExp(tse.micro, "Genus"), "clr")))[,genus.sig]

corr.genus <- corr.test(x = df.sig.g.clr, 
                         y = p2PW.sig,
                         use = "pairwise",method="spearman",adjust="BH",
                         alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
significant_cor <- corr.genus$r 
significant_cor[corr.genus$p.adj> 0.05] <- NA # Remove rows and columns with less than 5 (inclusive) significant
#significant_cor <- significant_cor[, colSums(!is.na(significant_cor)) > 6]
significant_cor <- significant_cor[rowSums(!is.na(significant_cor)) > 5, colSums(!is.na(significant_cor)) > 7]
corr.genus.filt <- list()
corr.genus.filt$r <- corr.genus$r[rownames(significant_cor), colnames(significant_cor)]
corr.genus.filt$padj <- corr.genus$p.adj[rownames(significant_cor), colnames(significant_cor)]

r<- corr.genus.filt$r
p <- corr.genus.filt$padj

cl <- rbind(df.genus.tax,
            df.genus.adj.tax) %>%
  distinct(Genus, .keep_all =TRUE)
cl <- cl %>% filter(Genus %in% rownames(corr.genus.filt$r))
ha2 = rowAnnotation(Phylum = cl$Phylum, show_annotation_name = F,
                    simple_anno_size = unit(0.3,"cm"), 
                    col = list(Phylum = col2[unique(cl$Phylum)]))
p.genus_PW <- Heatmap(r, name = "Spearman's \u03C1",
                       left_annotation = ha2, 
                       cluster_rows = T, cluster_columns = T, row_km = 3, column_km = 3,
                       row_names_side = "left", column_names_side = "top",
                       show_column_dend = F, show_row_dend = F,
                       row_names_gp = gpar(fontface = "italic",fontsize = 16),
                       column_names_gp = gpar(fontsize = 14),
                       row_title_rot = 90,
                       heatmap_legend_param = list(legend_gp = gpar(fontsize = 14)),
                       cluster_column_slices =F, column_title=NULL,
                       cell_fun = cell_fun, 
                       col = col_fun(9), column_names_rot = 45, 
                       width = ncol(r)*unit(7, "mm"), 
                       height = nrow(r)*unit(6, "mm")) +
  rowAnnotation(block = anno_block(gp = gpar(fill = 2:4, col = NA),
                                   labels = c("1", "2", "3")),
                show_annotation_name = F, width = unit(0.5, "cm"))
print(p.genus_PW)


corr.genus.PW.r <- corr.genus[["r"]] %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.genus.PW.r) <- c("Genus", "p2PW", "Spearman_rho")
corr.genus.PW.p <- corr.genus[["p"]] %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.genus.PW.p) <- c("Genus", "p2PW", "p_val")
corr.genus.PW.p.adj <- corr.genus[["p.adj"]] %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.genus.PW.p.adj) <- c("Genus", "p2PW", "p_adj")
corr.genus.PW.summary<-cbind(corr.genus.PW.r, corr.genus.PW.p, corr.genus.PW.p.adj) %>%
  dplyr::select(c(1:3,6,9)) %>%
  dplyr::left_join(p2PW_desc, by = c("p2PW" = "pathway")) %>%
  relocate(description, .after = "p2PW") %>%
  arrange(p_adj) %>%
  filter(p_adj < 0.05)
file.sig <- paste0("results/PICRUSt2/", 
                   paste("Genus","p2PW", "corr.microb.p2.summary.csv", sep = "_"))
write.csv(corr.genus.PW.summary, file.sig, row.names = F)


####KO correlation to family####
p2KO_desc <- rbind(KO_global, KO_adj_global) %>% 
  dplyr::select(`function`, `description`) %>%
  unique()

corr.family.KO <- corr.test(x = df.sig.f.clr, 
                             y = p2KO.sig,
                             use = "pairwise",method="spearman",adjust="BH",
                             alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corr.family.KO.r <- corr.family.KO[["r"]] %>% as.data.frame %>% 
      tibble::rownames_to_column() %>% 
      tidyr::pivot_longer(-rowname)
colnames(corr.family.KO.r) <- c("Family", "p2KO", "Spearman_rho")
corr.family.KO.p <- corr.family.KO[["p"]] %>% as.data.frame %>% 
      tibble::rownames_to_column() %>% 
      tidyr::pivot_longer(-rowname)
colnames(corr.family.KO.p) <- c("Family", "p2KO", "p_val")
corr.family.KO.p.adj <- corr.family.KO[["p.adj"]] %>% as.data.frame %>% 
      tibble::rownames_to_column() %>% 
      tidyr::pivot_longer(-rowname)
colnames(corr.family.KO.p.adj) <- c("Family", "p2KO", "p_adj")
corr.family.KO.summary<-cbind(corr.family.KO.r, corr.family.KO.p, corr.family.KO.p.adj) %>%
  dplyr::select(c(1:3,6,9)) %>%
  dplyr::left_join(p2KO_desc, by = c("p2KO" = "function")) %>%
  relocate(description, .after = "p2KO") %>%
  arrange(p_adj) %>%
  filter(p_adj < 0.05)
file.sig <- paste0("results/PICRUSt2/", 
                       paste("Family","p2KO", "corr.microb.p2.summary.csv", sep = "_"))
write.csv(corr.family.KO.summary, file.sig, row.names = F)

####KO correlation to genus####
corr.genus.KO <- corr.test(x = df.sig.g.clr, 
                        y = p2KO.sig,
                        use = "pairwise",method="spearman",adjust="BH",
                        alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
corr.genus.KO.r <- corr.genus.KO[["r"]] %>% as.data.frame %>%
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.genus.KO.r) <- c("Genus", "p2KO", "Spearman_rho")
corr.genus.KO.p <- corr.genus.KO[["p"]] %>% as.data.frame %>%
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.genus.KO.p) <- c("Genus", "p2KO", "p_val")
corr.genus.KO.p.adj <- corr.genus.KO[["p.adj"]] %>% as.data.frame %>%
  tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
colnames(corr.genus.KO.p.adj) <- c("Genus", "p2KO", "p_adj")
corr.genus.KO.summary<-cbind(corr.genus.KO.r, corr.genus.KO.p, corr.genus.KO.p.adj) %>%
  dplyr::select(c(1:3,6,9)) %>%
  dplyr::left_join(p2KO_desc, by = c("p2KO" = "function")) %>%
  relocate(description, .after = "p2KO") %>%
  arrange(p_adj) %>%
  filter(p_adj < 0.05)
file.sig <- paste0("results/PICRUSt2/", 
                   paste("Genus","p2KO", "corr.microb.p2.summary.csv", sep = "_"))
write.csv(corr.genus.KO.summary, file.sig, row.names = F)
