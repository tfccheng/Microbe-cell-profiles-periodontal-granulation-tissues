library("data.table")   
library("phyloseq")
library("tidyverse")
library(readxl)
library(dplyr)
library(TreeSummarizedExperiment)
library(mia)
library(ANCOMBC)
library(KEGGREST)
library(writexl)
library(MicrobiomeProfiler)
library(clusterProfiler)
library(ggpubr)

setwd("")

p2_KO = "KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"
p2_PW = "pathways_out/path_abun_unstrat_descrip.tsv.gz"

p2KO_des = as.data.frame(fread(p2_KO)) %>%
  tibble::column_to_rownames("function")

p2KO = as.matrix(p2KO_des[,-1]) %>% round() %>% as.data.frame()
write.csv(p2KO, "D:/Research Data/Group Members/TFC/GT 16S/personal analysis/R/analysis/results/PICRUSt2/p2KO.csv")

p2PW_des = as.data.frame(fread(p2_PW)) %>%
  tibble::column_to_rownames("pathway")

p2PW = as.matrix(p2PW_des[,-1]) %>% round() %>% as.data.frame()
write.csv(p2PW, "D:/Research Data/Group Members/TFC/GT 16S/personal analysis/R/analysis/results/PICRUSt2/p2PW.csv")

setwd("D:/Research Data/Group Members/TFC/GT 16S/personal analysis/R/analysis/")
meta <- readRDS("data preparation/meta_tab.finalest.rds") %>% 
  tibble::rownames_to_column("sample_name") %>%
  dplyr::arrange(sample_name)

p2.PW <- TreeSummarizedExperiment(assays = S4Vectors::SimpleList(counts = as.matrix(p2PW)),
                                  colData = meta)
p2.KO <- TreeSummarizedExperiment(assays = S4Vectors::SimpleList(counts = as.matrix(p2KO)),
                                  colData = meta)
output.PW = ancombc2(data = p2.PW, assay_name = "counts", 
                  fix_formula = "Sample", 
                  p_adj_method = "fdr", 
                  pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                  group = "Sample", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.2, n_cl = 12, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                  iter_control = list(tol = 1e-2, max_iter = 20,
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 20),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(-1, 0, 0, 0,
                                                                1, -1, 0, 0,
                                                                0, 1, -1, 0,
                                                                0, 0, 1, -1),
                                                              nrow = 4,
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))
View(output.PW$res)
output.KO = ancombc2(data = p2.KO, assay_name = "counts", 
                     fix_formula = "Sample", 
                     p_adj_method = "fdr", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                     group = "Sample", struc_zero = TRUE, neg_lb = TRUE,
                     alpha = 0.2, n_cl = 12, verbose = TRUE,
                     global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                     iter_control = list(tol = 1e-2, max_iter = 20,
                                         verbose = TRUE),
                     em_control = list(tol = 1e-5, max_iter = 20),
                     lme_control = lme4::lmerControl(),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                     trend_control = list(contrast = list(matrix(c(-1, 0, 0, 0,
                                                                   1, -1, 0, 0,
                                                                   0, 1, -1, 0,
                                                                   0, 0, 1, -1),
                                                                 nrow = 4,
                                                                 byrow = TRUE)),
                                          node = list(2, 2),
                                          solver = "ECOS",
                                          B = 100))
View(output.KO$res)

output.PW.adj = ancombc2(data = p2.PW, assay_name = "counts", 
                         fix_formula = "Sample + Sex + Age + Obese + Smoker + DM",  
                         p_adj_method = "fdr", pseudo_sens = TRUE,
                         prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                         group = "Sample", struc_zero = TRUE, neg_lb = TRUE,
                         alpha = 0.2, n_cl = 12, verbose = TRUE,
                         global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                         iter_control = list(tol = 1e-2, max_iter = 20,
                                             verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 20),
                         lme_control = lme4::lmerControl(),
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                         trend_control = list(contrast = list(matrix(c(-1, 0, 0, 0,
                                                                       1, -1, 0, 0,
                                                                       0, 1, -1, 0,
                                                                       0, 0, 1, -1),
                                                                     nrow = 4,
                                                                     byrow = TRUE)),
                                              node = list(2, 2),
                                              solver = "ECOS",
                                              B = 100))
View(output.PW.adj$res)
output.KO.adj = ancombc2(data = p2.KO, assay_name = "counts", 
                         fix_formula = "Sample + Sex + Age + Obese + Smoker + DM", 
                         p_adj_method = "fdr", pseudo_sens = TRUE,
                         prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                         group = "Sample", struc_zero = TRUE, neg_lb = TRUE,
                         alpha = 0.2, n_cl = 12, verbose = TRUE,
                         global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = F,
                         iter_control = list(tol = 1e-2, max_iter = 20,
                                             verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 20),
                         lme_control = lme4::lmerControl(),
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                         trend_control = list(contrast = list(matrix(c(-1, 0, 0, 0,
                                                                       1, -1, 0, 0,
                                                                       0, 1, -1, 0,
                                                                       0, 0, 1, -1),
                                                                     nrow = 4,
                                                                     byrow = TRUE)),
                                              node = list(2, 2),
                                              solver = "ECOS",
                                              B = 100))
View(output.KO.adj$res)


####plot####
#####PW####
#####Global####
PW_global = output.PW$res_global
PW_res = output.PW$res
df_Sample = PW_res %>%
  dplyr::select(taxon, contains("Sample")) 
df_fig_PW_global = df_Sample %>%
  dplyr::left_join(PW_global %>%
                     dplyr::transmute(taxon, 
                                      diff_Sample = diff_abn, 
                                      passed_ss = passed_ss)) %>%
  dplyr::filter(diff_Sample == 1) %>%
  dplyr::mutate(lfc_PT = lfc_SamplePT,
                lfc_RT = lfc_SampleRT,
                lfc_ST = lfc_SampleST,
                color = ifelse(passed_ss == 1, "aquamarine3", "black")) %>%
  dplyr::transmute(taxon,
                   `PT - GT` = round(lfc_PT, 2),
                   `RT - GT` = round(lfc_RT, 2), 
                   `ST - GT` = round(lfc_ST, 2), 
                   color = color) %>%
  tidyr::pivot_longer(cols = `PT - GT`:`ST - GT`, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_PW_global$group = factor(df_fig_PW_global$group, 
                             levels = c("PT - GT",
                                        "RT - GT",
                                        "ST - GT"))

lo = floor(min(df_fig_PW_global$value))
up = ceiling(max(df_fig_PW_global$value))
mid = 0
fig_PW_global = df_fig_PW_global %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color = df_fig_PW_global %>%
                                     dplyr::distinct(taxon, color) %>%
                                     .$color))
fig_PW_global

PW_global %>%
  dplyr::filter(diff_abn == 1) %>%
  dplyr::rename(pathway = taxon) %>%
  dplyr::left_join(p2PW_des %>% 
                     tibble::rownames_to_column("pathway") %>%
                     dplyr::select(pathway, description), 
                   by = "pathway") %>%
  relocate(description, .after = pathway) %>%
  write_xlsx("results/PICRUSt2/PW_global_desc.xlsx")

PW_adj_global <- output.PW.adj$res_global
PW_adj_global %>%
  dplyr::filter(diff_abn == 1) %>%
  dplyr::rename(pathway = taxon) %>%
  dplyr::left_join(p2PW_des %>% 
                     tibble::rownames_to_column("pathway") %>%
                     dplyr::select(pathway, description), 
                   by = "pathway") %>%
  relocate(description, .after = pathway) %>%
  write_xlsx("results/PICRUSt2/PW_adj_global_desc.xlsx")
#####Dunn####
PW_dunn = output.PW$res_dunn

df_fig_PW_dunn1 = PW_dunn %>%
  dplyr::filter(diff_SamplePT == 1 | 
                  diff_SampleRT == 1 |
                  diff_SampleST == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_SamplePT == 1, 
                              round(lfc_SamplePT, 2), 0),
                lfc2 = ifelse(diff_SampleRT == 1, 
                              round(lfc_SampleRT, 2), 0),
                lfc3 = ifelse(diff_SampleST == 1, 
                              round(lfc_SampleST, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_PW_dunn2 = PW_dunn %>%
  dplyr::filter(diff_SamplePT == 1 | 
                  diff_SampleRT == 1 |
                  diff_SampleST == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_SamplePT == 1 & diff_SamplePT == 1, 
                              "#631879", "black"),
                lfc2 = ifelse(passed_ss_SampleRT == 1 & diff_SampleRT == 1, 
                              "#631879", "black"),
                lfc3 = ifelse(passed_ss_SampleST == 1 & diff_SampleST == 1, 
                              "#631879", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_PW_dunn = df_fig_PW_dunn1 %>%
  dplyr::left_join(df_fig_PW_dunn2, by = c("taxon", "group"))

df_fig_PW_dunn$group = recode(df_fig_PW_dunn$group, 
                           `lfc1` = "PT - GT",
                           `lfc2` = "RT - GT",
                           `lfc3` = "ST - GT")
df_fig_PW_dunn$group = factor(df_fig_PW_dunn$group, 
                           levels = c("PT - GT",
                                      "RT - GT",
                                      "ST - GT"))

lo = floor(min(df_fig_PW_dunn$value))
up = ceiling(max(df_fig_PW_dunn$value))
mid = 0
fig_PW_dunn = df_fig_PW_dunn %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#00468B", high = "#AD002A", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "lfc") +
  geom_text(aes(group, taxon, label = value, color = color), size = 6) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to GT\n(unadjusted)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))
fig_PW_dunn
PW_dunn_desc <- PW_dunn %>%
  dplyr::filter(diff_SamplePT == 1 | 
                  diff_SampleRT == 1 |
                  diff_SampleST == 1) %>%
  dplyr::rename(pathway = taxon) %>%
  dplyr::left_join(p2PW_des %>% 
                     tibble::rownames_to_column("pathway") %>%
                     dplyr::select(pathway, description), 
                   by = "pathway") %>%
  relocate(description, .after = pathway)
write_xlsx(PW_dunn_desc, "results/PICRUSt2/PW_dunn_desc.xlsx")
#####Adj####
PW_dunn.adj = output.PW.adj$res_dunn

df_fig_PW_dunn1.adj = PW_dunn.adj %>%
  dplyr::filter(diff_SamplePT == 1 | 
                  diff_SampleRT == 1 |
                  diff_SampleST == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_SamplePT == 1, 
                              round(lfc_SamplePT, 2), 0),
                lfc2 = ifelse(diff_SampleRT == 1, 
                              round(lfc_SampleRT, 2), 0),
                lfc3 = ifelse(diff_SampleST == 1, 
                              round(lfc_SampleST, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_PW_dunn2.adj = PW_dunn.adj %>%
  dplyr::filter(diff_SamplePT == 1 | 
                  diff_SampleRT == 1 |
                  diff_SampleST == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_SamplePT == 1 & diff_SamplePT == 1, 
                              "#631879", "black"),
                lfc2 = ifelse(passed_ss_SampleRT == 1 & diff_SampleRT == 1, 
                              "#631879", "black"),
                lfc3 = ifelse(passed_ss_SampleST == 1 & diff_SampleST == 1, 
                              "#631879", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_PW_dunn.adj = df_fig_PW_dunn1.adj %>%
  dplyr::left_join(df_fig_PW_dunn2.adj, by = c("taxon", "group"))

df_fig_PW_dunn.adj$group = recode(df_fig_PW_dunn.adj$group, 
                              `lfc1` = "PT - GT",
                              `lfc2` = "RT - GT",
                              `lfc3` = "ST - GT")
df_fig_PW_dunn.adj$group = factor(df_fig_PW_dunn.adj$group, 
                              levels = c("PT - GT",
                                         "RT - GT",
                                         "ST - GT"))

lo = floor(min(df_fig_PW_dunn.adj$value))
up = ceiling(max(df_fig_PW_dunn.adj$value))
mid = 0
fig_PW_dunn.adj = df_fig_PW_dunn.adj %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#00468B", high = "#AD002A", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "lfc") +
  geom_text(aes(group, taxon, label = value, color = color), size = 6) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to GT\n(covariate adjusted)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))
fig_PW_dunn.adj

PW_dunn.adj_desc <- PW_dunn.adj %>%
  dplyr::filter(diff_SamplePT == 1 | 
                  diff_SampleRT == 1 |
                  diff_SampleST == 1) %>%
  dplyr::rename(pathway = taxon) %>%
  dplyr::left_join(p2PW_des %>% 
                     tibble::rownames_to_column("pathway") %>%
                     dplyr::select(pathway, description), 
                   by = "pathway") %>%
  relocate(description, .after = pathway)
write_xlsx(PW_dunn.adj_desc, "results/PICRUSt2/PW_dunn.adj_desc.xlsx")

ggarrange(fig_PW_dunn, fig_PW_dunn.adj, common.legend = T)

####KEGG####
#####Global####
KO_global = output.KO$res_global
KO_global %>%
  dplyr::filter(diff_abn == 1) %>%
  dplyr::rename(`function` = taxon) %>%
  dplyr::left_join(p2KO_des %>% 
                     tibble::rownames_to_column("function") %>%
                     dplyr::select(`function`, description), 
                   by = "function") %>%
  relocate(description, .after = `function`) %>%
  write_xlsx("results/PICRUSt2/KO_global_desc.xlsx")

KO_adj_global <- output.KO.adj$res_global
KO_adj_global %>%
  dplyr::filter(diff_abn == 1) %>%
  dplyr::rename(`function` = taxon) %>%
  dplyr::left_join(p2KO_des %>% 
                     tibble::rownames_to_column("function") %>%
                     dplyr::select(`function`, description), 
                   by = "function") %>%
  relocate(description, .after = `function`) %>%
  write_xlsx("results/PICRUSt2/KO_adj_global_desc.xlsx")

#Dunn
KO_dunn_desc <- output.KO$res_dunn %>%
  # dplyr::filter(diff_SamplePT == 1 | 
  #                 diff_SampleRT == 1 |
  #                 diff_SampleST == 1) %>%
  dplyr::rename(`function` = taxon) %>%
  dplyr::left_join(p2KO_des %>% 
                     tibble::rownames_to_column("function") %>%
                     dplyr::select(`function`, description), 
                   by = "function") %>%
  relocate(description, .after = `function`)
write_xlsx(KO_dunn_desc, "results/PICRUSt2/KO_dunn_desc.xlsx")
KO_dunn.adj_desc <- output.KO.adj$res_dunn %>%
  # dplyr::filter(diff_SamplePT == 1 | 
  #                 diff_SampleRT == 1 |
  #                 diff_SampleST == 1) %>%
  dplyr::rename(`function` = taxon) %>%
  dplyr::left_join(p2KO_des %>% 
                     tibble::rownames_to_column("function") %>%
                     dplyr::select(`function`, description), 
                   by = "function") %>%
  relocate(description, .after = `function`)
write_xlsx(KO_dunn.adj_desc, "results/PICRUSt2/KO_dunn.adj_desc.xlsx")

####KEGG####
library(MicrobiomeProfiler)
library(clusterProfiler)

KO_enrich_global <- enrichKO(dplyr::union((KO_global %>% 
                                           dplyr::filter(q_val < 0.2) %>%
                                           dplyr::filter(passed_ss == 1))$taxon, 
                                          #only those passed ss were used for enrichment
                                        (KO_adj_global %>% 
                                           dplyr::filter(q_val < 0.2) %>%
                                           dplyr::filter(passed_ss == 1))$taxon), 
                           qvalueCutoff = 0.05)
KO_enrich_global.sig <- KO_enrich_global@result %>%
  dplyr::filter(qvalue < 0.05)
write_xlsx(KO_enrich_global.sig, "results/PICRUSt2/KO_enrich_global_sig.xlsx")
dp_KO <- dotplot(KO_enrich_global, showCategory = 15) + ggtitle("All - GT (Global)")
KO_enrich_Dunn.PT <- enrichKO(dplyr::union((KO_dunn_desc %>%
                                dplyr::filter(q_SamplePT < 0.2) %>%
                                  dplyr::filter(passed_ss_SamplePT == 1))$`function` ,
                                #only those passed ss were used for enrichment
                               (KO_dunn.adj_desc %>%
                                dplyr::filter(q_SamplePT < 0.2) %>%
                                  dplyr::filter(passed_ss_SamplePT == 1))$`function`) ,
                           qvalueCutoff = 0.05)
dp_KO.PT <-dotplot(KO_enrich_Dunn.PT, showCategory = 15) + ggtitle("PT - GT (Dunn)")
KO_enrich_Dunn.RT <- enrichKO(dplyr::union((KO_dunn_desc %>%
                                              dplyr::filter(q_SampleRT < 0.2) %>%
                                              dplyr::filter(passed_ss_SampleRT == 1))$`function` ,
                                           (KO_dunn.adj_desc %>%
                                              dplyr::filter(q_SampleRT < 0.2) %>%
                                              dplyr::filter(passed_ss_SampleRT == 1))$`function`) ,
                              qvalueCutoff = 0.05)
dp_KO.RT <- dotplot(KO_enrich_Dunn.RT, showCategory = 15) + ggtitle("RT - GT (Dunn)")
write_xlsx(KO_enrich_Dunn.RT@result %>%
             dplyr::filter(qvalue < 0.05), "results/PICRUSt2/KO_enrich_Dunn.RT.sig.xlsx")
KO_enrich_Dunn.ST <- enrichKO(dplyr::union((KO_dunn_desc %>%
                                              dplyr::filter(q_SampleST < 0.2) %>%
                                              dplyr::filter(passed_ss_SampleST == 1))$`function` ,
                                           (KO_dunn.adj_desc %>%
                                              dplyr::filter(q_SampleST < 0.2) %>%
                                              dplyr::filter(passed_ss_SampleST == 1))$`function`) , 
                              qvalueCutoff = 0.05)
dp_KO.ST <- dotplot(KO_enrich_Dunn.ST, showCategory = 15) + ggtitle("ST - GT (Dunn)")
write_xlsx(KO_enrich_Dunn.ST@result %>%
             dplyr::filter(qvalue < 0.05), "results/PICRUSt2/KO_enrich_Dunn.ST.sig.xlsx")

ggarrange(dp_KO, dp_KO.RT, dp_KO.ST, ncol = 3)

