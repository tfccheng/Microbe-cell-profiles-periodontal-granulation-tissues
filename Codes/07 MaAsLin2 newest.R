library(Maaslin2)
library(dplyr)
library(mia)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(tidyr)
library(writexl)
library(ComplexHeatmap)

setwd("")

set.seed(123)
tse.micro.BE.new1 <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds")
pseq.micro.BE.new1 <- makePhyloseqFromTreeSummarizedExperiment(tse.micro.BE.new1)

meta_data <- colData(tse.micro.BE.new1) %>% data.frame() %>% 
  mutate(Obese = case_when(Obese == 1 ~ "Normal",
                           Obese == 2 ~ "Obese"),
         Obese = factor(Obese))


tse.micro.BE.new1 <- tse.micro.BE.new1 %>%
  transformAssay(assay.type = "counts", method = "relabundance") 


asv.micro <- t(as.data.frame(assay(tse.micro.BE.new1, "relabundance")))
asv.micro.fam <- t(as.data.frame(assay(altExp(tse.micro.BE.new1, "Family"), "relabundance"))) %>%
  as.data.frame() %>% dplyr::select(-Mitochondria) # remove mitochondria
asv.micro.genus <- t(as.data.frame(assay(altExp(tse.micro.BE.new1, "Genus"), "relabundance"))) %>%
  as.data.frame() 
asv.micro.sp <- t(as.data.frame(assay(altExp(tse.micro.BE.new1, "Species"), "relabundance"))) %>%
  as.data.frame() 
tax.family <- data.frame(rowData(altExp(tse.micro.BE.new1, "Family"))) %>%
  filter(Family != "Mitochondria") # Remove mitochondria
tax.genus <- data.frame(rowData(altExp(tse.micro.BE.new1, "Genus")))
tax.sp <- data.frame(rowData(altExp(tse.micro.BE.new1, "Species")))

####level with periodontal effects####

#####family####
fit_data_ASV.fam <- Maaslin2(
  input_data = asv.micro.fam,
  input_metadata = meta_data,
  output = "results/MaAsLin2/output_Family",
  transform = "NONE", 
  analysis_method = "CPLM",
  fixed_effects = c("Sample"),
  reference = "Sample,GT",  
  normalization = "TSS",
  standardize = TRUE,
  min_prevalence = 0.05,
  cores= 2
)
maaslin2_family <-fit_data_ASV.fam


fit_data_ASV.fam.adj <- Maaslin2(
  input_data = asv.micro.fam,
  input_metadata = meta_data,
  output = "results/MaAsLin2/output_Family_adj",
  transform = "NONE", 
  analysis_method = "CPLM",
  fixed_effects = c("Sample", "Sex", "Age", "Obese", "Smoker", "DM"),
  reference = "Sample,GT",  
  normalization = "TSS",
  standardize = TRUE,
  min_prevalence = 0.05,
  cores= 2
)
maaslin2_family.adj <-fit_data_ASV.fam.adj

fit_data_ASV.genus <- Maaslin2(
  input_data = asv.micro.genus,
  input_metadata = meta_data,
  output = "results/MaAsLin2/output_Genus",
  transform = "NONE", 
  analysis_method = "CPLM",
  fixed_effects = c("Sample"),
  reference = "Sample,GT",  
  normalization = "TSS",
  standardize = TRUE,
  min_prevalence = 0.05,
  cores= 2
)
maaslin2_genus <-fit_data_ASV.genus

fit_data_ASV.genus.adj <- Maaslin2(
  input_data = asv.micro.genus,
  input_metadata = meta_data,
  output = "results/MaAsLin2/output_Genus_adj",
  transform = "NONE", 
  analysis_method = "CPLM",
  fixed_effects = c("Sample", "Sex", "Age", "Obese", "Smoker", "DM"),
  reference = "Sample,GT",  
  normalization = "TSS",
  standardize = TRUE,
  min_prevalence = 0.05,
  cores= 2
)
maaslin2_genus.adj <-fit_data_ASV.genus.adj



fit_data_ASV.sp <- Maaslin2(
  input_data = asv.micro.sp,
  input_metadata = meta_data,
  output = "results/MaAsLin2/output_Species",
  transform = "NONE", 
  analysis_method = "CPLM",
  fixed_effects = c("Sample"),
  reference = "Sample,GT",  
  normalization = "TSS",
  standardize = TRUE,
  min_prevalence = 0.05,
  cores= 2
)
maaslin2_sp <-fit_data_ASV.sp

fit_data_ASV.sp.adj <- Maaslin2(
  input_data = asv.micro.sp,
  input_metadata = meta_data,
  output = "results/MaAsLin2/output_Species_adj",
  transform = "NONE", 
  analysis_method = "CPLM",
  fixed_effects = c("Sample", "Sex", "Age", "Obese", "Smoker", "DM"),
  reference = "Sample,GT",  
  normalization = "TSS",
  standardize = TRUE,
  min_prevalence = 0.05,
  cores= 2
)
maaslin2_sp.adj <-fit_data_ASV.sp.adj



####plot####
#####Family####
df.family <- maaslin2_family$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  # No asterisk for p-value >= 0.05
  ))
df.family$upper_ci <- df.family$coef + 1.96*df.family$stderr
df.family$lower_ci <- df.family$coef - 1.96*df.family$stderr
write_xlsx(df.family %>% dplyr::rename(Family = feature) %>% 
             dplyr::left_join(tax.family, by = "Family"), "Figures/MaAsLin2/df.family.tax.xlsx")
df.family.adj <- maaslin2_family.adj$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  filter(metadata == "Sample") %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  
  ))
df.family.adj$upper_ci <- df.family.adj$coef + 1.96*df.family.adj$stderr
df.family.adj$lower_ci <- df.family.adj$coef - 1.96*df.family.adj$stderr
df.family.adj.mod <- df.family.adj %>%
  mutate( feature = case_when(feature == 'Family.XI' ~ 'Family XI',
                              feature == 'SC.I.84' ~ 'SC-I-84',
                              feature == 'Rhizobiales.Incertae.Sedis' ~ 'Rhizobiales Incertae Sedis',
                              feature == 'Bacteroidales.Incertae.Sedis' ~ 'Bacteroidales Incertae Sedis',
                              .default = feature))
write_xlsx(df.family.adj.mod %>% dplyr::rename(Family = feature) %>% 
             dplyr::left_join(tax.family, by = "Family"), "Figures/MaAsLin2/df.family.adj.tax.xlsx")
df.family.adj.all <- maaslin2_family.adj$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  
  ))
df.family.adj.all$upper_ci <- df.family.adj.all$coef + 1.96*df.family.adj.all$stderr
df.family.adj.all$lower_ci <- df.family.adj.all$coef - 1.96*df.family.adj.all$stderr

# manually correct feature names
df.family.adj.all.mod <- df.family.adj.all %>%
  mutate( feature = case_when(feature == 'Rhizobiales.Incertae.Sedis' ~ 'Rhizobiales Incertae Sedis',
                              feature == 'Bacteroidetes.vadinHA17' ~ 'Bacteroidetes vadinHA17',
                              feature == 'Acidobacteriaceae..Subgroup.1.' ~ 'Acidobacteriaceae (Subgroup 1)',
                              feature == 'Amb.16S.1323' ~ 'Amb-16S-1323',
                              feature == 'Bacteroidales.Incertae.Sedis' ~ 'Bacteroidales Incertae Sedis',
                              feature == 'Family.XI' ~ 'Family XI',
                              feature == 'JG30.KF.AS9' ~ 'JG30-KF-AS9',
                              feature == 'SC.I.84' ~ 'SC-I-84',
                              feature == 'Unknown.Family' ~ 'Unknown Family',
                              feature == 'X67.14' ~ '67-14',
                              feature == 'TRA3.20' ~ 'TRA3-20',
                              feature == 'UCG.010' ~ 'UCG-010',
                              feature == 'X.Eubacterium..coprostanoligenes.group' ~ '[Eubacterium] coprostanoligenes group',
                              .default = feature))

write_xlsx(df.family.adj.all.mod %>% dplyr::rename(Family = feature) %>% 
             dplyr::left_join(tax.family, by = "Family"), "Figures/MaAsLin2/df.family.adj.all.tax.xlsx")

summary(maaslin2_family.adj$results$N.not.zero) #check the mean of N.not.zero

p.coef.family <- ggdotplot(data = df.family %>% filter(`N.not.zero` > 7), x = "feature", y = "coef", # sample number of RT is 8
                           color = "value", fill = "value", rotate = T, legend = "none") +
  labs(x= "", y = "Effect size") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, color = value),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.75) +
  scale_color_manual(name = "value", values = c("#00a087", "#e64b35", "#631879"))+
  scale_fill_manual(name = "value", values = c("#00a087", "#e64b35", "#631879")) +
  geom_text(aes(y= 5, label = significance, color = value), size = 8, nudge_x = -0.1) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))
#plot for adj further control by q<0.1
p.coef.family.adj <- ggdotplot(data = df.family.adj.mod %>% filter(`N.not.zero` > 7) %>% filter(qval < 0.1),
                               x = "feature", y = "coef", 
                           color = "value", fill = "value", rotate = T, size = 0.75, legend = "none") +
  labs(x= "", y = "Covariate-adjusted Effect size") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, color = value),
                position = position_dodge(width = 0.8), width = 0.5, linewidth = 0.75) +
  scale_color_manual(name = "value", values = c("#00a087", "#e64b35", "#631879"))+
  scale_fill_manual(name = "value", values = c("#00a087", "#e64b35", "#631879")) +
  geom_text(aes(y=7.5, label = significance, color = value), size = 8, nudge_x = -0.4) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))


p.rel.family <- ggboxplot(melt(cbind(log2( data.frame(asv.micro.fam)[,unique( df.family %>% 
                                                                               filter(`N.not.zero` > 7) %>%
                                                                               pull(feature) )]+ 1E-6), meta_data[,1:2])),x = "variable", y = "value",
                          color = "Sample", fill = "Sample", alpha = 0.2, rotate = T, legend = "none",
                          add = "jitter") +
  labs(x ="", y="log2(Rel. abundance)") + 
  scale_color_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  scale_fill_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  stat_compare_means(aes(group = Sample), method = "kruskal.test", label = "p.signif", vjust = 0.5, size = 8, hide.ns = T) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))
 

ggarrange(p.rel.family, 
          p.coef.family + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank()),
          p.coef.family.adj, ncol = 3, widths = c(2,1,2))
  
#####Genus####
df.genus <- maaslin2_genus$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  # No asterisk for p-value >= 0.05
  )) 
  
df.genus$upper_ci <- df.genus$coef + 1.96*df.genus$stderr
df.genus$lower_ci <- df.genus$coef - 1.96*df.genus$stderr
df.genus.mod <- df.genus %>%
  mutate( feature = case_when(#feature == 'X.Eubacterium..saphenum.group' ~ '[Eubacterium] saphenum group',
                              #feature == 'X.Clostridium..innocuum.group' ~ '[Clostridium] innocuum group',
                              #feature == 'Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium' ~ 'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium',
                              feature == 'Burkholderia.Caballeronia.Paraburkholderia' ~ 'Burkholderia-Caballeronia-Paraburkholderia',
                              #feature == 'Methylobacterium.Methylorubrum' ~ 'Methylobacterium-Methylorubrum',
                              feature == 'X.Eubacterium..hallii.group' ~ '[Eubacterium] hallii group',
                              #feature == 'X.Eubacterium..xylanophilum.group' ~ '[Eubacterium] xylanophilum group',
                              feature == 'BCf9.17.termite.group' ~ 'BCf9-17 termite group',
                              .default = feature))
write_xlsx(df.genus.mod %>% dplyr::rename(Genus = feature) %>% 
             dplyr::left_join(tax.genus, by = "Genus"), "Figures/MaAsLin2/df.genus.tax.xlsx")
df.genus.adj <- maaslin2_genus.adj$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  filter(metadata == "Sample") %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  
  ))
df.genus.adj$upper_ci <- df.genus.adj$coef + 1.96*df.genus.adj$stderr
df.genus.adj$lower_ci <- df.genus.adj$coef - 1.96*df.genus.adj$stderr
df.genus.adj.mod <- df.genus.adj %>%
  mutate( feature = case_when(feature == 'X.Eubacterium..saphenum.group' ~ '[Eubacterium] saphenum group',
                              feature == 'X.Clostridium..innocuum.group' ~ '[Clostridium] innocuum group',
                              feature == 'Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium' ~ 'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium',
                              feature == 'Family.XIII.UCG.001' ~ 'Family XIII UCG-001',
                              feature == 'Burkholderia.Caballeronia.Paraburkholderia' ~ 'Burkholderia-Caballeronia-Paraburkholderia',
                              feature == 'X.Eubacterium..yurii.group' ~ '[Eubacterium] yurii group',
                              # feature == 'Methylobacterium.Methylorubrum' ~ 'Methylobacterium-Methylorubrum',
                              feature == 'X.Eubacterium..hallii.group' ~ '[Eubacterium] hallii group',
                              # feature == 'X.Eubacterium..xylanophilum.group' ~ '[Eubacterium] xylanophilum group',
                              # feature == 'BCf9.17.termite.group' ~ 'BCf9-17 termite group',
                              .default = feature))
write_xlsx(df.genus.adj.mod %>% dplyr::rename(Genus = feature) %>% 
             dplyr::left_join(tax.genus, by = "Genus"), "Figures/MaAsLin2/df.genus.adj.tax.xlsx")

df.genus.adj.all <- maaslin2_genus.adj$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  
  ))
df.genus.adj.all$upper_ci <- df.genus.adj.all$coef + 1.96*df.genus.adj.all$stderr
df.genus.adj.all$lower_ci <- df.genus.adj.all$coef - 1.96*df.genus.adj.all$stderr
df.genus.adj.all.mod <- df.genus.adj.all %>%
  mutate( feature = case_when(feature == 'X.Eubacterium..saphenum.group' ~ '[Eubacterium] saphenum group',
                              feature == 'X.Clostridium..innocuum.group' ~ '[Clostridium] innocuum group',
                              feature == 'Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium' ~ 'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium',
                              feature == 'Burkholderia.Caballeronia.Paraburkholderia' ~ 'Burkholderia-Caballeronia-Paraburkholderia',
                              #feature == 'Methylobacterium.Methylorubrum' ~ 'Methylobacterium-Methylorubrum',
                              #feature == 'X.Eubacterium..hallii.group' ~ '[Eubacterium] hallii group',
                              #feature == 'X.Eubacterium..xylanophilum.group' ~ '[Eubacterium] xylanophilum group',
                              #feature == 'BCf9.17.termite.group' ~ 'BCf9-17 termite group',
                              feature == 'Candidatus.Solibacter' ~ 'Candidatus Solibacter',
                              feature == 'Candidatus.Udaeobacter' ~ 'Candidatus Udaeobacter',
                              feature == 'Clostridium.sensu.stricto.1' ~ 'Clostridium sensu stricto 1',
                              feature == 'Defluviitaleaceae.UCG.011' ~ 'Defluviitaleaceae UCG-011',
                              feature == 'Erysipelotrichaceae.UCG.003' ~ 'Erysipelotrichaceae UCG-003',
                              feature == 'Family.XIII.AD3011.group' ~ 'Family XIII AD3011 group',
                              feature == 'Family.XIII.UCG.001' ~ 'Family XIII UCG-001',
                              feature == 'Incertae.Sedis' ~ 'Incertae Sedis',
                              feature == 'IS.44' ~ 'IS-44',
                              feature == 'NK4A214.group' ~ 'NK4A214 group',
                              feature == 'Prevotellaceae.Ga6A1.group' ~ 'Prevotellaceae Ga6A1 group',
                              feature == 'Rikenellaceae.RC9.gut.group' ~ 'Rikenellaceae RC9 gut group',
                              feature == 'X.Eubacterium..yurii.group' ~ '[Eubacterium] yurii group',
                              #feature == 'X.Ruminococcus..torques.group' ~ '[Ruminococcus] torques group',
                              .default = feature))
write_xlsx(df.genus.adj.all.mod %>% dplyr::rename(Genus = feature) %>% 
             dplyr::left_join(tax.genus, by = "Genus"), "Figures/MaAsLin2/df.genus.adj.all.tax.xlsx")

summary(maaslin2_genus.adj$results$N.not.zero) #check the mean of N.not.zero

p.coef.genus <- ggdotplot(data = df.genus.mod %>% filter(`N.not.zero` > 7), x = "feature", y = "coef", 
                           color = "value", fill = "value", rotate = T, legend = "none") +
  labs(x= "", y = "Effect size") +
  scale_x_discrete(labels = function(x) stringr::str_trunc(x, 25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, color = value),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.75) +
  scale_color_manual(name = "value", values = c("#00a087", "#e64b35", "#631879"))+
  scale_fill_manual(name = "value", values = c("#00a087", "#e64b35", "#631879")) +
  geom_text(aes(y=8, label = significance, color = value), size = 8, nudge_x = -0.15)+
  theme(axis.text.y = element_text(face="italic", size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))

p.coef.genus.adj <- ggdotplot(data = df.genus.adj.mod %>% filter(`N.not.zero` > 7) %>% filter(qval < 0.1), 
                              x = "feature", y = "coef", 
                               color = "value", fill = "value", rotate = T, size = 0.75, legend = "none") +
  labs(x= "", y = "Covariate-adjusted Effect size") +
  scale_x_discrete(labels = function(x) stringr::str_trunc(x, 25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, color = value),
                position = position_dodge(width = 0.8), width = 0.5, linewidth = 0.75) +
  scale_color_manual(name = "value", values = c("#00a087", "#e64b35", "#631879"))+
  scale_fill_manual(name = "value", values = c("#00a087", "#e64b35", "#631879")) +
  geom_text(aes(y=7.5, label = significance, color = value), size = 8, nudge_x = -0.4) +
  theme(axis.text.y = element_text(face="italic", size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))


p.rel.genus <- ggboxplot(melt(cbind(log2( data.frame(asv.micro.genus)[,unique(df.genus%>% 
                                                                                filter(`N.not.zero` > 7) %>%
                                                                                pull(feature) )]+ 1E-6), meta_data[,1:2])),x = "variable", y = "value",
                          color = "Sample", fill = "Sample", alpha = 0.2, rotate = T, legend = "none",
                          add = "jitter") +
  labs(x ="", y="log2(Rel. abundance)") + 
  scale_x_discrete(labels = function(x) stringr::str_trunc(x, 25)) +
  scale_color_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  scale_fill_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  stat_compare_means(aes(group = Sample), method = "kruskal.test", label = "p.signif", vjust = 0.5, size = 8, hide.ns = T)+
  theme(axis.text.y = element_text(face="italic", size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))


ggarrange(p.rel.genus, 
          p.coef.genus + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank()),
          p.coef.genus.adj, ncol = 3, widths = c(2,1,2))

#####Species####
df.sp <- maaslin2_sp$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  # No asterisk for p-value >= 0.05
  ))
df.sp$upper_ci <- df.sp$coef + 1.96*df.sp$stderr
df.sp$lower_ci <- df.sp$coef - 1.96*df.sp$stderr
df.sp.mod <- df.sp %>%
  mutate( feature = case_when(feature == 'Methylobacterium.Methylorubrum_extorquens' ~ 'Methylobacterium-Methylorubrum_extorquens',
                              feature == 'Methylobacterium.Methylorubrum_brachiatum' ~ 'Methylobacterium-Methylorubrum_brachiatum',
                              feature == 'X.Eubacterium..hallii_hallii' ~ '[Eubacterium] hallii_hallii',
                              .default = feature))
write_xlsx(df.sp.modify %>% dplyr::rename(Species = feature) %>% 
             dplyr::left_join(tax.sp, by = "Species"), "Figures/MaAsLin2/df.sp.tax.xlsx")
df.sp.adj <- maaslin2_sp.adj$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  filter(metadata == "Sample") %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  
  ))
df.sp.adj$upper_ci <- df.sp.adj$coef + 1.96*df.sp.adj$stderr
df.sp.adj$lower_ci <- df.sp.adj$coef - 1.96*df.sp.adj$stderr
df.sp.adj.mod <- df.sp.adj %>%
  mutate( feature = case_when(feature == 'X.Eubacterium..saphenum_saphenum' ~ '[Eubacterium] saphenum_saphenum',
                              .default = feature))
write_xlsx(df.sp.adj %>% dplyr::rename(Species = feature) %>% 
             dplyr::left_join(tax.sp, by = "Species"), "Figures/MaAsLin2/df.sp.adj.tax.xlsx")

df.sp.adj.all <- maaslin2_sp.adj$results %>%
  filter(qval < 0.25) %>%
  arrange(coef) %>%
  mutate(significance = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01 ~ "**",
    qval < 0.05 ~ "*",
    TRUE ~ ""  
  ))
df.sp.adj.all$upper_ci <- df.sp.adj.all$coef + 1.96*df.sp.adj.all$stderr
df.sp.adj.all$lower_ci <- df.sp.adj.all$coef - 1.96*df.sp.adj.all$stderr
df.sp.adj.all.mod <- df.sp.adj.all %>%
  mutate( feature = case_when(feature == 'X.Eubacterium..saphenum_saphenum' ~ '[Eubacterium] saphenum_saphenum',
                              feature == 'Clostridium.sensu.stricto.1_butyricum' ~ 'Clostridium sensu stricto 1_butyricum',
                              feature == 'Escherichia.Shigella_coli' ~ 'Escherichia-Shigella_coli',
                              .default = feature))
write_xlsx(df.sp.adj.all %>% dplyr::rename(Species = feature) %>% 
             dplyr::left_join(tax.sp, by = "Species"), "Figures/MaAsLin2/df.sp.adj.all.tax.xlsx")

summary(maaslin2_sp.adj$results$N.not.zero) #check the mean of N.not.zero

p.coef.sp <- ggdotplot(data = df.sp.mod %>% filter(`N.not.zero` > 7) %>% filter(qval < 0.1) %>% 
                         filter(feature != "Bacillus_virus"), 
                       x = "feature", y = "coef", 
                           color = "value", fill = "value", rotate = T, legend = "none") +
  labs(x= "", y = "Effect size") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, color = value),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.75) +
  scale_color_manual(name = "value", values = c("#00a087", "#e64b35", "#631879"))+
  scale_fill_manual(name = "value", values = c("#00a087", "#e64b35", "#631879")) +
  geom_text(aes(y=11, label = significance, color = value), size = 8, nudge_x = -0.1)+
  theme(axis.text.y = element_text(face="italic", size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))

p.coef.sp.adj <- ggdotplot(data = df.sp.adj %>% filter(`N.not.zero` > 7) %>% filter(qval < 0.1) %>% 
                             filter(feature != "Bacillus_virus"), 
                           x = "feature", y = "coef", 
                               color = "value", fill = "value", rotate = T, size = 0.75, legend = "none") +
  labs(x= "", y = "Covariate-adjusted Effect size") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci, color = value),
                position = position_dodge(width = 0.8), width = 0.5, linewidth = 0.75) +
  scale_color_manual(name = "value", values = c("#00a087", "#e64b35", "#631879"))+
  scale_fill_manual(name = "value", values = c("#00a087", "#e64b35", "#631879")) +
  geom_text(aes(y=10, label = significance, color = value), size = 8)+
  theme(axis.text.y = element_text(face="italic", size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))


p.rel.sp <- ggboxplot(melt(cbind(log2( data.frame(asv.micro.sp)[,unique(df.sp%>% 
                                                                          filter(`N.not.zero` > 7) %>% filter(qval < 0.1) %>% 
                                                                          filter(feature != "Bacillus_virus") %>%
                                                                          pull(feature) )]+ 1E-6), meta_data[,1:2])),x = "variable", y = "value",
                          color = "Sample", fill = "Sample", alpha = 0.2, rotate = T, legend = "none",
                          add = "jitter") +
  labs(x ="", y="log2(Rel. abundance)") + 
  scale_color_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  scale_fill_manual(name = "Sample", values = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  stat_compare_means(aes(group = Sample), method = "kruskal.test", label = "p.signif", vjust = 0.5, size = 8, hide.ns = T)+
  theme(axis.text.y = element_text(face="italic", size = 14),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16))


ggarrange(p.rel.sp, 
          p.coef.sp + theme(axis.title.y = element_blank(),
                                axis.text.y = element_blank()),
          p.coef.sp.adj, ncol = 3, widths = c(2,1,2))


####draw heatmap####
library(openxlsx)
family.sig <- read.xlsx("Figures/MaAsLin2/result summary.xlsx", sheet = "Family") %>%
  filter(`N.not.zero` > 7)
family.sig.adj <- read.xlsx("Figures/MaAsLin2/result summary.xlsx", sheet = "Family.adj") %>%
  filter(metadata == "Sample") %>%
  filter(`N.not.zero` > 7) %>% filter(qval < 0.1)
family.union <- union(unique(family.sig$Family), unique(family.sig.adj$Family))
# genus.sig <- read.xlsx("Figures/MaAsLin2/result summary.xlsx", sheet = "Genus") %>%
#   filter(`N.not.zero` > 7)
# genus.sig.adj <- read.xlsx("Figures/MaAsLin2/result summary.xlsx", sheet = "Genus.adj") %>%
#   filter(metadata == "Sample") %>%
#   filter(`N.not.zero` > 7)
# genus.union <- union(unique(genus.sig$Genus), unique(genus.sig.adj$Genus))
# species.sig <- read.xlsx("Figures/MaAsLin2/result summary.xlsx", sheet = "Species") %>%
#   filter(`N.not.zero` > 7)
# species.sig.adj <- read.xlsx("Figures/MaAsLin2/result summary.xlsx", sheet = "Species.adj") %>%
#   filter(metadata == "Sample") %>%
#   filter(`N.not.zero` > 7)
# species.union <- union(unique(species.sig$Species), unique(species.sig.adj$Species))

fam.Z.sig.select <- as.data.frame(assay(altExp(tse.micro.BE.new1, "Family")%>%
                                         transformAssay(assay.type = "counts", method = "clr", pseudocount=1, 
                                                        MARGIN = "samples") %>%
                                         transformAssay(assay.type = "clr", method = "z", 
                                                        name = "clr_z", 
                                                        MARGIN = "features"), "clr_z"))[family.union,] %>% t() %>%
  cbind(meta_data %>% dplyr::select(Sample)) %>% #bind meta data of Sample
  arrange(Sample) 
mean.fam.rel.sig.select <- as.data.frame(assay(altExp(tse.micro.BE.new1, "Family")%>%
                                                 transformAssay(assay.type = "counts", MARGIN = "features", 
                                                                method = "relabundance"), #transform to relabundance by each taxa
                                               "relabundance"))[family.union,] %>% t() %>%
  cbind(meta_data %>% dplyr::select(Sample)) %>% #bind meta data of Sample
  arrange(Sample) %>%
  dplyr::group_by(Sample) %>% #calculate mean by sample groups for every taxa
  summarise_at(vars(all_of(family.union)), mean) %>% 
  tibble::column_to_rownames("Sample") %>% t()


df.fam.Z.sig.select <- fam.Z.sig.select  %>% # order row by Sample
  dplyr::select(-Sample) %>% t()

#####heatmap####
sample_col <- meta_data %>% arrange(Sample) %>% # order by Sample levels
  dplyr::select(Sample)

class.family <- tax.family[family.union,] %>% dplyr::pull(Class)
class.genus <- tax.genus[genus.union,] %>% dplyr::pull(Class)
class.sp <- tax.sp[species.union,] %>% dplyr::pull(Class)
cl_all <- union(class.family,
                   union(class.genus,
                         class.sp))


require(RColorBrewer)

col_fun = colorRampPalette(colors = rev(brewer.pal(11,"RdBu")))


col1 <- setNames(c("#0f77c1", "#00a087", "#e64b35", "#631879"), c("GT", "PT", "RT", "ST"))
col2 <- setNames(c(paletteer::paletteer_d("ggsci::category20_d3", 
                                          length(cl_all))), c(cl_all))


ha1 = HeatmapAnnotation(Sample = sample_col$Sample, 
                        show_annotation_name = F, show_legend = T,
                        simple_anno_size = unit(0.3,"cm"),
                        col = list(`Sample` = col1))
ha2 = rowAnnotation(Class = class.family, show_annotation_name = F,
                    simple_anno_size = unit(0.3,"cm"), 
                    col = list(Class = col2[class.family]))
ha3 <- rowAnnotation(block = anno_block(gp = gpar(fill = 2:5, col = NA),
                                        labels = c("1", "2", "3", "4")),
                     show_annotation_name = F, width = unit(0.3, "cm"))
ha4 <- rowAnnotation(`Mean\nrelabundance` = anno_barplot(mean.fam.rel.sig.select , 
                                                   gp = gpar(fill = c("#0f77c1", "#00a087", "#e64b35", "#631879")), 
                                                   width = unit(1, "cm")))
p.family <- ha4 + Heatmap(df.fam.Z.sig.select, name = "z score (clr)",
                     top_annotation = ha1, right_annotation = ha2, left_annotation = ha3,
                     cluster_rows = T, cluster_columns = T, 
                     column_split = sample_col$Sample, cluster_column_slices = T,
                     row_names_side = "right", column_names_side = "bottom",
                     show_column_names = F, row_km = 4, 
                     row_title = "", row_title_side = "right",
                     column_title=NULL,
                     col = col_fun(9),
                     width = ncol(df.fam.Z.sig.select)*unit(2, "mm"), 
                     height = nrow(df.fam.Z.sig.select)*unit(4, "mm")) 
  
p.family



ha2 = rowAnnotation(Class = class.genus, show_annotation_name = F,
                    simple_anno_size = unit(0.3,"cm"), 
                    col = list(Class = col2[class.genus]))

ha4 <- rowAnnotation(`Mean\nrelabundance` = anno_barplot(mean.genus.rel.sig.select , 
                                                         gp = gpar(fill = c("#0f77c1", "#00a087", "#e64b35", "#631879")), 
                                                         width = unit(1, "cm")))
p.genus <- ha4 + Heatmap(df.genus.Z.sig.select, name = "z score (clr)",
                          top_annotation = ha1, right_annotation = ha2, left_annotation = ha3,
                          cluster_rows = T, cluster_columns = T, 
                          column_split = sample_col$Sample, cluster_column_slices = T,
                          row_names_side = "right", column_names_side = "bottom",
                          show_column_names = F, row_km = 4, 
                          row_names_gp = gpar(fontface = "italic"),
                          row_title = "", row_title_side = "right",
                          column_title=NULL,
                          col = col_fun(9),
                          width = ncol(df.genus.Z.sig.select)*unit(2, "mm"), 
                          height = nrow(df.genus.Z.sig.select)*unit(4, "mm")) 

p.genus




ha2 = rowAnnotation(Class = class.sp, show_annotation_name = F,
                    simple_anno_size = unit(0.3,"cm"), 
                    col = list(Class = col2[class.sp]))

ha4 <- rowAnnotation(`Mean\nrelabundance` = anno_barplot(mean.species.rel.sig.select , 
                                                         gp = gpar(fill = c("#0f77c1", "#00a087", "#e64b35", "#631879")), 
                                                         width = unit(1, "cm")))
p.species <- ha4 + Heatmap(df.species.Z.sig.select, name = "z score (clr)",
                         top_annotation = ha1, right_annotation = ha2, left_annotation = ha3,
                         cluster_rows = T, cluster_columns = T, 
                         column_split = sample_col$Sample, cluster_column_slices = T,
                         row_names_side = "right", column_names_side = "bottom",
                         show_column_names = F, row_km = 4, 
                         row_names_gp = gpar(fontface = "italic"),
                         row_title = "", row_title_side = "right",
                         column_title=NULL,
                         col = col_fun(9),
                         width = ncol(df.species.Z.sig.select)*unit(2, "mm"), 
                         height = nrow(df.species.Z.sig.select)*unit(4, "mm")) 

p.species



