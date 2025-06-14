library(phyloseq)
library(mia)
library(dplyr)
library("readxl")
library(mixOmics) # import the mixOmics library
library(gridExtra)
library(readxl)
library(writexl)
library(ggpubr)
#for Windows
setwd("")

set.seed(123)

tse.micro.BE.new1 <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds")




X.samples.BE.new1 <- tse.micro.BE.new1 %>%
  subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05) %>% #filter with 0.05% detection & 5% prevalence
  #altExp("Genus") %>%
  transformAssay(method = "relabundance") %>%
  transformAssay(assay.type = "relabundance", method = "rclr") %>% #transform counts to robust clr
  assay("rclr") %>%
  t()

sample_meta <- data.frame(colData(tse.micro.BE.new1))
Y.samples <- sample_meta$Sample %>%
  factor(levels = c("GT", "PT", "RT", "ST"))
levels(Y.samples) <- c("Granulation", "Periodontium", "Root", "Socket")




#############################################################################
####sPLS-DA of ASV level####



# pca 

pca.BE.new1 = pca(X.samples.BE.new1, ncomp = 10) 

#plot(pca)  # barplot of the eigenvalues 
plotIndiv(pca.BE.new1, group = Y.samples, ind.names = FALSE, # plot the samples 
          legend = TRUE,  ellipse = T, title = 'PCA, comp 1 - 2') 




# initial sPLS-DA

splsda.BE.new1 <- splsda(X.samples.BE.new1, near.zero.var = T, Y.samples, ncomp = 10)  
 
# plot the samples projected onto the first two components of the PLS-DA subspace
# par(mfrow=c(2, 2))

plotIndiv(splsda.BE.new1 , comp = 1:2, 
          group = Y.samples, ind.names = FALSE, 
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# # use the max.dist measure to form decision boundaries between classes based on PLS-DA data
# background = background.predict(splsda, comp.predicted=2, dist = "max.dist")
# 
# # plot the samples projected onto the first two components of the PLS-DA subspace
# plotIndiv(splsda, comp = 1:2,
#           group = Y, ind.names = FALSE, # colour points by class
#           background = background, # include prediction background for each class
#           legend = TRUE, title = " (b) PLSDA with prediction background")

#tuning sPLS-DA
# undergo performance evaluation in order to tune the number of components to use

perf.splsda.BE.new1 <- perf(splsda.BE.new1, validation = "Mfold", 
                       folds = 8, nrepeat = 10, # use repeated cross-validation
                       progressBar = TRUE, cpus = 14, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components

plot(perf.splsda.BE.new1, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.BE.new1$choice.ncomp

# grid of possible keepX values that will be tested for each component
list.keepX <- c(seq(20, 300, 20), seq(350, 1500, 50))

# undergo the tuning process to determine the optimal number of variables

tune.splsda.BE.new1 <- tune.splsda(X.samples.BE.new1, Y.samples, ncomp = 2, 
                              validation = 'Mfold',
                              folds = 8, nrepeat = 50, # use repeated cross-validation
                              dist = 'max.dist', # use max.dist measure
                              measure = "BER", # use balanced error rate of dist measure
                              test.keepX = list.keepX,
                              cpus = 14, progressBar = T) # allow for parallel computing


plot(tune.splsda.BE.new1, col = color.jet(2)) # plot output of variable number tuning
tune.splsda.BE.new1$choice.ncomp$ncomp
optimal.ncomp.BE.new1 <- tune.splsda.BE.new1$choice.ncomp$ncomp
optimal.keepX.BE.new1 <- tune.splsda.BE.new1$choice.keepX 



# final sPLS-DA

final.splsda.BE.new1 <- splsda(X.samples.BE.new1, Y.samples, 
                          ncomp = 2, 
                          keepX = optimal.keepX.BE.new1)
plot.final.splsda.BE.new1 <- plotIndiv(final.splsda.BE.new1, comp = c(1,2), 
                                  group = Y.samples, ind.names = FALSE, 
                                  ellipse = TRUE, legend = TRUE, 
                                  abline = T, col = c ("#0f77c1", "#00a087", "#e64b35", "#631879"),
                                  title = 'Tissues')








#############################################################################
####sPLS-DA of genus level####

X.samples.BE.new1.g <- tse.micro.BE.new1 %>%
  subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05) %>% #filter with 0.05% detection & 5% prevalence
  altExp("Genus") %>%
  transformAssay(method = "relabundance") %>%
  transformAssay(assay.type = "relabundance", method = "rclr") %>% #transform counts to robust clr
  assay("rclr") %>%
  t()


# pca 

pca.BE.new1.g = pca(X.samples.BE.new1.g, ncomp = 10) 

#plot(pca)  # barplot of the eigenvalues 

plotInd.pca.BE.new1.g <- plotIndiv(pca.BE.new1.g, group = Y.samples, ind.names = FALSE, # plot the samples 
          legend = TRUE,  ellipse = T, pch = 19,
          col.per.group = c ("#0f77c1", "#00a087", "#e64b35", "#631879"), 
          title = 'PCA, comp 1 - 2') 


biplot.pca.BE.new1.g <- biplot(pca.BE.new1.g, cex = 1, ind.names = F, 
                               group = Y.samples, var.names.size = 3,# colour by sample class
                              legend.title = '', var.names.col = "black", 
                              col.per.group = c ("#0f77c1", "#00a087", "#e64b35", "#631879"), 
                              title = 'PCA comp 1 - 2')
pl.biplot.pca.BE.new1.g <- biplot.pca.BE.new1.g +
  stat_ellipse(data = data.frame(pca.BE.new1.g$variates$X), 
               aes(PC1, PC2, color = Y.samples), type = "t", level = 0.95,
               lwd = 1) +
  ggrepel::geom_text_repel(max.overlaps = 50)+ 
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.position = "bottom",
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        aspect.ratio = 1)
pl.plotInd.pca.BE.new1.g <- plotInd.pca.BE.new1.g$graph +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1)


# initial sPLS-DA

splsda.BE.new1.g <- splsda(X.samples.BE.new1.g, near.zero.var = T, Y.samples, ncomp = 10)  

# plot the samples projected onto the first two components of the PLS-DA subspace
# par(mfrow=c(2, 2))

plotIndiv(splsda.BE.new1.g , comp = 1:2, 
          group = Y.samples, ind.names = FALSE, 
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# # use the max.dist measure to form decision boundaries between classes based on PLS-DA data
# background = background.predict(splsda, comp.predicted=2, dist = "max.dist")
# 
# # plot the samples projected onto the first two components of the PLS-DA subspace
# plotIndiv(splsda, comp = 1:2,
#           group = Y, ind.names = FALSE, # colour points by class
#           background = background, # include prediction background for each class
#           legend = TRUE, title = " (b) PLSDA with prediction background")

#tuning sPLS-DA
# undergo performance evaluation in order to tune the number of components to use

perf.splsda.BE.new1.g <- perf(splsda.BE.new1.g, validation = "Mfold", 
                            folds = 8, nrepeat = 10, # use repeated cross-validation
                            progressBar = TRUE, cpus = 14, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components

plot(perf.splsda.BE.new1.g, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.BE.new1.g$choice.ncomp

# grid of possible keepX values that will be tested for each component
list.keepX.g <- c(seq(20, 400, 20))

# undergo the tuning process to determine the optimal number of variables

tune.splsda.BE.new1.g <- tune.splsda(X.samples.BE.new1.g, Y.samples, ncomp = 3, 
                                   validation = 'Mfold',
                                   folds = 8, nrepeat = 50, # use repeated cross-validation
                                   dist = 'max.dist', # use max.dist measure
                                   measure = "BER", # use balanced error rate of dist measure
                                   test.keepX = list.keepX.g,
                                   cpus = 14, progressBar = T) # allow for parallel computing

plot(tune.splsda.BE.new1.g, col = color.jet(3)) # plot output of variable number tuning
tune.splsda.BE.new1.g$choice.ncomp$ncomp
optimal.ncomp.BE.new1.g <- tune.splsda.BE.new1.g$choice.ncomp$ncomp
optimal.keepX.BE.new1.g <- tune.splsda.BE.new1.g$choice.keepX 



# final sPLS-DA

final.splsda.BE.new1.g <- splsda(X.samples.BE.new1.g, Y.samples, 
                                 ncomp = 2, 
                                 keepX = optimal.keepX.BE.new1.g)
plot.final.splsda.BE.new1.g <- plotIndiv(final.splsda.BE.new1.g, comp = c(1,2), 
                                       group = Y.samples, ind.names = FALSE, 
                                       ellipse = TRUE, legend = TRUE, 
                                       abline = T, col = c("#0f77c1", "#00a087", "#e64b35", "#631879"),
                                       title = '')
pl.plot.final.splsda.BE.new1.g <- plot.final.splsda.BE.new1.g$graph +
                                        theme(strip.background = element_blank(),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              axis.text=element_text(size=18),
                                              axis.title=element_text(size=18),
                                              legend.position = "bottom",
                                              legend.title = element_blank(),
                                              legend.text = element_text(size=16),
                                              aspect.ratio = 1)




legend=list(legend = levels(Y.samples), # set of classes
            col = c("#0f77c1", "#00a087", "#e64b35", "#631879"), # set of colours
            title = "Tissues", # legend title
            cex = 0.7) # legend size
color.sample <- Y.samples
levels(color.sample) <- c("#0f77c1", "#00a087", "#e64b35", "#631879")

cim_result <- cim(final.splsda.BE.new1.g, comp = c(1,2), 
    title = "Tissues",
    margins = c(2, 20), row.cex = 1,
    row.sideColors = color.sample, 
    row.names = T, col.names = T, transpose = T,
    legend = legend)
plotLoadings.final.splsda.BE.new1.g <- plotLoadings(final.splsda.BE.new1.g, comp = 1, size.name = 1,
                                          contrib = 'max', method = 'median', show.ties = F,
                                          title = "Tissue")

selectVar(final.splsda.BE.new1.g, comp =1)
df.pld.splsda.BE.new1.g <- plotLoadings.final.splsda.BE.new1.g$X %>%
  tibble::rownames_to_column("Genus")
write_xlsx(df.pld.splsda.BE.new1.g, 
           "Figures/sPLS-DA/df.pld.splsda.BE.new1.g.xlsx")
pl.pld.splsda.BE.new1.g <- ggplot(data=df.pld.splsda.BE.new1.g, aes(x=importance, y = reorder(Genus, -abs(importance)), fill=GroupContrib)) +
  geom_bar(stat="identity")+
  labs(title="", x = "Loading", y = "", fill = "Tissue") +
  scale_fill_manual(values=c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  scale_x_continuous(limits=c(-0.2, 0.2), expand = c(0,0)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        title=element_text(size=16),
        axis.text.x =element_text(size=14),
        axis.text.y =element_text(size=14, face = "italic"),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        #axis.ticks.y = element_blank(),
        #legend.position='none',
        legend.title = element_blank(), 
        legend.text = element_text(size=14),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pl.pld.splsda.BE.new1.g


ggarrange(pl.biplot.pca.BE.new1.g, 
          pl.plot.final.splsda.BE.new1.g)



####micro ratio####
tse.genus <-tse.micro.BE.new1 %>%
  subsetByPrevalentFeatures(detection = 1E-5, prevalence = 0.05) %>% #filter with 0.05% detection & 5% prevalence
  altExp("Genus")

df.PorRot <- t(as.data.frame(assay(tse.genus, "counts")))[, c("Porphyromonas","Rothia")] %>% 
  as.data.frame() %>%
  dplyr::mutate(PR_ratio = log2((Porphyromonas+1) / Rothia ))
write_xlsx(df.PorRot %>% tibble::rownames_to_column("ID"), "Figures/sPLS-DA/df.PorRot.xlsx")

p.g.PorRot <- ggplot(data.frame(df.PorRot,
                                sample_meta), 
                     aes(x=Sample, y=PR_ratio, color = Sample, fill = Sample)) + 
  geom_boxplot(notch=F, alpha = 0.2, outlier.shape = NA) +
  geom_jitter(shape=16, position=position_jitter(0.2))+ 
  scale_fill_manual(values=c("#0f77c1", "#00a087", "#e64b35", "#631879"), name = "") +
  scale_color_manual(values=c("#0f77c1", "#00a087", "#e64b35", "#631879"), name = "") +
  theme_classic() + labs(x = "", y = "Porphyromonas/Rothia\nlog2(count ratio)") + 
  theme(aspect.ratio= 1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(), legend.position = "none") +
  geom_pwc(aes(group = Sample), method = "wilcox.test", p.adjust.method = "fdr",label = "p.format", hide.ns = T) +
  stat_kruskal_test(aes(group = Sample), label.y = 20)


library(ROCR)
library(rcompanion)

temp <- data.frame(df.PorRot, sample_meta)
a<-temp %>% filter(Sample %in% c("GT", "PT")) 
a$Sample <- factor(a$Sample, levels = c("GT", "PT"))
pred1 <- prediction(a$PR_ratio,a$Sample)
roc.perf.PR.PT = performance(pred1, measure = "tpr", x.measure = "fpr")
auc.perf.PR.PT = performance(pred1, measure = "auc")
auc.perf.PR.PT@y.values[[1]]
a<-temp %>% filter(Sample %in% c("GT", "RT")) 
a$Sample <- factor(a$Sample, levels = c("GT", "RT"))
pred1 <- prediction(a$PR_ratio,a$Sample)
roc.perf.PR.RT = performance(pred1, measure = "tpr", x.measure = "fpr")
auc.perf.PR.RT = performance(pred1, measure = "auc")
auc.perf.PR.RT@y.values[[1]]
a<-temp %>% filter(Sample %in% c("GT", "ST")) 
a$Sample <- factor(a$Sample, levels = c("GT", "ST"))
pred1 <- prediction(a$PR_ratio,a$Sample)
roc.perf.PR.ST = performance(pred1, measure = "tpr", x.measure = "fpr")
auc.perf.PR.ST = performance(pred1, measure = "auc")
auc.perf.PR.ST@y.values[[1]]

df <- data.frame(curve=as.factor(rep(c(1,2,3), c(length(roc.perf.PR.PT@x.values[[1]]),
                                                 length(roc.perf.PR.RT@x.values[[1]]),
                                                 length(roc.perf.PR.ST@x.values[[1]])),                                     )), 
                 falsepositive=c(roc.perf.PR.PT@x.values[[1]],
                                 roc.perf.PR.RT@x.values[[1]],
                                 roc.perf.PR.ST@x.values[[1]]),
                 truepositive=c(roc.perf.PR.PT@y.values[[1]],
                                roc.perf.PR.RT@y.values[[1]],
                                roc.perf.PR.ST@y.values[[1]]))
plt <- ggplot(df, aes(x=falsepositive, y=truepositive, color=curve)) + 
  geom_line(size = .75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(name = "AUC", values=c("#00a087", "#e64b35", "#631879"),
                     labels = c("PT - GT: 0.379", 
                                "RT - GT: 0.754",
                                "ST - GT: 0.816")) +
  labs(x = "False-positive rate", y = "True-positive rate", title = "Porphyromonas/Rothia") +
  theme_bw()+
  theme(aspect.ratio= 1,
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.8, 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggarrange(p.g.PorRot, plt, labels = c("A", "B"))

####paired samples####
meta_pair <- sample_meta %>%
  tibble::rownames_to_column("ID") %>%
  group_by(Case) %>%
  filter(n() > 1) %>% #select paired samples
  mutate(Case_Sample = paste(Case, Sample, sep = "_"))

meta_pair_GT_PT <- meta_pair %>%
  filter(Sample == "GT" | Sample == "PT")#select GT and PT sample
write.csv(meta_pair_GT_PT, "Figures/sPLS-DA/meta_pair_GT_PT.csv")
meta_pair_RT_ST <- meta_pair %>%
  filter(Sample == "RT" | Sample == "ST")#select RT and ST sample
write.csv(meta_pair_RT_ST, "Figures/sPLS-DA/meta_pair_RT_ST.csv")

meta_pair_GT_PT_select <- read.csv("Figures/sPLS-DA/meta_pair_GT_PT_select.csv") 
meta_pair_RT_ST_select <- read.csv("Figures/sPLS-DA/meta_pair_RT_ST_select.csv")

p1<-df.PorRot %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::right_join(meta_pair_GT_PT_select, by = "ID") %>%
  arrange(Sample, Case) %>%
  dplyr::select(PR_ratio, Sample) %>%
  ggpaired(data = ., x = "Sample", y = "PR_ratio", color = "Sample", 
           line.color = "gray", line.size = 0.4,
           palette = c("#0f77c1", "#00a087"), 
           xlab = "", ylab = "Porphyromonas/Rothia\nlog2(count ratio)", 
           legend.title = "", 
           theme = theme_classic()) + theme(legend.position = "none") +
  stat_compare_means(paired = TRUE)
  
p2<-df.PorRot %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::right_join(meta_pair_RT_ST_select, by = "ID") %>%
  arrange(Sample, Case) %>%
  dplyr::select(PR_ratio, Sample) %>% 
  ggpaired(data = ., x = "Sample", y = "PR_ratio", color = "Sample", 
           line.color = "gray", line.size = 0.4,
           palette = c("#e64b35", "#631879"), 
           xlab = "", ylab = "Porphyromonas/Rothia\nlog2(count ratio)", 
           legend.title = "",  
           theme = theme_classic()) + theme(legend.position = "none") +
  stat_compare_means(paired = TRUE)
ggarrange(p1, p2, labels = c("A", "B"))


#### GT with periodontal parameter####
meta_w_Perio <- read.csv("meta.finalest with Perio.csv") 
meta_w_Perio_GT <- meta_w_Perio %>%
  filter(Sample == "GT") 
df.PorRot_GT <- df.PorRot %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::right_join(meta_w_Perio_GT, by = "ID") 

cor.PorRot_Perio_GT <- psych::corr.test(x = df.PorRot_GT$PR_ratio,
          y = df.PorRot_GT[, c("FMPS", "FMBS" ,"PD" ,"CAL" ,"PD_4" ,"PD_6" ,"M_II" ,"FI_II")],        
          use = "pairwise",method="spearman",adjust="BH", 
          alpha=.05,ci=TRUE,minlength=5,normal=TRUE)

meta_w_Perio_RTST <- meta_w_Perio %>%
  filter(Sample %in% c("RT", "ST")) 

df.PorRot_RTST <- df.PorRot %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::right_join(meta_w_Perio_RTST, by = "ID") 

cor.PorRot_Perio_RTST <- psych::corr.test(x = df.PorRot_RTST$PR_ratio,
                                     y = df.PorRot_RTST[, c("FMPS", "FMBS" ,"PD" ,"CAL" ,"PD_4" ,"PD_6" ,"M_II" ,"FI_II")],        
                                     use = "pairwise",method="spearman",adjust="BH", 
                                     alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
df.PorRot_all <- df.PorRot %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::right_join(meta_w_Perio, by = "ID") 

cor.PorRot_Perio <- psych::corr.test(x = df.PorRot_all$PR_ratio,
                                        y = df.PorRot_all[, c("FMPS", "FMBS" ,"PD" ,"CAL" ,"PD_4" ,"PD_6" ,"M_II" ,"FI_II")],        
                                        use = "pairwise",method="spearman",adjust="BH", 
                                        alpha=.05,ci=TRUE,minlength=5,normal=TRUE)

meta_site <- read_xlsx("Sample vs local defect periodontal parameters new.xlsx", sheet = 1) %>%
  rbind(read_xlsx("Sample vs local defect periodontal parameters new.xlsx", sheet = 2)) %>%
  rbind(read_xlsx("Sample vs local defect periodontal parameters new.xlsx", sheet = 3)) %>%
  rbind(read_xlsx("Sample vs local defect periodontal parameters new.xlsx", sheet = 4)) %>%
  filter(!is.na(True.Name)) 

meta_pisa <- read.csv("sample data tables v2(PISA).csv")

meta.site.final <- sample_meta %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(meta_site, by = "True.Name") 

df.PorRot_site <- df.PorRot %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(meta.site.final, by = "ID") %>%
  dplyr::left_join(meta_pisa[, c("ID", "PISA")], by = "ID")

cor.PorRot_site <- psych::corr.test(x = df.PorRot_site$PR_ratio,
                                     y = df.PorRot_site[, c("Average PPD", "Average CAL", "PISA")],        
                                     use = "pairwise",method="spearman",adjust="BH", 
                                     alpha=.05,ci=TRUE,minlength=5,normal=TRUE)

p.PorRot_PD <- ggscatter(df.PorRot_site, x = "Average PPD", 
          y = "PR_ratio", 
          color = "Sample", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line", 
          add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.sep = "\n"), title = "",
          xlab = "Average PD (mm) for teeth involved", ylab = "Porphyromonas/Rothia\nlog2(count ratio)",
          palette = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  geom_vline(xintercept = 6, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") 
  
df.PorRot_site.pisa <- df.PorRot_site 
cor.PorRot_site.pisa <- psych::corr.test(x = df.PorRot_site.pisa$PR_ratio.x,
                                    y = df.PorRot_site.pisa[, c("PISA")],        
                                    use = "pairwise",method="spearman",adjust="BH", 
                                    alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
p.PorRot_PISA <- ggscatter(df.PorRot_site, x = "PISA", 
                         y = "PR_ratio", 
                         color = "Sample", shape = 21, size = 3, # Points color, shape and size
                         add = "reg.line", 
                         add.params = list(color = "#AD002A", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE, # Add confidence interval
                         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                         cor.coeff.args = list(method = "spearman", label.sep = "\n", label.y = 7.5), title = "",
                         xlab = "PISA (mm2) for teeth involved", ylab = "Porphyromonas/Rothia\nlog2(count ratio)",
                         palette = c("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  geom_vline(xintercept = 200, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") 


ggarrange(p.PorRot_PD, p.PorRot_PISA, labels = c("A", "B"), common.legend = T)
