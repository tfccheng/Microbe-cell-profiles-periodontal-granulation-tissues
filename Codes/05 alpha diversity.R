library(mia)
library(miaViz)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scater)
library(ggsignif)
library(ggpubr)
library(ggsci)
set.seed(123)

setwd("")
#load prefiltered data


tse.silva.ConQuR_BE_remove.new1 <- readRDS("data preparation/tse.silva.ConQuR_BE_remove.new1.rds")
pseq.ConQuR_BE_remove.new1 <- makePhyloseqFromTreeSummarizedExperiment(tse.silva.ConQuR_BE_remove.new1)


tab.alpha <-microbiome::alpha(pseq.ConQuR_BE_remove.new1, index = "all")

colData(tse.silva.ConQuR_BE_remove.new1) <- cbind(colData(tse.silva.ConQuR_BE_remove.new1), tab.alpha) 
tse.silva.ConQuR_BE_remove.new1 <- mia::estimateDiversity(tse.silva.ConQuR_BE_remove.new1, 
                                    abund_values = "counts",
                                    index = "faith", 
                                    name = "faith")


#all alpha indices

ggboxplot(reshape2::melt(data.frame(colData(tse.silva.ConQuR_BE_remove.new1)[,-c(1,3:23)])), 
          x = "variable", y = "value",
          color = "Sample", fill = "Sample", alpha = 0.2,
          add = "jitter") +
  #labs(x ="", y="log10(Counts)") + 
  scale_fill_manual(values =c ("#0f77c1", "#00a087", "#e64b35", "#631879")) +
  scale_color_manual(values = c("#0f77c1", "#00a087", "#e64b35", "#631879"))+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) +
  stat_compare_means(aes(group = Sample), method = "wilcox.test", label = "p.signif") 




#sig indices
library(rstatix)
df <-data.frame(colData(tse.silva.ConQuR_BE_remove.new1)[,-c(1,3:23)])

plots.source2 <- list()
ylabs = c("Observed", "Shannon", "Inverse Simpson", 
          "Berger-Parker index (dbp)")
i<-1
for (n in c("observed", "diversity_shannon",
            "diversity_inverse_simpson", 
            "dominance_dbp")){
  plots.source2[[n]] <- ggboxplot(df,
                                  x = "Sample", y = n,
                                  color = "Sample", fill = "Sample", alpha = 0.2,
                                  #notch = TRUE, 
                                  add = "jitter") +
    labs(y=ylabs[i]) + 
    geom_pwc(method = "wilcox_test", label = "p.format", p.adjust.method = "none", hide.ns = T)+
    # stat_pvalue_manual(stat.test, label = "p.adj")
    theme(axis.title.x = element_blank())+ 
    scale_fill_manual(values =c ("#0f77c1", "#00a087", "#e64b35", "#631879")) +
    scale_color_manual(values = c("#0f77c1", "#00a087", "#e64b35", "#631879"))
  i = i+1
}
ggpubr::ggarrange(plotlist = plots.source2, nrow = 2, ncol = 2, 
                  common.legend = TRUE, legend = "right")













