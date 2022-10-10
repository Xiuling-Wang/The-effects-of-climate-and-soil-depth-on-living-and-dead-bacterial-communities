library(ggplot2)
library(patchwork)
library(tidyverse)
library(introdataviz)
library(ggpubr)
library(ggsci)

alpha <- read.table('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/alpha_diversity.txt', 
                    header = TRUE,check.names = FALSE)
env <- read.csv('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv',
                  header = TRUE,check.names = FALSE)

env <- env %>% 
  filter(env$sample_id%in%alpha$sample_id)


alpha1 <- merge(alpha,env,
               by= "sample_id") %>% 
  select(site,pit,dna_type,layer_n,depth_minus,Shannon,Simpson,Pielous,goods_coverage,depth,observed_species,Chao1,ACE) %>% 
  na.omit() %>% 
  mutate(site = factor(site,
                       levels = c("AZ",
                                  "SG",
                                  "LC",
                                  "NB")),
         dna_type = factor(dna_type,
                         levels = c("iDNA","eDNA"))) 

#1. main figure  ------------------------------------------------------------

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/03 Alpha diversity/shannon.pdf'),
    width=8,height=6)

  ggplot(alpha1, aes(x = site,y = Shannon,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T,lwd=0.2) +
  geom_boxplot(width = .2, alpha = .6, show.legend = FALSE,lwd=0.2) +
  ylim(2,9) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  ggtitle("") +
  ylab("Shannon (H)") + 
  xlab("Sites")+
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  )
  
  
  
  #theme(legend.position = "none")

  dev.off()
  
?geom_boxplot()
# 2. test for main figure between sites ---------------------------------------------------
  my_comparisons = list(c("AZ","SG"),c("AZ","LC"),
                        c("AZ","NB"),c("SG","LC"),
                        c( "SG","NB"),c("LC","NB"))
  
  pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/03 Alpha diversity/i_shannon.pdf'),
      width=8,height=6) 
  
  ggviolin(alpha1 %>% 
             filter(dna_type=="iDNA"),x="site",y="Shannon",fill = "site",
           palette = c('lancet'),
           add = "boxplot",
           add.params = list(fill="white"),
           trim = T)+
    stat_compare_means(comparisons = my_comparisons,label="p.signif")+
    ggtitle("idna") +
    ylab("Shannon (H)") + 
    xlab("Sites")
dev.off()

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/03 Alpha diversity/e_shannon.pdf'),
    width=8,height=6) 

ggviolin(alpha1 %>% 
           filter(dna_type=="eDNA"),x="site",y="Shannon",fill = "site",
         palette = c('lancet'),
         add = "boxplot",
         add.params = list(fill="white"),
         trim = T)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")+
  ggtitle("edna") +
  ylab("Shannon (H)") + 
  xlab("Sites")
dev.off()


# 3. figure in appendix (3 layer)--------------------------------------------------------

shannon_layer1 <- ggplot(alpha1 %>% 
         filter(layer_n==1), aes(x = site,y = Shannon,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE) +
  ylim(2,8) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  theme_minimal()+
  ggtitle("") +
  ylab("Shannon (H) 0-40") + 
  xlab("Sites")+
  theme(legend.position = "none")

shannon_layer2 <- ggplot(alpha1 %>% 
         filter(layer_n==2), aes(x = site,y = Shannon,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE) +
  ylim(2,8) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  theme_minimal()+
  ggtitle("") +
  ylab("Shannon (H) 40-120") + 
  xlab("Sites")+
  theme(legend.position = "none")

shannon_layer3 <- ggplot(alpha1 %>% 
         filter(layer_n==3), aes(x = site,y = Shannon,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE) +
  ylim(2,8) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  theme_minimal()+
  ggtitle("") +
  ylab("Shannon (H) 120-200") + 
  xlab("Sites")+
  theme(legend.position = "none")

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/03 Alpha diversity/shannon_differentlayer.pdf'),
    width=12,height=6) 
shannon_layer1+shannon_layer2+shannon_layer3
dev.off()


# test for appendix -------------------------------------------------------
layeri1 <- ggviolin(alpha1 %>% 
           filter(dna_type=="iDNA",layer_n==1),x="site",y="Shannon",fill = "site",
         palette = c('lancet'),
         add = "boxplot",
         add.params = list(fill="white"),
         trim = T)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")+
theme_minimal()+
  ggtitle("idna") +
  ylab("Shannon (H) 0-40") + 
  xlab("Sites")+
  theme(legend.position = "none")


layeri2 <- ggviolin(alpha1 %>% 
                     filter(dna_type=="iDNA",layer_n==2),x="site",y="Shannon",fill = "site",
                   palette = c('lancet'),
                   add = "boxplot",
                   add.params = list(fill="white"),
                   trim = T)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")+
  theme_minimal()+
  ggtitle("idna") +
  ylab("Shannon (H) 40-120") + 
  xlab("Sites")+
  theme(legend.position = "none")

layeri3 <- ggviolin(alpha1 %>% 
                     filter(dna_type=="iDNA",layer_n==3),x="site",y="Shannon",fill = "site",
                   palette = c('lancet'),
                   add = "boxplot",
                   add.params = list(fill="white"),
                   trim = T)+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")+
  theme_minimal()+
  ggtitle("idna") +
  ylab("Shannon (H) 120-200") + 
  xlab("Sites")+
  theme(legend.position = "none")


pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/03 Alpha diversity/i3layer_shannon.pdf'),
    width=8,height=6) 
layeri1+layeri2+layeri3
dev.off()










pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/03 Alpha diversity/i3layer_shannon.pdf'),
    width=8,height=6) 

# Pielous -------------------------------------------------------------------

pielous_layer1 <- ggplot(alpha1 %>% 
                           filter(layer_n==1), aes(x = site,y = Pielous,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE) +
  ylim(0.5,1) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  theme_minimal()+
  ggtitle("") +
  ylab("Pielou (J)") + 
  xlab("Sites")+
  theme(legend.position = "none")


pielous_layer2 <- ggplot(alpha1 %>% 
                           filter(layer_n==2), aes(x = site,y = Pielous,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE) +
  ylim(0.5,1) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  theme_minimal()+
  ggtitle("") +
  ylab("Pielou (J)") + 
  xlab("Sites")+
  theme(legend.position = "none")

pielous_layer3 <- ggplot(alpha1 %>% 
                           filter(layer_n==3), aes(x = site,y = Pielous,fill = dna_type))+
  geom_split_violin(alpha = .4, trim = T) +
  geom_boxplot(width = .1, alpha = .6, show.legend = FALSE) +
  ylim(0.5,1) +
  stat_compare_means(aes(group = dna_type), label = "p.signif")+
  scale_fill_manual(values = c("dodgerblue2", "darkorange"), name = "DNA Types") +
  theme_minimal()+
  ggtitle("") +
  ylab("Pielou (J)") + 
  xlab("Sites")+
  theme(legend.position = "none")
