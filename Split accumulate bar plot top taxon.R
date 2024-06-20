library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(forcats)
library(stringr)
library(RColorBrewer)
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1,
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  header = TRUE)
# clean raw data ----------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id %in% colnames(asv))

filt1 <- colnames(asv) %in% env1$sample_id
asv1 <- asv[,filt1]
asv1$asv_id <- rownames(asv1)

tax1 <- tax %>%
  filter(tax$asv_id %in% asv1$asv_id)

# 01. top Phylum ----------------------------------------------------------
top_phylum <- asv1 %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Phylum) %>%
  summarise(PhylumAbun = sum(Abun)) %>%
  mutate(RelAbun = 100*PhylumAbun/sum(PhylumAbun)) %>%
  ungroup() %>%
  arrange(desc(RelAbun))

# 01. top Phylum ----------------------------------------------------------


top11phylum <- top_phylum$Phylum[1:13]
unique(top_phylum$Phylum)
# tidy plot data ----------------------------------------------------------
plot_data1 <- asv1 %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(site,
           dna_type,
           depths_cm_re,
           Phylum) %>%
  summarise(PhylumAbun = sum(Abun)) %>%
  ungroup %>%
  group_by(site,dna_type) %>%
  mutate(RelAbun = 100*PhylumAbun/sum(PhylumAbun)) %>%
  ungroup()
# calculate other phylum(except top) -------------------------------------
plot_data2 <- plot_data1 %>%
  filter(!Phylum %in% top11phylum) %>%
  group_by(site,
           depths_cm_re,
           dna_type) %>%
  summarise(Other = sum(RelAbun)) %>%
  ungroup()
# rebuild a new data frame for Other Phylum
plot_data3 <- plot_data2 %>%
  rename(RelAbun = Other)
plot_data3$Phylum <- "Other"
plot_data3$PhylumAbun <- 1

plot_data4 <- rbind(plot_data1,
                    plot_data3)
# plot bar plot -----------------------------------------------------------
top12phylum <- c(top11phylum,'Other')

plot <- plot_data4 %>% 
  filter(Phylum %in% top12phylum)
rev(top12phylum)

plot$Phylum <- factor(plot$Phylum,
                     levels = c("p__NA","Other","p__GAL15","p__Myxococcota","p__Methylomirabilota",
                                "p__Bacteroidota","p__Gemmatimonadota",
                                 "p__Firmicutes","p__Planctomycetota","p__Verrucomicrobiota","p__Acidobacteriota",
                                 "p__Chloroflexi","p__Proteobacteria","p__Actinobacteriota"),
                     labels = c("Unassigned","Other Phyla","GAL15","Myxococcota","Methylomirabilota",
                                "Bacteroidota","Gemmatimonadota",
                                "Firmicutes","Planctomycetota","Verrucomicrobiota","Acidobacteriota",
                                "Chloroflexi","Proteobacteria","Actinobacteriota"))

plot$depths_cm_re <- factor(plot$depths_cm_re,
                            levels = c("01_180_200","02_160_180", "03_140_160","04_120_140",
                                       "05_100_120","06_80_100","07_60_80", "08_40_60",
                                       "09_20_40","10_10_20","11_5_10","12_0_5"),
                            labels = c("180-200","160-180","140-160","120-140","100-120",
                                       "80-100","60-80","40-60","20-40","10-20","5-10","0-5"))
plot$site <- factor(plot$site,
                    levels = c("AZ",
                               "SG",
                               "LC",
                               "NB"))

paired_col <- brewer.pal(12, "Paired")
new_col <- c("white","darkgray",paired_col)
# AZ --------------------------------------------------------------------
p1 <- ggplot(plot %>%
               filter(site =='AZ', 
                      dna_type =='iDNA'),
             aes(x = RelAbun, 
                 y = depths_cm_re, 
                 fill = Phylum)) +
  geom_bar(stat = "identity",
           position = "fill", 
           color = "black", 
           width = 0.7,
           size = 0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",
                  label=function(x){paste0(100*x)})+
  background_grid(major = "x", 
                  minor = "none")+
  ggtitle("AZ iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))
p1

p2 <- ggplot(plot %>%
               filter(site =='AZ', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",
                     label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("AZ eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))

p2
# SG -------------------------------------------------------------------
p3 <- ggplot(plot %>%
               filter(site =='SG', dna_type=='iDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("SG iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))
p3

p4 <- ggplot(plot %>%
               filter(site =='SG', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("SG eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))
p4
#LC ----------------------------------------------------------------------
p5 <- ggplot(plot %>%
               filter(site =='LC', dna_type=='iDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("LC iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))

p5

p6 <- ggplot(plot %>%
               filter(site =='LC', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("LC eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))

p6
#NA ------------------------------------------------------------------
p7 <- ggplot(plot %>%
               filter(site =='NB', dna_type=='iDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("NA iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))
p7

p8 <- ggplot(plot %>%
               filter(site =='NB', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Phylum)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("NA eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        #legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  guides(fill = guide_legend(reverse=TRUE))
p8

all <- (p1+p2|p3+p4)/
  (p5+p6|p7+p8)
all
ggsave(all, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/04 Species structure barplot/top12_Phyla.pdf",
       width = 60, 
       height = 40, 
       units = "cm")
# legend ------------------------------------------------------------------
plegend <- ggplot(plot %>%
                    filter(site =='NB', 
                           dna_type =='eDNA'),
                  aes(x = RelAbun, 
                      y = depths_cm_re, 
                      fill=Phylum)) +
  geom_bar(stat="identity",
           position="fill", 
           color="black", 
           width=0.7,
           size=0.15)+
  scale_fill_manual(values = new_col)+
  theme(legend.position="bottom")+
  guides(fill = guide_legend(reverse=TRUE))

plegend
ggsave(legend, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/legendsV.pdf",
       width = 40, 
       height = 30, 
       units = "cm")

# 02. top Class -----------------------------------------------------------
# Load data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  header = TRUE)
# clean raw data ----------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id%in%colnames(asv))
tax1 <- tax %>% 
  filter(tax$asv_id%in%asv$asv_id)

top_class <- asv %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Class) %>%
  summarise(ClassAbun = sum(Abun)) %>%
  mutate(RelAbun = 100*ClassAbun/sum(ClassAbun)) %>%
  ungroup() %>%
  arrange(desc(RelAbun))

top10class <- top_class$Class[1:11]
# tidy plot data ----------------------------------------------------------
plot_data1 <- asv %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(site,
           dna_type,
           depths_cm_re,
           Class) %>%
  summarise(ClassAbun = sum(Abun)) %>%
  ungroup %>%
  group_by(site,dna_type) %>%
  mutate(RelAbun = 100*ClassAbun/sum(ClassAbun)) %>%
  ungroup()
# calculate other class(except top) -------------------------------------
plot_data2 <- plot_data1 %>%
  filter(!Class%in%top10class) %>%
  group_by(site,
           depths_cm_re,
           dna_type) %>%
  summarise(Other = sum(RelAbun)) %>%
  ungroup()
# rebuild a new data frame for Other Class
plot_data3 <- plot_data2 %>%
  rename(RelAbun = Other)
plot_data3$Class <- "Other"
plot_data3$ClassAbun <- 1

plot_data4 <- rbind(plot_data1,
                    plot_data3)
# plot bar plot -----------------------------------------------------------
top11class <- c(top10class,'Other')

plot <- plot_data4 %>% 
  filter(Class %in% top11class)

plot$Class <- factor(plot$Class,
                     levels = c("Other","c__Vicinamibacteria","c__NA","c__Bacilli",            
                                "c__Acidobacteriae","c__Planctomycetes","c__Verrucomicrobiae",  
                                "c__AD3","c__Gammaproteobacteria","c__Alphaproteobacteria",
                                "c__Actinobacteria","c__Thermoleophilia"),
                     labels = c("Other","Vicinamibacteria","unassigned","Bacilli",            
                                "Acidobacteriae","Planctomycetes","Verrucomicrobiae",  
                                "AD3","Gammaproteobacteria","Alphaproteobacteria",
                                "Actinobacteria","Thermoleophilia"))

plot$depths_cm_re <- factor(plot$depths_cm_re,
                         levels = c("01_180_200","02_160_180", "03_140_160","04_120_140",
                                    "05_100_120","06_80_100","07_60_80", "08_40_60",
                                    "09_20_40","10_10_20","11_5_10","12_0_5"),
                         labels = c("180-200","160-180","140-160","120-140","100-120",
                                    "80-100","60-80","40-60","20-40","10-20","5-10","0-5"))
plot$site <- factor(plot$site,
                    levels = c("AZ",
                               "SG",
                               "LC",
                               "NB"))
# AZ --------------------------------------------------------------------
p1 <- ggplot(plot %>%
               filter(site =='AZ', 
                      dna_type=='iDNA'),
             aes(x=RelAbun, 
                 y=depths_cm_re, 
                 fill=Class)) +
  geom_bar(stat="identity",
           position="fill", 
           color="black", 
           width=0.7,
           size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",
                  label=function(x){paste0(100*x)})+
  background_grid(major = "x", 
                  minor = "none")+
  ggtitle("AZ iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p1

p2 <- ggplot(plot %>%
               filter(site =='AZ', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("AZ eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p2

az <- p1+pcraz+p2

ggsave(az, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/accumlatebarAZ.pdf",
       width = 15, 
       height = 20, 
       units = "cm")

# SG -------------------------------------------------------------------
p3 <- ggplot(plot %>%
               filter(site =='SG', dna_type=='iDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("SG iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p3

p4 <- ggplot(plot %>%
               filter(site =='SG', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("SG eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p4

sg <- p3+pcrsg+p4
sg

# LC ----------------------------------------------------------------------
p5 <- ggplot(plot %>%
               filter(site =='LC', dna_type=='iDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("LC iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p5

p6 <- ggplot(plot %>%
               filter(site =='LC', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("LC eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p6

lc <- p5+pcrlc+p6
lc
ggsave(lc, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/accumlatebarLC.pdf",
       width = 15, 
       height = 20, 
       units = "cm")

# NA ------------------------------------------------------------------

p7 <- ggplot(plot %>%
               filter(site =='NB', dna_type=='iDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_reverse(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("NA iDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size=16,vjust=1,hjust = 1),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
p7

p8 <- ggplot(plot %>%
               filter(site =='NB', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("NA eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.position="none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

p8
nb <- p7+pcrna+p8
nb
ggsave(nb, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/accumlatebarNA.pdf",
       width = 15, 
       height = 20, 
       units = "cm")

all <- (p1+pcraz+p2|p3+pcrsg+p4)/
  (p5+pcrlc+p6|p7+pcrna+p8)
all
ggsave(all, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/accumlatebarallclass.pdf",
       width = 60, 
       height = 40, 
       units = "cm")

#legend
plegend <- ggplot(plot %>%
               filter(site =='NB', dna_type=='eDNA'),
             aes(x=RelAbun, y=depths_cm_re, fill=Class)) +
  geom_bar(stat="identity",position="fill", 
           color="black", width=0.7,size=0.15)+
  scale_fill_manual(values = new_col)+
  scale_x_continuous(position="top",label=function(x){paste0(100*x)})+
  background_grid(major = "x", minor = "none")+
  ggtitle("NA eDNA") +
  xlab("Relative Abundance (%)") + 
  ylab("Depth (cm)")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=16, angle=45),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        legend.position="bottom",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
plegend

legend123 <- plegend+pcrlegend
legend123
ggsave(legend123, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/legendsV.pdf",
       width = 40, 
       height = 30, 
       units = "cm")

#03. top genus --------------------------------------------------------
top_genus <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = 100*GenusAbun/sum(GenusAbun)) %>%
  ungroup()%>%
  arrange(desc(RelAbun))

top10genus <- top_genus$Genus[1:11]
# otherclass <- top_class$Class[12:141]
# filt1=top_class$RelAbun
# filt2=filt1[22:141]
# otherclasssum <- sum(filt2)
# tidy plot data ----------------------------------------------------------
plot_data1 <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun')%>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(site,dna_type,depth,Genus) %>%
  summarise(GenusAbun = sum(Abun)) %>%
  ungroup %>%
  group_by(site, dna_type)%>%
  mutate(RelAbun = 100*GenusAbun/sum(GenusAbun)) %>%
  ungroup()

# calculate other class(except top20) -------------------------------------
plot_data2 <- plot_data1%>%
  filter(!Genus%in%top10genus)%>%
  group_by(site,depth,dna_type) %>%
  summarise(Other = sum(RelAbun)) %>%
  ungroup()

# rebuild a new data frame for Other Class
plot_data3 <- plot_data2%>%
  rename(RelAbun= Other)

plot_data3$Genus <- "Other"
plot_data3$GenusAbun <- 1

plot_data4 <- rbind(plot_data1,plot_data3)
# plot bar plot -----------------------------------------------------------
top11genus <- c(top10genus,'Other')
top10nonagenus <- top10genus[2:11]

plot <- plot_data4 %>% 
  filter(Genus%in%
           top10nonagenus)

plot$Genus <- factor(plot$Genus,
                     levels = c("g__Solirubrobacter","g__Delftia",
                                "g__Nocardioides","g__Conexibacter",        
                                "g__Rubrobacter","g__RB41",
                                "g__Halomonas","g__Gaiella",             
                                "g__Streptomyces","g__Candidatus_Udaeobacter"),
                     labels = c("Solirubrobacter","Delftia",
                                "Nocardioides","Conexibacter",        
                                "Rubrobacter","RB41",
                                "Halomonas","Gaiella",             
                                "Streptomyces","Candidatus_Udaeobacter"))

#AZ -------------------------------------------------------------------
p01 <- ggplot(plot %>% 
               filter(site =='AZ', dna_type=='iDNA'),
             aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("AZ iDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))
p01
p02 <- ggplot(plot %>% 
               filter(site =='AZ', dna_type=='eDNA'),
             aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_continuous()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("AZ eDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

#SG -----------------------------------------
p03 <- ggplot(plot %>% 
                filter(site =='SG', dna_type=='iDNA'),
              aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("SG iDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))

p04 <- ggplot(plot %>% 
                filter(site =='SG', dna_type=='eDNA'),
              aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_continuous()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("SG eDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")
#LC --------------------------------------------------------------------
p05 <- ggplot(plot %>% 
                filter(site =='LC', dna_type=='iDNA'),
              aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("LC iDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))

p06 <- ggplot(plot %>% 
                filter(site =='LC', dna_type=='eDNA'),
              aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_continuous()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("LC eDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

#NA -------------------------------------------------------------------
p07 <- ggplot(plot %>% 
                filter(site =='NB', dna_type=='iDNA'),
              aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("NA iDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5))

p08 <- ggplot(plot %>% 
                filter(site =='NB', dna_type=='eDNA'),
              aes(y=RelAbun, x=depth, fill=Genus)) +
  geom_bar(stat="identity",position="stack", color="black", width=4.5,size=0.15)+
  scale_fill_manual(values = new_col)+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE)+
  scale_x_reverse()+
  scale_y_continuous()+
  theme(panel.background = element_rect(fill='transparent'))+
  background_grid(major = "xy", minor = "none")+
  theme_minimal()+
  ggtitle("NA eDNA") +
  ylab("Relative Abundance") + 
  xlab("Depth (cm)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

allgenus <- (p01+p02+p03+p04)|(p05+p06+p07+p08)

allgenus <- plot_grid(p01+p02,p03+p04,p05+p06,p07+p08,
                     labels=c("a","b","c","d"), ncol = 2,align = "h")
save_plot("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/56.pdf",
          allgenus,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)
          
ggsave(allgenus, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/accumlatebarallgenus.pdf",
       width = 30, height = 40, units = "cm")

# 计算不同样地top phylum/class/genus --------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1,header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  header = TRUE)
# clean raw data ----------------------------------------------------------
env1 <- env %>%
  filter(site=="LC",#replace
         dna_typen==0,
         env$sample_id%in%colnames(asv))

filt1 <- colnames(asv)%in%env1$sample_id

asv1 <- asv[,filt1]
asv1$asv_id <- rownames(asv1)

tax1 <- tax %>%
  filter(tax$asv_id%in%asv1$asv_id)
# get top 10 CLASS list first ---------------------------------------------
top_class <- asv1 %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus) %>%#replace
  summarise(ClassAbun = sum(Abun)) %>%
  mutate(RelAbun = 100*ClassAbun/sum(ClassAbun)) %>%
  ungroup() %>%
  arrange(desc(RelAbun))

top10class <- top_class$Class[1:11]