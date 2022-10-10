library(tidyverse)
library(magrittr)
library(VennDiagram)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
# Load data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header = TRUE,check.names = FALSE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                header = TRUE,check.names = FALSE)
#
env <- env %>% 
  filter(depth<200)
all_depths <- unique(env$depths_cm)
depth1 <- all_depths[1]
data_list <- list()
# 1. AZ  all depth----------------------------------------------------------------------
for(depth1 in all_depths){
  tmp <- env %>% 
    filter(site=="AZ") %>% 
    filter(depths_cm==depth1) %>% 
    select(depths_cm,dna_type,sample_id)
  
  eDNA_sample <- tmp %>% 
    filter(dna_type=='eDNA') %>% 
    pull(sample_id)
  
  iDNA_sample <- tmp  %>% 
    filter(dna_type=='iDNA') %>% 
    pull(sample_id)
  
  eDNA_sample <- eDNA_sample[eDNA_sample %in% names(asv)]
  
  iDNA_sample <- iDNA_sample[iDNA_sample %in% names(asv)]
  
  eDNA_asv <- asv$asv_id[asv[eDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  iDNA_asv <- asv$asv_id[asv[iDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  intersect_num <- length(intersect(eDNA_asv, iDNA_asv))
  
  iDNA_num <- length(iDNA_asv)
  eDNA_num <- length(eDNA_asv)
  
  data_list[[depth1]] <- c(iDNA = iDNA_num, 
                           eDNA = eDNA_num,
                           iDNAandeDNA = intersect_num)
}
data_list
names(data_list)
a <-  Reduce(rbind, data_list)
rownames(a) = names(data_list)
a
a2 <- sweep(a, 1, rowSums(a),'/')*100 
a2 <- as.data.frame(a2)
depth_list <- rownames(a2)
a3 <- cbind(depth_list,a2)

a3long <- a3 %>% 
  pivot_longer(cols = iDNA:iDNAandeDNA,
               names_to = "area",
               values_to = "percent")
a3long$area <- factor(a3long$area,
                      levels = c("eDNA","iDNAandeDNA","iDNA"),
                      labels = c("eDNA","Share","iDNA"))
# c("12_0_5","11_5_10","10_10_20",
#   "09_20_40","08_40_60","07_60_80",
#   "06_80_100","05_100_120","04_120_140",
#   "03_140_160","02_160_180","01_180_200")
#c("0_5","5_10","10_20","20_40","40_60","60_80","80_100","100_120","120_140","140_160","160_180","180_200")
a3long$depth_list <- factor(a3long$depth_list,
                          levels = rev(c("12_0_5","11_5_10","10_10_20",
                                         "09_20_40","08_40_60","07_60_80",
                                         "06_80_100","05_100_120","04_120_140",
                                         "03_140_160","02_160_180","01_180_200")),
                          labels = c("180-200","160-180","140-160","120-140",
                                     "100-120","80-100","60-80","40-60",
                                     "20-40","10-20","5-10","0-5"))

AZ <- ggplot(data=a3long,
             aes(x=percent,
                 y=depth_list,
                 fill=area))+
  geom_col(position='stack', 
           width = 0.6)+
  scale_x_continuous(position="top")+
  scale_fill_manual(values = c("darkorange","gray","dodgerblue2"), 
                    name = "ASVs distribution") +
  ggtitle("AZ") +
  ylab("Depth (cm)") + 
  xlab("Relative Abundance of ASVs (%)")+
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.title.x = element_text(size = 8)
  )

AZ
# 2.  SG in all depth -----------------------------------------------------
for(depth1 in all_depths){
  tmp <- env %>% 
    filter(site=="SG") %>% 
    filter(depths_cm==depth1) %>% 
    select(depths_cm,dna_type,sample_id)
  
  eDNA_sample <- tmp %>% 
    filter(dna_type=='eDNA') %>% 
    pull(sample_id)
  
  iDNA_sample <- tmp  %>% 
    filter(dna_type=='iDNA') %>% 
    pull(sample_id)
  
  eDNA_sample <- eDNA_sample[eDNA_sample %in% names(asv)]
  
  iDNA_sample <- iDNA_sample[iDNA_sample %in% names(asv)]
  
  eDNA_asv <- asv$asv_id[asv[eDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  iDNA_asv <- asv$asv_id[asv[iDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  intersect_num <- length(intersect(eDNA_asv, iDNA_asv))
  
  iDNA_num <- length(iDNA_asv)
  eDNA_num <- length(eDNA_asv)
  
  data_list[[depth1]] <- c(iDNA = iDNA_num, 
                           eDNA = eDNA_num,
                           iDNAandeDNA = intersect_num)
}
data_list
names(data_list)
a <-  Reduce(rbind, data_list)
rownames(a) = names(data_list)
a
a2 <- sweep(a, 1, rowSums(a),'/')*100 
a2 <- as.data.frame(a2)
depth_list <- rownames(a2)
a3 <- cbind(depth_list,a2)

a3long <- a3 %>% 
  pivot_longer(cols = iDNA:iDNAandeDNA,
               names_to = "area",
               values_to = "percent")
a3long$area <- factor(a3long$area,
                      levels = c("eDNA","iDNAandeDNA","iDNA"),
                      labels = c("eDNA","Share","iDNA"))
# c("12_0_5","11_5_10","10_10_20",
#   "09_20_40","08_40_60","07_60_80",
#   "06_80_100","05_100_120","04_120_140",
#   "03_140_160","02_160_180","01_180_200")
#c("0_5","5_10","10_20","20_40","40_60","60_80","80_100","100_120","120_140","140_160","160_180","180_200")
a3long$depth_list <- factor(a3long$depth_list,
                            levels = rev(c("12_0_5","11_5_10","10_10_20",
                                           "09_20_40","08_40_60","07_60_80",
                                           "06_80_100","05_100_120","04_120_140",
                                           "03_140_160","02_160_180","01_180_200")),
                            labels = c("180-200","160-180","140-160","120-140",
                                       "100-120","80-100","60-80","40-60",
                                       "20-40","10-20","5-10","0-5"))

SG <- ggplot(data=a3long,
             aes(x=percent,
                 y=depth_list,
                 fill=area))+
  geom_col(position='stack', 
           width = 0.6)+
  scale_x_continuous(position="top")+
  scale_fill_manual(values = c("darkorange","gray","dodgerblue2"), 
                    name = "ASVs distribution") +
  ggtitle("SG") +
  ylab("Depth (cm)") + 
  xlab("Relative Abundance of ASVs (%)")+
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'),
        
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  )
SG

# 3. LC at all depth ------------------------------------------------------
for(depth1 in all_depths){
  tmp <- env %>% 
    filter(site=="LC") %>% 
    filter(depths_cm==depth1) %>% 
    select(depths_cm,dna_type,sample_id)
  
  eDNA_sample <- tmp %>% 
    filter(dna_type=='eDNA') %>% 
    pull(sample_id)
  
  iDNA_sample <- tmp  %>% 
    filter(dna_type=='iDNA') %>% 
    pull(sample_id)
  
  eDNA_sample <- eDNA_sample[eDNA_sample %in% names(asv)]
  
  iDNA_sample <- iDNA_sample[iDNA_sample %in% names(asv)]
  
  eDNA_asv <- asv$asv_id[asv[eDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  iDNA_asv <- asv$asv_id[asv[iDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  intersect_num <- length(intersect(eDNA_asv, iDNA_asv))
  
  iDNA_num <- length(iDNA_asv)
  eDNA_num <- length(eDNA_asv)
  
  data_list[[depth1]] <- c(iDNA = iDNA_num, 
                           eDNA = eDNA_num,
                           iDNAandeDNA = intersect_num)
}
data_list
names(data_list)
a <-  Reduce(rbind, data_list)
rownames(a) = names(data_list)
a
a2 <- sweep(a, 1, rowSums(a),'/')*100 
a2 <- as.data.frame(a2)
depth_list <- rownames(a2)
a3 <- cbind(depth_list,a2)

a3long <- a3 %>% 
  pivot_longer(cols = iDNA:iDNAandeDNA,
               names_to = "area",
               values_to = "percent")
a3long$area <- factor(a3long$area,
                      levels = c("eDNA","iDNAandeDNA","iDNA"),
                      labels = c("eDNA","Share","iDNA"))
# c("12_0_5","11_5_10","10_10_20",
#   "09_20_40","08_40_60","07_60_80",
#   "06_80_100","05_100_120","04_120_140",
#   "03_140_160","02_160_180","01_180_200")
#c("0_5","5_10","10_20","20_40","40_60","60_80","80_100","100_120","120_140","140_160","160_180","180_200")
a3long$depth_list <- factor(a3long$depth_list,
                            levels = rev(c("12_0_5","11_5_10","10_10_20",
                                           "09_20_40","08_40_60","07_60_80",
                                           "06_80_100","05_100_120","04_120_140",
                                           "03_140_160","02_160_180","01_180_200")),
                            labels = c("180-200","160-180","140-160","120-140",
                                       "100-120","80-100","60-80","40-60",
                                       "20-40","10-20","5-10","0-5"))

LC <- ggplot(data=a3long,
             aes(x=percent,
                 y=depth_list,
                 fill=area))+
  geom_col(position='stack', 
           width = 0.6)+
  scale_x_continuous(position="top")+
  scale_fill_manual(values = c("darkorange","gray","dodgerblue2"), 
                    name = "ASVs distribution") +
  ggtitle("LC") +
  ylab("Depth (cm)") + 
  xlab("Relative Abundance of ASVs (%)")+
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'),
        
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  )
LC

# 4. NA all depth ---------------------------------------------------------

for(depth1 in all_depths){
  tmp <- env %>% 
    filter(site=="NB") %>% 
    filter(depths_cm==depth1) %>% 
    select(depths_cm,dna_type,sample_id)
  
  eDNA_sample <- tmp %>% 
    filter(dna_type=='eDNA') %>% 
    pull(sample_id)
  
  iDNA_sample <- tmp  %>% 
    filter(dna_type=='iDNA') %>% 
    pull(sample_id)
  
  eDNA_sample <- eDNA_sample[eDNA_sample %in% names(asv)]
  
  iDNA_sample <- iDNA_sample[iDNA_sample %in% names(asv)]
  
  eDNA_asv <- asv$asv_id[asv[eDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  iDNA_asv <- asv$asv_id[asv[iDNA_sample] %>% 
                           apply(1, sum) > 0]
  
  intersect_num <- length(intersect(eDNA_asv, iDNA_asv))
  
  iDNA_num <- length(iDNA_asv)
  eDNA_num <- length(eDNA_asv)
  
  data_list[[depth1]] <- c(iDNA = iDNA_num, 
                           eDNA = eDNA_num,
                           iDNAandeDNA = intersect_num)
}
data_list
names(data_list)
a <-  Reduce(rbind, data_list)
rownames(a) = names(data_list)
a
a2 <- sweep(a, 1, rowSums(a),'/')*100 
a2 <- as.data.frame(a2)
depth_list <- rownames(a2)
a3 <- cbind(depth_list,a2)

a3long <- a3 %>% 
  pivot_longer(cols = iDNA:iDNAandeDNA,
               names_to = "area",
               values_to = "percent")
a3long$area <- factor(a3long$area,
                      levels = c("eDNA","iDNAandeDNA","iDNA"),
                      labels = c("eDNA","Share","iDNA"))
# c("12_0_5","11_5_10","10_10_20",
#   "09_20_40","08_40_60","07_60_80",
#   "06_80_100","05_100_120","04_120_140",
#   "03_140_160","02_160_180","01_180_200")
#c("0_5","5_10","10_20","20_40","40_60","60_80","80_100","100_120","120_140","140_160","160_180","180_200")
a3long$depth_list <- factor(a3long$depth_list,
                            levels = rev(c("12_0_5","11_5_10","10_10_20",
                                           "09_20_40","08_40_60","07_60_80",
                                           "06_80_100","05_100_120","04_120_140",
                                           "03_140_160","02_160_180","01_180_200")),
                            labels = c("180-200","160-180","140-160","120-140",
                                       "100-120","80-100","60-80","40-60",
                                       "20-40","10-20","5-10","0-5"))

NB <- ggplot(data=a3long,
             aes(x=percent,
                 y=depth_list,
                 fill=area))+
  geom_col(position='stack', 
           width = 0.6)+
  scale_x_continuous(position="top")+
  scale_fill_manual(values = c("darkorange","gray","dodgerblue2"), 
                    name = "ASVs distribution") +
  ggtitle("NA") +
  ylab("Depth (cm)") + 
  xlab("Relative Abundance of ASVs (%)")+
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'),
        
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
        #legend.key = element_rect(fill = "transparent")
  )
NB

fourvenn <- AZ+SG+LC+NB+plot_layout(ncol=4)
fourvenn

ggsave(fourvenn, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/00 Veendiagram/barvennall.pdf",
       width = 24, height = 12, units = "cm")

legend <- ggplot(data=a3long,
             aes(x=percent,
                 y=depth_list,
                 fill=area))+
  geom_col(position='stack', 
           width = 0.6)+
  scale_x_continuous(position="top")+
  scale_fill_manual(values = c("darkorange","gray","dodgerblue2"), 
                    name = "ASVs distribution") +
  ggtitle("NA") +
  ylab("Depth (cm)") + 
  xlab("Relative Abundance of ASVs (%)")+
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'),
        
        plot.title = element_text(hjust = 0.5),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.key = element_rect(fill = "transparent")
  )
legend

ggsave(legend, filename = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/00 Veendiagram/legend.pdf",
       width = 24, height = 12, units = "cm")
