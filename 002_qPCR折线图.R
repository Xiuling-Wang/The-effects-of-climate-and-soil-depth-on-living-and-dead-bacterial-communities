library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(scales)
env<- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                       header = TRUE)
env1<- env %>%
  filter(depth<200)

genecopydata <- dplyr::select(env1,
                              site,
                              dna_type,
                              depths_cm,
                              genecopy_r1,
                              genecopy_r2,
                              genecopy_r3) %>% 
  na.omit()

genelong <- pivot_longer(genecopydata,
                         cols= c(genecopy_r1,genecopy_r2,genecopy_r3),
                         names_to = "genecopy_r",values_to = "genecopy")
head(genelong,n=18)
table(genelong$site)

genelong %>% 
  group_by(site,dna_type,depths_cm) %>% 
  summarise(mean_value=mean(genecopy)) %>% 
  mutate(site = factor(site,
                     levels = c("AZ",
                                "SG",
                                "LC",
                                "NB")),
         depths_cm = factor(depths_cm,
                            levels = c("01_180_200","02_160_180", "03_140_160","04_120_140",
                                       "05_100_120","06_80_100","07_60_80", "08_40_60",
                                       "09_20_40","10_10_20","11_5_10","12_0_5"),
                            labels = c("180-200","160-180","140-160","120-140","100-120",
                                       "80-100","60-80","40-60","20-40","10-20","5-10","0-5"))) -> df01

genelong %>% 
  mutate(site = factor(site,
                       levels = c("AZ",
                                  "SG",
                                  "LC",
                                  "NB")),
         depths_cm = factor(depths_cm,
                            levels = c("01_180_200","02_160_180", "03_140_160","04_120_140",
                                       "05_100_120","06_80_100","07_60_80", "08_40_60",
                                       "09_20_40","10_10_20","11_5_10","12_0_5"),
                            labels = c("180-200","160-180","140-160","120-140","100-120",
                                       "80-100","60-80","40-60","20-40","10-20","5-10","0-5"))) -> df02

# 1. all sites ------------------------------------------------------------
pcr <- ggplot(data=df01,
              aes(x=mean_value,
                  y=depths_cm))+
  geom_point(data = df02,
             aes(x= genecopy, 
                 y=depths_cm,
                 color=dna_type),
             size=3,stat = "identity")+
  geom_path(aes(group = dna_type,
                linetype = dna_type))+
  scale_x_continuous(position = "top",trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "t")+
  facet_wrap(~site,
             scales = "free_x",
             strip.position = "b")+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_color_manual(values=c('darkorange','dodgerblue2'))+
  ggtitle("QPCR") +
  ylab("Depth (cm)") + 
  xlab("genecopy number/g")+
  theme_bw()+
  theme(legend.position = "bottom")

pcr
# won't add error bar in the plot due to 10^x data does fit well

# 2. for eatch site ----------------------------------------------------------

# 3. az --------------------------------------------------------------------
pcraz <- ggplot(data=df01%>% 
         filter(site=="AZ"),
         aes(x=mean_value,
             y=depths_cm))+
  theme_bw()+
  geom_point(data = df02%>% 
               filter(site=="AZ"),aes(x= genecopy, y=depths_cm,color=dna_type),size=3, stat = "identity")+
  geom_path(aes(group=dna_type,linetype= dna_type))+
  scale_x_continuous(position="top",trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "t")+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_color_manual(values=c('darkorange','dodgerblue2'))+
  ggtitle("") +
  ylab("") + 
  xlab("genecopy number/g")+
  theme(legend.position = "none",
        legend.background = element_blank(),axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
pcraz

# 4.SG --------------------------------------------------------------------

pcrsg <- ggplot(data=df01%>% 
                  filter(site=="SG"),
                aes(x=mean_value,
                    y=depths_cm))+
  theme_bw()+
  geom_point(data = df02%>% 
               filter(site=="SG"),aes(x= genecopy, y=depths_cm,color=dna_type),size=3, stat = "identity")+
  geom_path(aes(group=dna_type,linetype= dna_type))+
  scale_x_continuous(position="top",trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "t")+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_color_manual(values=c('darkorange','dodgerblue2'))+
  ggtitle("") +
  ylab("") + 
  xlab("genecopy number/g")+
  theme(legend.position = "none",
        legend.background = element_blank(),axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
pcrsg

# 5. LC -------------------------------------------------------------------

pcrlc <- ggplot(data=df01%>% 
                  filter(site=="LC"),
                aes(x=mean_value,
                    y=depths_cm))+
  theme_bw()+
  geom_point(data = df02%>% 
               filter(site=="LC"),aes(x= genecopy, y=depths_cm,color=dna_type),size=3, stat = "identity")+
  geom_path(aes(group=dna_type,linetype= dna_type))+
  scale_x_continuous(position="top",trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "t")+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_color_manual(values=c('darkorange','dodgerblue2'))+
  ggtitle("") +
  ylab("") + 
  xlab("genecopy number/g")+
  theme(legend.position = "none",
        legend.background = element_blank(),axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
pcrlc

# 6. NA -------------------------------------------------------------------

pcrna <- ggplot(data=df01%>% 
                  filter(site=="NB"),
                aes(x=mean_value,
                    y=depths_cm))+
  theme_bw()+
  geom_point(data = df02%>% 
               filter(site=="NB"),aes(x= genecopy, y=depths_cm,color=dna_type),size=3, stat = "identity")+
  geom_path(aes(group=dna_type,linetype= dna_type))+
  scale_x_continuous(position="top",trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "t")+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_color_manual(values=c('darkorange','dodgerblue2'))+
  ggtitle("") +
  ylab("") + 
  xlab("genecopy number/g")+
  theme(legend.position = "none",
        legend.background = element_blank(),axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
pcrna


#legend
pcrlegend <- ggplot(data=df01%>% 
                  filter(site=="AZ"),
                aes(x=mean_value,
                    y=depths_cm))+
  theme_bw()+
  geom_point(data = df02%>% 
               filter(site=="AZ"),aes(x= genecopy, y=depths_cm,color=dna_type),size=3, stat = "identity")+
  geom_path(aes(group=dna_type,linetype= dna_type))+
  scale_x_continuous(position="top",trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "t")+
  scale_linetype_manual(values=c("twodash", "solid"))+
  scale_color_manual(values=c('darkorange','dodgerblue2'))+
  ggtitle("") +
  ylab("") + 
  xlab("genecopy number/g")+
  theme(legend.position = "bottom",
        legend.background = element_blank(),axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
pcrlegend
