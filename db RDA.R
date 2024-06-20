library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(rdacca.hp)
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names =1,
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                row.names =1,
                header = TRUE)
# test data and hier part env ---------------------------------------------------
# 0. for all -------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv)) %>% 
  select("site","dna_type","depth","pH","Conductivity",
         "moisture","C","N","CN","Feo","Mno","Alo",
         "Sio","NH4","NO3","Po","Pi") %>%
  na.omit() %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]
asv3 <- t(asv2) %>% 
  as.data.frame()
asv_hel_all <- decostand(asv3,
                         method = "hellinger") 
sum(rownames(asv_hel_all) == rownames(env1))#check 
# mod1 <- cca(asv_hel_all ~ depth + pH + Conductivity + moisture+
#               N + C + CN+
#              Feo + Alo + Mno + Sio +
#               NH4 + NO3 + Po + Pi,
#            data = env1)
# vif.cca(mod1)
# mod2 <- cca(asv_hel_all ~ depth + pH + Conductivity + moisture + 
#               CN + Feo + Alo + Mno + Sio + NH4 + NO3 + Po + Pi,
#             data = env1,
#             scale = T)
# vif.cca(mod2)# ok
# summary(mod2)
# mod3 <- dbrda(asv_hel_all ~ depth+ pH + Conductivity + moisture +
#                 N + C + CN+ Feo + Alo + Mno + Sio+
#                 NH4 + NO3 + Po + Pi,
#              data = env1,distance = "bray",scale = T)
# vif.cca(mod3)
mod4 <- dbrda(asv_hel_all ~ depth + pH + Conductivity + moisture + CN+
                Feo + Alo +  Mno + Sio + NH4 + NO3 + Po + Pi,
              data = env1,
              distance = "bray",
              scale = T)
vif.cca(mod4)# ok
RsquareAdj(mod4)

abc<- envfit(mod4 ~ depth + pH + Conductivity + moisture + CN+
         Feo + Alo +  Mno + Sio + NH4 + NO3 + Po + Pi,
       data = env1,
       permu = 999)# all work. YES
abc
p.adjust(abc$vectors$pvals,
         method = 'bonferroni')#校正p值，用来判断
#
adonis2(asv3 ~ dna_type + depth + site, 
        data = env1, 
        permutations = 999)# nice!
#层次分割
tmp <- c("depth","pH","Conductivity","moisture","CN",
         "Feo","Alo","Mno","Sio","NH4","NO3","Po","Pi")
env_hp <- env1[,tmp]

dbRDAhp <- rdacca.hp(vegdist(asv_hel_all),
                     env_hp,
                     method = "dbRDA",
                     type = "adjR2")

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/06_RDA_CCA/i&e_dbrda_hp.pdf'),
    width = 8,
    height = 6)

plot(dbRDAhp, 
     plot.perc = TRUE)

dev.off()
# 1. for iDNA ------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv)) %>% 
  filter(dna_typen == 0) %>% 
  select("site","dna_type","depth","pH","Conductivity","moisture","C","N",
         "CN","Feo","Mno","Alo","Sio","NH4","NO3","Po","Pi") %>% 
  na.omit() %>% 
  data.frame()

filt3 <-  colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt3]
fil4 <- rowSums(asv1) > 0
asv2 = asv1[fil4,]
asv3 <- t(asv2) %>% 
  as.data.frame()
asv_hel_i <- decostand(asv3,
                       method = "hellinger") 
sum(rownames(asv_hel_i) == rownames(env1))#check 

# mod_i1 <- cca(asv3~ pH+Conductivity+moisture+
#                N+C+CN+Feo+Alo+Mno+Sio+Po+Pi+NH4+NO3,
#              data=env1)
# vif.cca(mod_i1)
# 
# mod_i2 <- cca(asv3~ pH+Conductivity+moisture+
#                CN+Feo+Alo+Mno+Sio+Po+Pi+NH4+NO3,
#              data=env1)
# vif.cca(mod_i2)#Al should remove???
# 
# mod_i3 <- dbrda(asv3 ~ depth + pH+Conductivity+moisture+N+C+CN+
#                 Feo+Alo+Mno+Sio+NH4+NO3+Po+Pi,
#               data=env1,distance = "bray")
# vif.cca(mod_i3)

mod_i4 <- dbrda(asv_hel_i~ depth + pH + Conductivity + moisture + CN+
                Feo + Alo +  Mno + Sio + NH4 + NO3 + Po + Pi,
              data = env1,distance = "bray")
vif.cca(mod_i4)#ok
RsquareAdj(mod_i4)

abc_i <- envfit(mod_i4 ~ depth + pH + Conductivity + moisture + CN+
         Feo + Alo +  Mno + Sio + NH4 + NO3 + Po + Pi,
       data = env1,permu=999)# 
abc_i
p.adjust(abc_i$vectors$pvals,
         method = 'bonferroni')
# depth           pH Conductivity     moisture           CN          Feo          Alo          Mno 
# 0.039        0.013        0.117        0.013        0.338        0.013        0.013        0.013 
# Sio          NH4          NO3           Po           Pi 
# 0.013        1.000        1.000        0.013        0.013 

adonis2(asv_hel_i ~ depth + site, 
        data = env1, 
        permutations=999)# 
#层次分割
tmp <- c("depth","pH","Conductivity","moisture","CN",
         "Feo","Alo","Mno","Sio","NH4","NO3","Po","Pi")
env_hp <- env1[,tmp]
dbRDAhp_i <- rdacca.hp(vegdist(asv_hel_i),
                     env_hp,
                     method = "dbRDA",
                     type = "adjR2")
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/06_RDA_CCA/i_dbrda_hp.pdf'),
    width=8,
    height=6)

plot(dbRDAhp_i, 
     plot.perc = TRUE)

dev.off()
#2. for eDNA ------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv)) %>% 
  filter(dna_typen == 1) %>% 
  select("site","dna_type","depth","pH","Conductivity","moisture","C","N",
         "CN","Feo","Mno","Alo","Sio","NH4","NO3","Po","Pi") %>% 
  na.omit() %>% 
  data.frame()

filt3 = colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt3]
fil4 <-  rowSums(asv1)>0
asv2 = asv1[fil4,]
asv3 <- t(asv2) %>% 
  as.data.frame()
asv_hel_e <- decostand(asv3,
                       method = "hellinger") 
sum(rownames(asv_hel_e) == rownames(env1))#check 
# mod_e1 <- cca(asv_hel_e~ pH+Conductivity+moisture+
#                 N+C+CN+Feo+Alo+Mno+Sio+Po+Pi+NH4+NO3,
#               data=env1)
# vif.cca(mod_e1)
# 
# mod_e2 <- cca(asv_hel_e~ pH+Conductivity+moisture+
#                 CN+Feo+Alo+Mno+Sio+Po+Pi+NH4+NO3,
#               data=env1)
# vif.cca(mod_e2)#ok
# 
# mod_e3 <- dbrda(asv_hel_e ~ pH+Conductivity+moisture+N+C+CN+
#                   Feo+Alo+Mno+Sio+NH4+NO3+Po+Pi,
#                 data=env1,distance = "bray")
# vif.cca(mod_e3)
mod_e4 <- dbrda(asv_hel_e ~ depth + pH + Conductivity + moisture + CN+
                  Feo + Alo +  Mno + Sio + NH4 + NO3 + Po + Pi,
                data = env1,
                distance = "bray")
vif.cca(mod_e4)#ok
RsquareAdj(mod_e4)
abc_e <- envfit(mod_e4 ~ depth + pH + Conductivity + moisture + CN+
         Feo + Alo +  Mno + Sio + NH4 + NO3 + Po + Pi,
       data = env1,permu=999)#
abc_e
p.adjust(abc_e$vectors$pvals,
         method = 'bonferroni')
# depth           pH Conductivity     moisture           CN          Feo          Alo          Mno 
# 1.000        0.013        0.052        0.013        0.234        0.013        0.013        0.013 
# Sio          NH4          NO3           Po           Pi 
# 0.013        0.091        1.000        0.013        0.013 
#
adonis2(asv_hel_e ~ depth + site, 
        data = env1, 
        permutations=999)# 
#层次分割
tmp <- c("depth", "pH","Conductivity","moisture","CN",
         "Feo","Alo","Mno","Sio","NH4","NO3","Po","Pi")
env_hp <- env1[,tmp]

dbRDAhp_e <- rdacca.hp(vegdist(asv_hel_e),
                       env_hp,
                       method = "dbRDA",
                       type = "adjR2")
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/06_RDA_CCA/e_dbrda_hp.pdf'),
    width = 8,
    height = 6)

plot(dbRDAhp_e, 
     plot.perc = TRUE)

dev.off()
# plot dbRDA --------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names =1,
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv", 
                row.names =1,
                header = TRUE)
# 0. all DNA-------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv)) %>% 
  filter(depth<200) %>% 
  select("site","dna_type","depth","depth_order",
         "pH","Conductivity","moisture","CN",
         "Feo","Sio","Mno","Alo",
         "NH4","NO3","Po","Pi") %>% 
  na.omit()
filt1 <- colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt1]
asv2 <-asv1 %>% 
  filter(rowSums(.)>0) %>% 
  t(.) %>% 
  as.data.frame()

asv_hel_plot <- decostand(asv2, method = "hellinger") 

dbRDA <- dbrda(asv_hel_plot ~ depth + pH + Conductivity+
                 moisture + CN+
                 Feo + Sio + Mno + Alo+
                 NH4 + NO3 + Po + Pi, 
               data = env1, distance = "bray")

RsquareAdj(dbRDA) 
# $r.squared
# [1] 0.3152576
# 
# $adj.r.squared
# [1] 0.2749786

0.08664/0.3152576#27
0.06608/0.3152576#21
dbRDA1 <- summary(dbRDA,
                  scaling = 1)
dbRDA1
plot(dbRDA)

dbRDA1$sites %>% 
  as.data.frame() -> data.site

identical(rownames(env1),row.names(data.site))

data.site %>% 
  mutate(site = env1$site,
         depth_order = env1$depth_order,
         dna_type = env1$dna_type) -> data.site01

dbRDA1$biplot
dbRDA1$biplot[,1:2] %>% 
  as.data.frame() -> data.biplot

data.biplot %>%
  mutate(dbRDA1=dbRDA1*2,
         dbRDA2=dbRDA2*2) -> data.biplot01

table(data.site01$site)
#plot
data.site01$site <- factor(data.site01$site,
                              levels = c('AZ','SG','LC','NB'))
data.site01$dna_type<-factor(data.site01$dna_type,
                         levels = c("iDNA","eDNA"))
data.site01$depth_order <- factor(data.site01$depth_order,
                                  levels = c("1","2","3","4","5","6","7","8","9","10","11","12"),
                                  labels = c("0-5 cm","5-10 cm","10-20 cm","20-40 cm","40-60 cm","60-80 cm","80-100 cm",
                                             "100-120 cm","120-140 cm","140-160 cm","160-180 cm","180-200 cm"))
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/06_RDA_CCA/i&e_dbrda.pdf'),
    width=8,
    height=6)

ggplot() +
  geom_point(data = data.site01,
             aes(dbRDA1,dbRDA2,
                 color = site,
                 shape = dna_type,
                 alpha = depth_order),
             size = 4)+
  stat_ellipse(data = data.site01,
               aes(dbRDA1,dbRDA2,
                   color=site, 
                   linetype = dna_type),
               geom = "polygon",
               alpha = 0)+
  geom_segment(data = data.biplot01,
               aes(x = 0, 
                   y = 0, 
                   xend = dbRDA1, 
                   yend = dbRDA2), 
               arrow = arrow(angle = 22.5,
                             length = unit(0.35,"cm"),
                             type = "closed"),
               linetype = 1, 
               size = 0.6,
               colour = "black")+
  labs(x = "dbRDA1 (27 %)",
       y = "dbRDA2 (21 %)")+
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Set1")) +
  geom_text_repel(data = data.biplot01,
                  aes(dbRDA1,dbRDA2,
                      label=row.names(data.biplot01)),
                  size = 4,
                  family = "serif",
                  color="black")+
  scale_alpha_ordinal(range = c(0.1,0.75))+
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12))+
  ggtitle("dbrda R2 = 0.3152576 adj R2 = 0.2749786")

dev.off()

# 2. iDNA --------------------------------------------------------------
env_i <- env %>%
  filter(rownames(env) %in% colnames(asv),
         dna_typen == 0) %>% 
  select("site","dna_typen","depth","pH","Conductivity","moisture",
         "CN","Feo","Alo","Sio","Mno","NH4","NO3","Po","Pi",
         "depth_order") %>% 
  na.omit() %>% 
  as.data.frame()

filt1 = colnames(asv) %in% rownames(env_i)
asv_i = asv[,filt1]
filt2 <-  rowSums(asv_i) > 0
asv_iDNA <-  asv_i[filt2,]
t_asv_iDNA <- t(asv_iDNA) %>% 
  as.data.frame()
asv_iDNA_hel <- decostand(t_asv_iDNA,
                          method = "hellinger")

dbRDA_i <- dbrda(asv_iDNA_hel ~ depth + pH + Conductivity + moisture +
                   CN + Feo + Sio + Mno + Alo + NH4 + NO3 + Po + Pi,
                 data = env_i, 
                 distance = "bray")

dbRDA_isum <- summary(dbRDA_i,
                      scaling = 1)
dbRDA_isum
RsquareAdj(dbRDA_i) 
# $r.squared
# [1] 0.3378432
# 
# $adj.r.squared
# [1] 0.2573943
0.08982/0.3378432#27
0.06463/0.3378432#19

dbRDA_isum$sites %>% 
  head()

iDNA_site <- as.data.frame(dbRDA_isum$sites)

identical(rownames(env_i),row.names(iDNA_site))

iDNA_site %>% 
  mutate(site = env_i$site,
         depth_order = env_i$depth_order) -> iDNA_site01

dbRDA_isum$biplot

dbRDA_isum$biplot[,1:2] %>% 
  as.data.frame() -> dbRDA_i_biplot

dbRDA_i_biplot %>% 
  mutate(dbRDA1=dbRDA1*2,
         dbRDA2=dbRDA2*2) -> dbRDA_i_biplot01

table(iDNA_site01$site)

iDNA_site01$site <- factor(iDNA_site01$site,
                           levels = c('AZ','SG','LC','NB'))
iDNA_site01$depth_order <- factor(iDNA_site01$depth_order,
                           levels = c("1","2","3","4","5","6","7","8","9","10","11","12"),
                           labels = c("0-5 cm","5-10 cm","10-20 cm","20-40 cm","40-60 cm","60-80 cm","80-100 cm",
                                      "100-120 cm","120-140 cm","140-160 cm","160-180 cm","180-200 cm"))

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/06_RDA_CCA/i_dbrda V5.pdf'),
    width=8,
    height=6)

ggplot() +
  geom_point(data = iDNA_site01,
             aes(dbRDA1,dbRDA2,
                 color = site,
                 alpha = depth_order),
             size = 3,
             shape = 17)+
  stat_ellipse(data = iDNA_site01,
               aes(dbRDA1,dbRDA2,
                   col = site),
               linetype = 2,
               level = 0.95, 
               show.legend = F)+
  geom_segment(data = dbRDA_i_biplot01,
               aes(x = 0,
                   y = 0,
                   xend = dbRDA1,
                   yend = dbRDA2),
               arrow = arrow(angle=18,
                             length = unit(0.25,"cm"),
                             type = "closed"),
               linetype = 1,
               size = 0.6,
               colour = "black")+
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Set1")) +
  labs(x = "RDA1 (27 %)",
       y = "RDA2 (19 %)")+
  geom_text_repel(data = dbRDA_i_biplot01,
                  aes(dbRDA1,dbRDA2,
                      label = row.names(dbRDA_i_biplot01)),
                  size = 4,
                  family = "serif",
                  color="black")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  ggtitle("iDNA R2 = 0.3378432 and adjR2 = 0.2573943")

dev.off()

# 3. eDNA --------------------------------------------------------------
env_e <- env %>%
  filter(rownames(env) %in% colnames(asv),
         dna_typen==1) %>% 
  filter(depth < 200) %>% 
  select("site","dna_typen","depth","pH","Conductivity","moisture",
         "CN","Feo","Alo","Sio","Mno","NH4","NO3","Po","Pi",
         "depth_order") %>% 
  na.omit()

colnames(env_e)

filt3 = colnames(asv) %in% rownames(env_e)
asv_e = asv[,filt3]
filt4 = rowSums(asv_e)>0
asv_eDNA=asv_e[filt4,]

t_asv_eDNA <- t(asv_eDNA) %>% 
  as.data.frame()
asv_eDNA_hel <- decostand(t_asv_eDNA,method = "hellinger")

dbRDA_e <- dbrda(asv_eDNA_hel ~ depth + pH+Conductivity+moisture+
                   CN+Feo+Sio+Mno+Alo
                 +NH4+NO3+Po+Pi,
                 data=env_e, distance = "bray")
dbRDA_esum <- summary(dbRDA_e,scaling = 1)
dbRDA_esum
RsquareAdj(dbRDA_e) 
# $r.squared
# [1] 0.3910828
# 
# $adj.r.squared
# [1] 0.3119235

0.1036/0.3910828#27
0.08803/0.3910828#23
dbRDA_esum$sites %>% 
  head()

eDNA_site <- as.data.frame(dbRDA_esum$sites)

identical(rownames(env_e),row.names(eDNA_site))

eDNA_site %>% 
  mutate(site=env_e$site,
         depth_order=env_e$depth_order) -> eDNA_site01

dbRDA_esum$biplot

dbRDA_esum$biplot[,1:2] %>% 
  as.data.frame() -> dbRDA_e_biplot

dbRDA_e_biplot %>% 
  mutate(dbRDA1=dbRDA1*2,
         dbRDA2=dbRDA2*2) -> dbRDA_e_biplot01

table(eDNA_site01$site)

eDNA_site01$site <- factor(eDNA_site01$site,
                           levels = c('AZ','SG','LC','NB'))
eDNA_site01$depth_order <- factor(eDNA_site01$depth_order,
                                  levels = c("1","2","3","4","5","6","7","8","9","10","11","12"),
                                  labels = c("0-5 cm","5-10 cm","10-20 cm","20-40 cm","40-60 cm","60-80 cm","80-100 cm",
                                             "100-120 cm","120-140 cm","140-160 cm","160-180 cm","180-200 cm"))

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/06_RDA_CCA/e_dbrda_V5.pdf'),
    width=8,
    height=6)

ggplot() +
  geom_point(data = eDNA_site01,
             aes(dbRDA1,dbRDA2,
                 color = site,
                 alpha = depth_order),
             size = 3,
             shape = 17)+
  stat_ellipse(data = eDNA_site01,
               aes(dbRDA1,dbRDA2,
                   col = site),
               linetype = 2,
               level = 0.95, 
               show.legend = F)+
  geom_segment(data = dbRDA_e_biplot01,
               aes(x = 0, 
                   y = 0, 
                   xend = dbRDA1, 
                   yend = dbRDA2), 
               arrow = arrow(angle=18,
                             length = unit(0.25,"cm"),
                             type = "closed"),
               linetype = 1, 
               size = 0.6,
               colour = "black")+
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Set1")) +
  labs(x = "RDA1 (27 %)",
       y = "RDA2 (23 %)")+
  geom_text_repel(data = dbRDA_e_biplot01,
                  aes(dbRDA1,dbRDA2,
                      label = row.names(dbRDA_e_biplot01)),
                  size = 4,
                  family = "serif",
                  color="black")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  ggtitle("eDNA R2 = 0.3910828 and adj R2 = 0.3119235")

dev.off()
