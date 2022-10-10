# library(remotes)
# remotes::install_github("GuangchuangYu/gglayer")
# library package ---------------------------------------------------------
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(cowplot)
# read data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                row.names =1,header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                row.names =1,header = TRUE)
colnames(env)


# 1. all data -------------------------------------------------------------
# clean data --------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env)%in%colnames(asv)) %>% 
  filter(depth<200) %>% #look at surface
  #filter(layer_n==1) %>% 
  select("site","dna_type","dna_typen","depth","pH","Conductivity","moisture",
         "N15","CN","Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>% 
  na.omit()

# depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi+dna_typen
# summary(env1)
# colnames(env1)

filt1 <- colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt1]
asv2 <-asv1 %>% 
  filter(rowSums(.)>0) %>% 
  t(.) %>% 
  as.data.frame()

dbRDA <- dbrda(asv2 ~ dna_typen+depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi,
                  data=env1, distance = "bray")
RsquareAdj(dbRDA) 
dbRDA1 <- summary(dbRDA,scaling = 1)
dbRDA1
0.07876/0.3209777
0.0550/0.3209777
dbRDA1$species %>% 
  head()
dbRDA1$sites %>% 
  head()
species <- as.data.frame(dbRDA1$species)

dbRDA1$sites %>% 
  as.data.frame() -> data.site

identical(rownames(env1),row.names(data.site))

data.site %>% 
  mutate(site=env1$site,
         dna_type=env1$dna_type) -> data.site01

dbRDA1$biplot

dbRDA1$biplot[,1:2] %>% 
  as.data.frame() -> data.biplot

data.biplot %>% 
  mutate(dbRDA1=dbRDA1*2,
         dbRDA2=dbRDA2*2) -> data.biplot01

table(data.site01$site)
# plot --------------------------------------------------------------------
data.site01$site <- factor(data.site01$site,
                              levels = c('AZ','SG','LC','NB'))
data.site01$dna_type<-factor(data.site01$dna_type,
                         levels = c("iDNA","eDNA"))

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/06_RDA_CCA/i&e_dbrda.pdf'),
              width=8,height=6)
ggplot() +
  geom_point(data = data.site01,
             aes(dbRDA1,dbRDA2,color=site,shape = dna_type),
             size=4,alpha=0.2)+
  stat_ellipse(data=data.site01,
               aes(dbRDA1,dbRDA2,color=site),
               geom = "polygon",
               alpha=0)+
  geom_segment(data = data.biplot01,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),
               linetype=1, 
               size=0.6,
               colour = "red")+
  labs(x="dbRDA1 (24.54%)",y="dbRDA2 (17.14%)")+
  geom_text_repel(data = data.biplot01,
                  aes(dbRDA1,dbRDA2,label=row.names(data.biplot01)),
                  color="blue")+
  theme_bw()+
  theme(axis.text.x= element_text(size=12),
      axis.text.y= element_text(size=12))+
  ggtitle("dbrda_0.3209777&0.2711412")

dev.off()

# 2. iDNA dbRDA --------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env)%in%colnames(asv),dna_typen==0) %>% 
  filter(depth<200) %>% 
  select("site","dna_typen","depth","pH","Conductivity","moisture",
         "N15","Feo","CN","FeoFed","Fed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>% 
  na.omit()

colnames(env1)

filt1=colnames(asv)%in%rownames(env1)

asv1=asv[,filt1]
filt2=rowSums(asv1)>0
asv1.5=asv1[filt2,]

asv2 <- t(asv1.5) %>% 
  as.data.frame()
colnames(env)
dbRDA <- dbrda(asv2 ~ depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi,
                  data=env1, distance = "bray")
dbRDA1 <- summary(dbRDA,scaling = 1)
dbRDA1
RsquareAdj(dbRDA) 
0.08991/0.3451308
0.05567/0.3451308
dbRDA1$species %>% 
  head()
dbRDA1$sites %>% 
  head()
species <- as.data.frame(dbRDA1$species)
dbRDA1$sites %>% 
  as.data.frame() -> data.site
identical(rownames(env1),row.names(data.site))
data.site %>% 
  mutate(site=env1$site) -> data.site01

dbRDA1$biplot

dbRDA1$biplot[,1:2] %>% 
  as.data.frame() -> data.biplot
data.biplot %>% 
  mutate(dbRDA1=dbRDA1*2,
         dbRDA2=dbRDA2*2) -> data.biplot01

table(data.site01$site)

data.site01$site <- factor(data.site01$site,
                           levels = c('AZ','SG','LC','NB'))

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/06_RDA_CCA/i_dbrda.pdf'),
    width=8,height=6)

ggplot() +
  geom_point(data = data.site01,
             aes(dbRDA1,dbRDA2,color=site),
             size=4,alpha=0.2)+
  stat_ellipse(data=data.site01,
               aes(dbRDA1,dbRDA2,color=site),
               geom = "polygon",
               alpha=0)+
  geom_segment(data = data.biplot01,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),
               linetype=1, 
               size=0.6,
               colour = "red")+
  labs(x="dbRDA1 (26.05%)",y="dbRDA2 (16.13%)")+
  geom_text_repel(data = data.biplot01,
                  aes(dbRDA1,dbRDA2,label=row.names(data.biplot01)),
                  color="blue")+
  theme_bw()+
  theme(axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12))+
  ggtitle("bray_dbrda_iDNA_r0.3451308adj.r0.251578")

dev.off()

# 3. eDNA dbRDA --------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env)%in%colnames(asv),dna_typen==1) %>%
  select("site","dna_typen","depth","pH","Conductivity","moisture",
         "N15","Feo","CN","FeoFed","Fed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi") %>% 
  filter(depth<200) %>% 
  na.omit()
colnames(env1)

filt1=colnames(asv)%in%rownames(env1)

asv1=asv[,filt1]
filt2=rowSums(asv1)>0
asv1.5=asv1[filt2,]

asv2 <- t(asv1.5) %>% 
  as.data.frame()
colnames(env)

dbRDA <- dbrda(asv2 ~ depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi,
                  data=env1, distance = "bray")
RsquareAdj(dbRDA) 
dbRDA1 <- summary(dbRDA,scaling = 1)
dbRDA1
0.09469/0.3896412
0.07766/0.3896412
dbRDA1$species %>% 
  head()
dbRDA1$sites %>% 
  head()
species <- as.data.frame(dbRDA1$species)
dbRDA1$sites %>% 
  as.data.frame() -> data.site
identical(rownames(env1),row.names(data.site))
data.site %>% 
  mutate(site=env1$site) -> data.site01

dbRDA1$biplot

dbRDA1$biplot[,1:2] %>% 
  as.data.frame() -> data.biplot
data.biplot %>% 
  mutate(dbRDA1=dbRDA1*2,
         dbRDA2=dbRDA2*2) -> data.biplot01

table(data.site01$site)

data.site01$site <- factor(data.site01$site,
                           levels = c('AZ','SG','LC','NB'))
pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/06_RDA_CCA/e_dbrda.pdf'),
    width=8,height=6)

ggplot() +
  geom_point(data = data.site01,
             aes(dbRDA1,dbRDA2,color=site),
             size=4,alpha=0.2)+
  stat_ellipse(data=data.site01,
               aes(dbRDA1,dbRDA2,color=site),
               geom = "polygon",
               alpha=0)+
  geom_segment(data = data.biplot01,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),
               linetype=1, 
               size=0.6,
               colour = "red")+
  labs(x="dbRDA1 (24.30%)",y="dbRDA2 (19.93%)")+
  geom_text_repel(data = data.biplot01,
                  aes(dbRDA1,dbRDA2,label=row.names(data.biplot01)),
                  color="blue")+
  theme_bw()+
  theme(axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12))+
  ggtitle("bray_dbrda_eDNA_0.3896412&0.2962189")

dev.off()
