# library package ---------------------------------------------------------
library(tidyverse)
library(psych)
library(pheatmap)

# load data ---------------------------------------------------------------

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt",
                header = TRUE)

# 1. i&e together ------------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id%in%colnames(asv))
tax1 <- tax %>% 
  filter(tax$asv_id%in%asv$asv_id)
#get top 20 genus list first###########
top_genus <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = GenusAbun/sum(GenusAbun)) %>%
  ungroup()%>%
  arrange(desc(RelAbun))
top20genus <- top_genus$Genus[2:21]

# tmp=top_genus[2:21,]
# sum(tmp$RelAbun)#0.233265,top20 gen count 23.33%,NA 59.23%

#group by sample id and genus
sample_genus_tab<- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus,sample_id)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = GenusAbun/sum(GenusAbun)) %>%
  ungroup()

plot01 <- sample_genus_tab %>%
  filter(Genus%in%top20genus)

plot01 %>% select(-GenusAbun) %>% 
  pivot_wider(names_from = sample_id,
              values_from = RelAbun)  -> genusname.wider

genusname.wider[is.na(genusname.wider)] <- 0

genusname.wider %>% 
  column_to_rownames("Genus") %>% 
  t() %>% 
  as.data.frame() -> genusname.wider.t

#
identical(env$sample_id,rownames(genusname.wider.t))#check False
env[match(rownames(genusname.wider.t),env$sample_id),] -> new.env#serve
identical(new.env$sample_id,rownames(genusname.wider.t))#check TRUE
##
colnames(new.env)
new.env.subset<-new.env %>% 
  select("dna_typen","depth","pH","Conductivity","moisture","CN","Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi")

# https://github.com/linhesun/bilibiliRlearning/tree/master/20210715
# https://www.bilibili.com/video/BV1vb4y1k7kv?spm_id_from=333.999.0.0

df.corr <- corr.test(new.env.subset, genusname.wider.t, use = 'pairwise', 
                 method = 'spearman', adjust = 'holm',
                 alpha = 0.05)
#p值是假设检验中假设零假设为真时
#观测到至少与实际观测样本Realization (probability)
#相同极端的样本的概率
relation <- df.corr$r 
p_value <- df.corr$p
pheatmap(df.corr$r)

df.corr$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>%  
  mutate(new_col=case_when(
    #value == NA ~ '',
    value < 0.001 ~ '***',
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    value >=0.05 & value < 0.1 ~ ".",
    value >=0.1 ~ '',
    TRUE ~ ''
  )) %>% 
  select(rowname,variable,new_col) %>% 
  pivot_wider(names_from = "variable",
              values_from = "new_col") %>% 
  column_to_rownames('rowname') -> p.signi

# pheatmap(df.corr$r,
#          display_numbers = p.signi)

#if you want change axis,run this one

t_relation <- t(df.corr$r)
t_p_value <- t(df.corr$p)
pheatmap(t_relation)

t_p_value %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>%  
  mutate(new_col=case_when(
    #value == NA ~ '',
    value < 0.001 ~ '***',
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    #value >=0.05 & value < 0.1 ~ ".",
    value >=0.1 ~ '',
    TRUE ~ ''
  )) %>% 
  select(rowname,variable,new_col) %>% 
  pivot_wider(names_from = "variable",
              values_from = "new_col") %>% 
  column_to_rownames('rowname') -> p.signi

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/07_cor heatmap with sig/top20gen_env.pdf'),
    width = 8,height = 6)

pheatmap(t_relation,
         display_numbers = p.signi
         # cluster_rows = FALSE,
         # cluster_cols = FALSE
         )
dev.off()

# 2. iDNA load data ---------------------------------------------------------------

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt",
                  header = TRUE)

# iDNA data clean ------------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id%in%colnames(asv)&dna_typen==0)
asv <- asv[,c(1,128:258)]#iDNA

tax1 <- tax %>% 
  filter(tax$asv_id%in%asv$asv_id)
#get top 20 genus list first###########
top_genus <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = GenusAbun/sum(GenusAbun)) %>%
  ungroup()%>%
  arrange(desc(RelAbun))
top20genus <- top_genus$Genus[2:21]

# tmp=top_genus[2:21,]
# sum(tmp$RelAbun)#0.233265,top20 gen count 23.33%,NA 59.23%

#group by sample id and genus
sample_genus_tab<- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus,sample_id)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = GenusAbun/sum(GenusAbun)) %>%
  ungroup()

plot01 <- sample_genus_tab %>%
  filter(Genus%in%top20genus)

plot01 %>% select(-GenusAbun) %>% 
  pivot_wider(names_from = sample_id,
              values_from = RelAbun)  -> genusname.wider

genusname.wider[is.na(genusname.wider)] <- 0

genusname.wider %>% 
  column_to_rownames("Genus") %>% 
  t() %>% 
  as.data.frame() -> genusname.wider.t

#
identical(env$sample_id,rownames(genusname.wider.t))#check False
env[match(rownames(genusname.wider.t),env$sample_id),] -> new.env#serve
identical(new.env$sample_id,rownames(genusname.wider.t))#check TRUE
##
colnames(new.env)
new.env.subset<-new.env %>% 
  select("depth","pH","Conductivity","moisture","CN",
         "Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi")

df.corr <- corr.test(new.env.subset, genusname.wider.t, use = 'pairwise', 
                     method = 'spearman', adjust = 'holm',
                     alpha = 0.05)
#p值是假设检验中假设零假设为真时
#观测到至少与实际观测样本Realization (probability)
#相同极端的样本的概率
relation <- df.corr$r 
p_value <- df.corr$p
pheatmap(df.corr$r)

df.corr$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>%  
  mutate(new_col=case_when(
    #value == NA ~ '',
    value < 0.001 ~ '***',
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    value >=0.05 & value < 0.1 ~ ".",
    value >=0.1 ~ '',
    TRUE ~ ''
  )) %>% 
  select(rowname,variable,new_col) %>% 
  pivot_wider(names_from = "variable",
              values_from = "new_col") %>% 
  column_to_rownames('rowname') -> p.signi

# pheatmap(df.corr$r,
#          display_numbers = p.signi)

#if you want change axis,run this one

t_relation <- t(df.corr$r)
t_p_value <- t(df.corr$p)
pheatmap(t_relation)

t_p_value %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>%  
  mutate(new_col=case_when(
    #value == NA ~ '',
    value < 0.001 ~ '***',
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    #value >=0.05 & value < 0.1 ~ ".",
    value >=0.1 ~ '',
    TRUE ~ ''
  )) %>% 
  select(rowname,variable,new_col) %>% 
  pivot_wider(names_from = "variable",
              values_from = "new_col") %>% 
  column_to_rownames('rowname') -> p.signi

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/07_cor heatmap with sig/top20gen_envidan.pdf'),
    width = 8,height = 6)

pheatmap(t_relation,
         display_numbers = p.signi
         # cluster_rows = FALSE,
         # cluster_cols = FALSE
)
dev.off()

# load data for eDNA  ---------------------------------------------------------------

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt",
                  header = TRUE)

# eDNA data clean ------------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id%in%colnames(asv)&dna_typen==1)
asv <- asv[,c(1:127)]#eDNA

tax1 <- tax %>% 
  filter(tax$asv_id%in%asv$asv_id)
#get top 20 genus list first###########
top_genus <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = GenusAbun/sum(GenusAbun)) %>%
  ungroup()%>%
  arrange(desc(RelAbun))
top20genus <- top_genus$Genus[2:21]

# tmp=top_genus[2:21,]
# sum(tmp$RelAbun)#0.233265,top20 gen count 23.33%,NA 59.23%

#group by sample id and genus
sample_genus_tab<- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Genus,sample_id)%>%
  summarise(GenusAbun = sum(Abun)) %>%
  mutate(RelAbun = GenusAbun/sum(GenusAbun)) %>%
  ungroup()

plot01 <- sample_genus_tab %>%
  filter(Genus%in%top20genus)

plot01 %>% select(-GenusAbun) %>% 
  pivot_wider(names_from = sample_id,
              values_from = RelAbun)  -> genusname.wider

genusname.wider[is.na(genusname.wider)] <- 0

genusname.wider %>% 
  column_to_rownames("Genus") %>% 
  t() %>% 
  as.data.frame() -> genusname.wider.t

#
identical(env$sample_id,rownames(genusname.wider.t))#check False
env[match(rownames(genusname.wider.t),env$sample_id),] -> new.env#serve
identical(new.env$sample_id,rownames(genusname.wider.t))#check TRUE
##
colnames(new.env)
new.env.subset<-new.env %>% 
  select("depth","pH","Conductivity","moisture","CN",
         "Fed","FeoFed","Fep","Sio","Mnd","NH4","NO3","Pt","Pi")

# https://github.com/linhesun/bilibiliRlearning/tree/master/20210715
# https://www.bilibili.com/video/BV1vb4y1k7kv?spm_id_from=333.999.0.0

df.corr <- corr.test(new.env.subset, genusname.wider.t, use = 'pairwise', 
                     method = 'spearman', adjust = 'holm',
                     alpha = 0.05)
#p值是假设检验中假设零假设为真时
#观测到至少与实际观测样本Realization (probability)
#相同极端的样本的概率
relation <- df.corr$r 
p_value <- df.corr$p
pheatmap(df.corr$r)

df.corr$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>%  
  mutate(new_col=case_when(
    #value == NA ~ '',
    value < 0.001 ~ '***',
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    value >=0.05 & value < 0.1 ~ ".",
    value >=0.1 ~ '',
    TRUE ~ ''
  )) %>% 
  select(rowname,variable,new_col) %>% 
  pivot_wider(names_from = "variable",
              values_from = "new_col") %>% 
  column_to_rownames('rowname') -> p.signi

# pheatmap(df.corr$r,
#          display_numbers = p.signi)

#if you want change axis,run this one

t_relation <- t(df.corr$r)
t_p_value <- t(df.corr$p)
pheatmap(t_relation)

t_p_value %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  reshape2::melt() %>%  
  mutate(new_col=case_when(
    #value == NA ~ '',
    value < 0.001 ~ '***',
    value >= 0.001 & value < 0.01 ~ "**",
    value >= 0.01 & value < 0.05 ~ "*",
    #value >=0.05 & value < 0.1 ~ ".",
    value >=0.1 ~ '',
    TRUE ~ ''
  )) %>% 
  select(rowname,variable,new_col) %>% 
  pivot_wider(names_from = "variable",
              values_from = "new_col") %>% 
  column_to_rownames('rowname') -> p.signi

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/07_cor heatmap with sig/top20gen_envedan.pdf'),
    width = 8,height = 6)

pheatmap(t_relation,
         display_numbers = p.signi
         # cluster_rows = FALSE,
         # cluster_cols = FALSE
)
dev.off()

