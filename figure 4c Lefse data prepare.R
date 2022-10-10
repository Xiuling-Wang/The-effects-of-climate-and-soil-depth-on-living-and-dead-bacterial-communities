library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")
library(phyloseq)
library(magrittr) 

#read data
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1, header = TRUE, sep='\t',check.names = FALSE)#经过抽平的数据
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header = TRUE, row.names = 1,sep = ",", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "")#所有样本的分组数据
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt", 
                  header = TRUE, row.names = 1,sep='\t',check.names = FALSE)##所有asv的taxnomy数据
tmp0 <- rownames(env) %in% colnames(asv)
env <- env[tmp0,]# make sure sample in env both belong to rarefied asv table
# 1. iDNA -----------------------------------------------------------------
env <- env %>%
  filter(dna_type=="iDNA",
         depth<20)
env$site <- factor(env$site,levels = c("AZ","SG","LC","NB"),
                   labels = c("4AZ","3SG","2LC","1NB"))
tmp1 <- colnames(asv) %in% rownames(env)
asv <- asv[,tmp1]

asv <- asv %>%
  filter(rowSums(asv)>0)

asv_long <- t(asv)
ra <- asv_long/rowSums(asv_long)
asv <- t(ra)
tax <- tax %>% 
  select(!c("Kingdom","Phylum"))%>%
  filter(!Class=="c__NA")

tmp2 <- rownames(tax) %in% rownames(asv)
tax <- tax[tmp2,]

tmp3 <- rownames(asv) %in% rownames(tax)
asv <- asv[tmp3,]

asv_tax <- merge(asv,tax,by = "row.names",all=F)#合并
head(asv_tax)

tax_levels <- names(tax)
tax_lst <- list()

for (idx in seq(4)){
HA_tax <- asv_tax %>%
  filter(!str_detect(!!sym(tax_levels[idx]), 'NA')) %>%
  group_by(!!sym(tax_levels[idx])) %>%
  summarise(across(where(is.numeric),sum))
HA_tax <- tax %>%
  filter(!str_detect(!!sym(tax_levels[idx]), 'NA')) %>%
  `[`(1:idx) %>%
  distinct() %>%
  unite(taxonomy,sep = '|', remove = FALSE) %>%
  select(!!sym(tax_levels[idx]), taxonomy) %>%
  left_join(HA_tax) %>%
  select(-!!sym(tax_levels[idx]))
tax_lst[[idx]] <- HA_tax 
}

# tax_lst[[7]] <- asv %>% 
#   as.data.frame() %>% 
#   mutate(taxonomy=rownames(.))

all_asv <- do.call(rbind, tax_lst)
rownames(all_asv)=all_asv[,1]
all_asv=all_asv[,-1]

group <- env %>% 
  dplyr::select(site) %>% #,layer,dna_type site_dna
  t() %>% 
  data.frame()
  
lefse_data <- rbind(group,all_asv)
head(lefse_data[,1:6])

write.table(lefse_data,"/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/lefse data/idna0-20.txt",
            sep="\t",row.names = T,col.names = F, quote=F)

# 2.  eDNA ----------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1, header = TRUE, sep='\t',check.names = FALSE)#经过抽平的数据
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header = TRUE, row.names = 1,sep = ",", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "")#所有样本的分组数据
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt", 
                  header = TRUE, row.names = 1,sep='\t',check.names = FALSE)##所有asv的taxnomy数据
tmp0 <- rownames(env) %in% colnames(asv)
env <- env[tmp0,]# make sure sample in env both belong to rarefied asv table

env <- env %>%
  filter(dna_type=="eDNA",
         depth<20)
env$site <- factor(env$site,levels = c("AZ","SG","LC","NB"),
                   labels = c("4AZ","3SG","2LC","1NB"))
tmp1 <- colnames(asv) %in% rownames(env)
asv <- asv[,tmp1]

asv <- asv %>%
  filter(rowSums(asv)>0)

asv_long <- t(asv)
ra <- asv_long/rowSums(asv_long)
asv <- t(ra)
tax <- tax %>% 
  select(!c("Kingdom","Phylum"))%>%
  filter(!Class=="c__NA")

tmp2 <- rownames(tax) %in% rownames(asv)
tax <- tax[tmp2,]

tmp3 <- rownames(asv) %in% rownames(tax)
asv <- asv[tmp3,]

asv_tax <- merge(asv,tax,by = "row.names",all=F)#合并
head(asv_tax)

tax_levels <- names(tax)
tax_lst <- list()

for (idx in seq(4)){
  HA_tax <- asv_tax %>%
    filter(!str_detect(!!sym(tax_levels[idx]), 'NA')) %>%
    group_by(!!sym(tax_levels[idx])) %>%
    summarise(across(where(is.numeric),sum))
  HA_tax <- tax %>%
    filter(!str_detect(!!sym(tax_levels[idx]), 'NA')) %>%
    `[`(1:idx) %>%
    distinct() %>%
    unite(taxonomy,sep = '|', remove = FALSE) %>%
    select(!!sym(tax_levels[idx]), taxonomy) %>%
    left_join(HA_tax) %>%
    select(-!!sym(tax_levels[idx]))
  tax_lst[[idx]] <- HA_tax 
}

# tax_lst[[7]] <- asv %>% 
#   as.data.frame() %>% 
#   mutate(taxonomy=rownames(.))

all_asv <- do.call(rbind, tax_lst)
rownames(all_asv)=all_asv[,1]
all_asv=all_asv[,-1]

group <- env %>% 
  dplyr::select(site) %>% #,layer,dna_type site_dna
  t() %>% 
  data.frame()

lefse_data <- rbind(group,all_asv)
head(lefse_data[,1:6])

write.table(lefse_data,"/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/lefse data/edna0-20cm.txt",
            sep="\t",row.names = T,col.names = F, quote=F)

# new method --------------------------------------------------------------

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1, header = TRUE, sep='\t',check.names = FALSE)#经过抽平的数据
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header = TRUE, row.names = 1,sep = ",", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "")#所有样本的分组数据
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt", 
                  header = TRUE, row.names = 1,sep='\t',check.names = FALSE)##所有asv的taxnomy数据
tmp0 <- rownames(env) %in% colnames(asv)
env <- env[tmp0,]# make sure sample in env both belong to rarefied asv table

env <- env %>%
  filter(dna_type=="eDNA",
         layer_n==3)
env$site <- factor(env$site,levels = c("AZ","SG","LC","NB"),
                   labels = c("4AZ","3SG","2LC","1NB"))
tmp1 <- colnames(asv) %in% rownames(env)
asv <- asv[,tmp1]
asv <- asv %>%
  filter(rowSums(asv)>0)

asv_long <- t(asv)
ra <- asv_long/rowSums(asv_long)
asv <- t(ra)

tax <- tax %>%
  filter(!Phylum=="p__NA")

tmp2 <- rownames(tax) %in% rownames(asv)
tax <- tax[tmp2,]

tmp3 <- rownames(asv) %in% rownames(tax)
asv <- asv[tmp3,]

asv_tax <- merge(asv,tax,by = "row.names",all=F)#合并
head(asv_tax)

tax_levels <- names(tax)
tax_lst <- list()

for (idx in seq(6)){
  HA_tax <- asv_tax %>%
    filter(!str_detect(!!sym(tax_levels[idx]), 'NA')) %>%
    group_by(!!sym(tax_levels[idx])) %>%
    summarise(across(where(is.numeric),sum))
  HA_tax <- tax %>%
    filter(!str_detect(!!sym(tax_levels[idx]), 'NA')) %>%
    `[`(1:idx) %>%
    distinct() %>%
    unite(taxonomy,sep = '|', remove = FALSE) %>%
    select(!!sym(tax_levels[idx]), taxonomy) %>%
    left_join(HA_tax) %>%
    select(-!!sym(tax_levels[idx]))
  tax_lst[[idx]] <- HA_tax 
}

# tax_lst[[7]] <- asv %>% 
#   as.data.frame() %>% 
#   mutate(taxonomy=rownames(.))

all_asv <- do.call(rbind, tax_lst)
rownames(all_asv)=all_asv[,1]
all_asv=all_asv[,-1]

group <- env %>% 
  dplyr::select(site) %>% #,layer,dna_type site_dna
  t() %>% 
  data.frame()

lefse_data <- rbind(group,all_asv)
head(lefse_data[,1:6])
ps_demo<-phyloseq(otu_table(comm,TRUE), tax_table(taxonomy), sample_data(treat))

