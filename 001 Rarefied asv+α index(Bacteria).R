library(tidyverse)
library(rtk)
library(vegan)
# load data
asv <- read.delim("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/file_A.txt",
                  header = TRUE,
                  row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_2024.csv",
                header = TRUE)
tax <- read.delim("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt", 
                  header = TRUE)
env200 <- env %>% 
  filter(depth < 201)

asv1 <- asv %>% 
  filter(Kingdom == 'Bacteria')%>%
  select(1:353)

filt1 = colnames(asv1) %in% env200$sample_id
asv2 = asv1[,filt1]  

identical(colnames(asv2),env200$sample_id)#check False
asv3 <- asv2[match(env200$sample_id,colnames(asv2))]
identical(env200$sample_id,colnames(asv3))#check
asv4 <- asv3 %>%
  data.frame() %>% 
  mutate(abc = rowSums(.))%>%
  filter(abc > 0) %>% 
  select(-abc)

col_sums <- colSums(asv4)
sorted_index <- order(col_sums, decreasing = TRUE)
sorted_col_sums <- col_sums[sorted_index]
print(sorted_col_sums)# 9435
# 01 Rarefied ----------------------------------------------------------------
asv_rarefied_9435<-rtk(input = asv4, 
                       depth = 9435, 
                       margin = 2, 
                       tmpdir = NULL, 
                       repeats = 100, 
                       ReturnMatrix = 1) 

asv5 <- asv_rarefied_9435$raremat[[1]]
filt2 <- rowSums(asv5) > 0
asv6 <- asv5[filt2,] %>% 
  as.matrix(.)
write.table(asv6, file = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt", 
            sep = '\t',
            col.names = TRUE,
            row.names = TRUE) 

# get env
env_rare <- env200 %>% 
  filter(sample_id %in% 
           colnames(asv6)) 
write.table(env_rare, file = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/rare_env_200.txt", 
            sep = '\t',col.names = TRUE,row.names = FALSE)
# get taxonomy
tax1 <- tax %>%
  filter(asv_id %in% rownames(asv6) )
write.table(tax1, file = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/tax_rare_200.txt", 
            sep="\t",row.names = F,quote=F) 


# 02 Alpha diversity ------------------------------------------------------
asv7 <- t(asv6)

Shannon <- diversity(asv7, index = 'shannon', base = exp(1))# or base = 2
#Gini-Simpson is normal Simpson,Classic Simpson index <- 1 - Gini_simpson
Gini_simpson  <- diversity(asv7, index = 'simpson')
observed_species <- rowSums(asv7 > 0)

Pielous_evenness_J <- Shannon/log(observed_species)

Chao1  <- estimateR(asv7)[2, ]
ACE  <- estimateR(asv7)[4, ]
goods_coverage <- 1 - rowSums(asv7 == 1) / rowSums(asv7)

pairs(cbind(Shannon,Gini_simpson,observed_species,Pielous_evenness_J,Chao1,ACE,goods_coverage), pch="+", col="red")

alpha <- list(Shannon=Shannon,Simpson=Gini_simpson,observed_species=observed_species,Pielous=Pielous_evenness_J,Chao1=Chao1,ACE=ACE,goods_coverage=goods_coverage)
write.table(alpha, "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/alpha_diversity.txt", 
            sep = "\t",
            row.names =TRUE, 
            col.names =TRUE, 
            quote = FALSE)
rm(list = ls())
dev.off()