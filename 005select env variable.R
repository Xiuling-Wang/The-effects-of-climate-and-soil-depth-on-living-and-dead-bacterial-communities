library(tidyverse)
library(vegan)
#如果一起做模型，效果不好，分开做，解释度都会提升
# 1. all data --------------------------------------------
# read table --------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                row.names =1,header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                row.names =1,header = TRUE)
# filter data -------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env)%in%colnames(asv)) %>% 
  select("dna_typen","depth","pH","Conductivity","moisture",
         "N","N15","C","CN","Feo","Fed","FeoFed","Fep",
         "Alo","Alp","Sio","NH4","NO3","Pt","Pi","Mnd","Mno","Mnp")%>% 
  na.omit() %>% 
  data.frame()

filt1 <- colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]
asv3 <- t(asv2) %>% 
  as.data.frame()

sum(rownames(asv3) == rownames(env1))#check 
mod <- cca(asv3~dna_typen+depth+pH+Conductivity+moisture+N+N15+C+CN+
             Feo+Fed+FeoFed+Fep+Alo+Alp+Mnd+Mno+Mnp+Sio+NH4+NO3+Pt+Pi,
           data=env1)

# check every single env --------------------------------------------------
mod <- cca(asv3~dna_typen+depth+pH+Conductivity+moisture+N+N15+C+CN+
               Feo+Fed+FeoFed+Fep+Alo+Alp+Mnd+Mno+Mnp+Sio+NH4+NO3+Pt+Pi,
             data=env1)


envfit(mod~dna_typen+depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi,
       data=env1,permu=999)

#just keep with sig


adonis2(asv3 ~ dna_typen, data=env1, permutations=99)
adonis2(asv3 ~ depth, data=env1, permutations=99)
adonis2(asv3 ~ pH, data=env1, permutations=99)
adonis2(asv3 ~ Conductivity, data=env1, permutations=99)
adonis2(asv3 ~ moisture, data=env1, permutations=99)
adonis2(asv3 ~ C, data=env1, permutations=99)
adonis2(asv3 ~ N, data=env1, permutations=99)
adonis2(asv3 ~ N15, data=env1, permutations=99)
adonis2(asv3 ~ CN, data=env1, permutations=99)
adonis2(asv3 ~ Feo, data=env1, permutations=99)
adonis2(asv3 ~ Fed, data=env1, permutations=99)
adonis2(asv3 ~ Alo, data=env1, permutations=99)
adonis2(asv3 ~ Alp, data=env1, permutations=99)
adonis2(asv3 ~ Mnd, data=env1, permutations=99)
adonis2(asv3 ~ Mno, data=env1, permutations=99)
adonis2(asv3 ~ Mnp, data=env1, permutations=99)
adonis2(asv3 ~ Sio, data=env1, permutations=99)
adonis2(asv3 ~ NH4, data=env1, permutations=99)
adonis2(asv3 ~ NO3, data=env1, permutations=99)
adonis2(asv3 ~ Pt, data=env1, permutations=99)
adonis2(asv3 ~ Pi, data=env1, permutations=99)

# vif test ----------------------------------------------------------------
mod <- dbrda(asv3~dna_typen+depth+pH+Conductivity+moisture+N+N15+C+CN+
             Feo+Fed+FeoFed+Fep+Alo+Alp+Mnd+Mno+Mnp+Sio+NH4+NO3+Pt+Pi,
             data=env1,distance = "bray")
summary(mod)
RsquareAdj(mod)

vif.cca(mod)

# dna_typen        depth           pH Conductivity     moisture            N          N15            C           CN 
# 1.018318     2.269696     4.137979     3.637432     3.758287    15.100905     2.399581    17.335071     2.292490 
# Feo          Fed       FeoFed          Fep          Alo          Alp          Mnd          Mno          Mnp 
# 11.135368     3.646109     6.587623     7.859695    26.709400    17.003423     1.925777    27.645277    20.961780 
# Sio          NH4          NO3           Pt           Pi 
# 3.109083     1.610932     1.806724     1.354687     3.329604 

#delete vif>10 
# got result --------------------------------------------------------------
dna_typen+depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi
#test result(remove coliner)
mod <- dbrda(asv3~dna_typen+depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi,
               data=env1,distance = "bray")
summary(mod)
RsquareAdj(mod)
vif.cca(mod)

# 2. test db_rda mod for iDNA ------------------------------------------------
# import data -------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names =1,header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                row.names =1,header = TRUE)
# filter data -------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env)%in%colnames(asv)) %>% 
  filter(dna_typen==0) %>% 
  select("depth","pH","Conductivity","moisture",
         "N","N15","C","CN","Feo","Fed","FeoFed","Fep",
         "Alo","Alp","Sio","NH4","NO3","Pt","Pi","Mnd","Mno","Mnp")%>% 
  na.omit() %>% 
  data.frame()

filt3 = colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt3]
fil4 = rowSums(asv1>0)
asv2 = asv1[fil4,]
asv3 <- t(asv2) %>% 
  as.data.frame()

sum(rownames(asv3) == rownames(env1))#check 
# check every single env --------------------------------------------------
# just keep with sig


adonis2(asv3 ~ depth, data=env1, permutations=99)
adonis2(asv3 ~ pH, data=env1, permutations=99)
adonis2(asv3 ~ Conductivity, data=env1, permutations=99)
adonis2(asv3 ~ moisture, data=env1, permutations=99)
adonis2(asv3 ~ C, data=env1, permutations=99)
adonis2(asv3 ~ N, data=env1, permutations=99)
adonis2(asv3 ~ N15, data=env1, permutations=99)
adonis2(asv3 ~ CN, data=env1, permutations=99)#no sig
adonis2(asv3 ~ Feo, data=env1, permutations=99)
adonis2(asv3 ~ Fed, data=env1, permutations=99)#no sig
adonis2(asv3 ~ Alo, data=env1, permutations=99)
adonis2(asv3 ~ Alp, data=env1, permutations=99)
adonis2(asv3 ~ Mnd, data=env1, permutations=99)
adonis2(asv3 ~ Mno, data=env1, permutations=99)
adonis2(asv3 ~ Mnp, data=env1, permutations=99)
adonis2(asv3 ~ Sio, data=env1, permutations=99)
adonis2(asv3 ~ NH4, data=env1, permutations=99)#no sig
adonis2(asv3 ~ NO3, data=env1, permutations=99)#no sig
adonis2(asv3 ~ Pt, data=env1, permutations=99)
adonis2(asv3 ~ Pi, data=env1, permutations=99)

mod_i <- cca(asv3~depth+pH+Conductivity+moisture+N+N15+C+
               Feo+FeoFed+Fep+Alo+Alp+Mnd+Mno+Mnp+Sio+Pt+Pi,
               data=env1)
vif.cca(mod_i)
# depth           pH Conductivity     moisture            N          N15            C          Feo 
# 1.794344     4.009088     4.133513     3.818880    11.546892     1.805880    11.351002     5.819695 
# FeoFed          Fep          Alo          Alp          Mnd          Mno          Mnp          Sio 
# 3.475032     6.453162    22.953499    11.432092     1.833886    27.380707    18.822111     2.817798 
# Pt           Pi 
# 1.206281     2.557517 
#result
depth+pH+Conductivity+moisture+N15+Feo+FeoFed+Fep+Alp+Mnd+Sio+Pt+Pi
#test result(remove coliner)
mod_i <- dbrda(asv3~depth+pH+Conductivity+moisture+N15+Feo+FeoFed+Fep+Alp+Mnd+Sio+Pt+Pi,
             data=env1,distance = "bray")
#分析中整体模型是否显著以及每个变量的重要性
permutest(mod_i,permu=999)
envfit(mod_i~depth+pH+Conductivity+moisture+N15+Feo+FeoFed+Fep+Mnd+Sio+Pt+Pi,
       data=env1,permu=999)

summary(mod_i)
RsquareAdj(mod_i)
vif.cca(mod_i)


mod_i <- cca(asv3~depth+pH+Conductivity+moisture+N15+Feo+FeoFed+Fep+Mnd+Sio+Pt+Pi,
               data=env1)



#3.  test db_rda mod for eDNA ------------------------------------------------
# import data -------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names =1,header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                row.names =1,header = TRUE)
# filter data -------------------------------------------------------------
env1 <- env %>%
  filter(rownames(env)%in%colnames(asv)) %>% 
  filter(dna_typen==1) %>% 
  select("depth","pH","Conductivity","moisture",
         "N","N15","C","CN","Feo","Fed","FeoFed","Fep",
         "Alo","Alp","Sio","NH4","NO3","Pt","Pi","Mnd","Mno","Mnp")%>% 
  na.omit() %>% 
  data.frame()

filt3 = colnames(asv)%in%rownames(env1)
asv1 <- asv[,filt3]
fil4 = rowSums(asv1>0)
asv2 = asv1[fil4,]
asv3 <- t(asv2) %>% 
  as.data.frame()

sum(rownames(asv3) == rownames(env1))#check 
# check every single env --------------------------------------------------
# just keep with sig
adonis2(asv3 ~ depth, data=env1, permutations=99)
adonis2(asv3 ~ pH, data=env1, permutations=99)
adonis2(asv3 ~ Conductivity, data=env1, permutations=99)
adonis2(asv3 ~ moisture, data=env1, permutations=99)
adonis2(asv3 ~ C, data=env1, permutations=99)
adonis2(asv3 ~ N, data=env1, permutations=99)
adonis2(asv3 ~ N15, data=env1, permutations=99)
adonis2(asv3 ~ CN, data=env1, permutations=99)
adonis2(asv3 ~ Feo, data=env1, permutations=99)
adonis2(asv3 ~ Fed, data=env1, permutations=99)
adonis2(asv3 ~ Alo, data=env1, permutations=99)
adonis2(asv3 ~ Alp, data=env1, permutations=99)
adonis2(asv3 ~ Mnd, data=env1, permutations=99)
adonis2(asv3 ~ Mno, data=env1, permutations=99)
adonis2(asv3 ~ Mnp, data=env1, permutations=99)
adonis2(asv3 ~ Sio, data=env1, permutations=99)
adonis2(asv3 ~ NH4, data=env1, permutations=99)
adonis2(asv3 ~ NO3, data=env1, permutations=99)
adonis2(asv3 ~ Pt, data=env1, permutations=99)
adonis2(asv3 ~ Pi, data=env1, permutations=99)#all hav sig

mod_e <- dbrda(asv3~depth+pH+Conductivity+moisture+N+N15+C+CN+
                 Feo+Fed+FeoFed+Fep+Alo+Alp+Mnd+Mno+Mnp+Sio+NH4+NO3+Pt+Pi,
               data=env1,distance = "bray")
summary(mod_e)
RsquareAdj(mod_e)

vif.cca(mod_e)
# depth           pH Conductivity     moisture            N          N15            C           CN 
# 2.292997     4.015903     3.198264     3.911302    18.045197     2.551121    16.518521     2.295329 
# Feo          Fed       FeoFed          Fep          Alo          Alp          Mnd          Mno 
# 12.644294     4.234230     7.220586     9.675572    33.649155    29.812333     2.023044    28.065596 
# Mnp          Sio          NH4          NO3           Pt           Pi 
# 23.247352     3.310152     1.732560     1.698673     1.552453     3.792609 
#result
depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Mnd+Sio+NH4+NO3+Pt+Pi

#test result(remove coliner)
mod_e <- dbrda(asv3~depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Mnd+Sio+NH4+NO3+Pt+Pi,
               data=env1,distance = "bray")
summary(mod_e)
RsquareAdj(mod_e)
vif.cca(mod_e)

summary
all:depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi+dna_typen
i  :depth+pH+Conductivity+moisture+N15+   Feo+FeoFed+Fep+Sio+Mnd+        Pt+Pi
e  :depth+pH+Conductivity+moisture+N15+CN+Fed+FeoFed+Fep+Sio+Mnd+NH4+NO3+Pt+Pi