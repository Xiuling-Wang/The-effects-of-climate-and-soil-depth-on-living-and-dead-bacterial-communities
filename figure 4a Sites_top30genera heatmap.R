library(pheatmap)
library(tidyverse)
library(magrittr)
library(forcats)
library(viridis)
library(vegan)
# Load data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/tax.txt",
                  header = TRUE)
# clean raw data ----------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id%in%colnames(asv))

tax1 <- tax %>% 
  filter(tax$asv_id%in%asv$asv_id)

# get top 30 Genus list first ---------------------------------------------
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

top30genus <- top_genus$Genus[2:31]
top30genus

plot_data1 <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun')%>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(site,dna_type,Class,Genus) %>% 
  summarise(GenusAbun = sum(Abun)) %>%
  ungroup %>%
  group_by(site,dna_type)%>%
  mutate(RelAbun = 100*GenusAbun/sum(GenusAbun)) %>%
  ungroup()
# merge site&dna_type and transfer wider data
plot_data2 <- plot_data1 %>%
  filter(Genus %in% top30genus) %>%
  mutate(tmp=paste0(.$site, .$dna_type)) %>%
  pivot_wider(id_cols = Genus, names_from = tmp, values_from = RelAbun)

# Method1 
tmp1 <- plot_data2$Genus
plot_data3 <- as.matrix(plot_data2[2:9])
rownames(plot_data3) <- tmp1

# Method2
plot_data3 <- as.matrix(plot_data2[2:9]) %>% `rownames<-`(plot_data2$Genus)

# plot 
row_anno <- tax %>% 
  select(Class, Genus) %>%
  distinct() %>%
  filter(Genus %in% top30genus)

#row_anno[row_anno$Class=='c__NA','Phylum'] <- 'Unassigned'
row_anno %<>% distinct()

tmp2 <- data.frame(row_anno$Class) 
rownames(tmp2) <- row_anno$Genus
row_anno <- tmp2
names(row_anno) <- 'Class'

my_colour = list(
  Class = RColorBrewer::brewer.pal(8, 'Paired')%>% `names<-`(unique(row_anno$Class))
)

drows = vegdist(plot_data3, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
dcols =vegdist(t(plot_data3), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

pdf(file = '/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/04_heatmap/site_top30genus.pdf',
    height = 8,
    width = 10)
pheatmap(plot_data3, 
                  scale='row', 
                  cutree_cols = 4,
                  cluster_rows = F,
                  cluster_cols = TRUE, 
                  #clustering_distance_rows = drows,
                  clustering_distance_cols = dcols,
                  clustering_method = "ward.D2",
                  annotation_row = row_anno, 
                  annotation_colors  = my_colour
                  )
dev.off()

# HZY script --------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(pheatmap)
library(forcats)

setwd('G:\\temp\\demo\\scripts')
## Load data
asv <- read.table("../data/raw/subsample_asv.txt", header = TRUE, sep='\t',check.names = FALSE)
env <- read.table("../data/raw/env_03_08.txt", header = TRUE, sep='\t',check.names = FALSE)
tax <- read.table("../data/raw/asv_taxnomy.txt", header=TRUE, sep='\t')

plot_data1 <- asv %>% 
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'Abun')%>% 
  filter(Abun>0) %>%
  left_join(tax) %>%
  left_join(env) %>%
  group_by(site, dna_type, Class) %>% # !! remove phylum
  summarise(ClassAbun = sum(Abun)) %>%
  ungroup() %>%
  group_by(site, dna_type)%>%
  mutate(RelAbun = ClassAbun/sum(ClassAbun)) %>%
  ungroup()

top20class <- c('Actinobacteria','Thermoleophilia','Gammaproteobacteria','Alphaproteobacteria',
                'AD3','Verrucomicrobiae','Bacilli','Planctomycetes','Unassigned','Acidobacteriae',
                'Vicinamibacteria','Acidimicrobiia','MB-A2-108','Blastocatellia','Gemmatimonadetes',
                'Rubrobacteria','Bacteroidia','Chloroflexia','Methylomirabilia','KD4-96')

plot_data2 <- plot_data1 %>%
  filter(Class %in% top20class) %>%
  mutate(tmp=paste0(.$site, .$dna_type)) %>%
  pivot_wider(id_cols = Class, names_from = tmp, values_from = RelAbun)

# Method1
tmp1 <- plot_data2$Class
plot_data3 <- as.matrix(plot_data2[2:9])
rownames(plot_data3) <- tmp1

# Method2
plot_data3 <- as.matrix(plot_data2[2:9]) %>% `rownames<-`(plot_data2$Class)

# plot 
row_anno <- tax %>% 
  select(Phylum, Class) %>%
  distinct() %>%
  filter(Class %in% top20class)
row_anno[row_anno$Class=='Unassigned','Phylum'] <- 'Other'
row_anno %<>% distinct()

tmp2 <- data.frame(row_anno$Phylum)
rownames(tmp2) <- row_anno$Class
row_anno <- tmp2
names(row_anno) <- 'Phylum'

my_colour = list(
  Phylum = RColorBrewer::brewer.pal(12, 'Paired')%>% `names<-`(unique(row_anno$Phylum))
)

p_tmp <- pheatmap(plot_data3, scale='row', cutree_cols = 4, annotation_row = row_anno, annotation_colors  = my_colour)
