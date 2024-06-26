library(vegan)   
library(ggplot2)
library(tidyverse)

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                row.names = 1, 
                header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_manage_220812.csv",
                header = TRUE)
# clean -------------------------------------------------------------------
env1 <- env %>% 
  filter(sample_id %in% colnames(asv))
asv1 <- data.frame(t(asv))

bray_dis <- vegdist(asv1, method = 'bray')
nmds_dis <- metaMDS(bray_dis, k = 2)# Nonmetric Multidimensional Scaling (NMDS)

nmds_dis_site <- data.frame(nmds_dis$points)
sample_id <- row.names(nmds_dis_site)

nmds_dis_site1 <- cbind(sample_id,nmds_dis_site)
nmds_dis_site2 <- left_join(nmds_dis_site1,env,by=c("sample_id"))

nmds_dis_site2$site <- factor(nmds_dis_site2$site,
                              levels = c('AZ','SG','LC','NB'))

nmds_dis_site2$dna_type <- factor(nmds_dis_site2$dna_type,
                              levels = c('iDNA','eDNA'))

nmds_dis_site2$depth_order <- factor(nmds_dis_site2$depth_order,
                                  levels = c("1","2","3","4","5","6","7","8","9","10","11","12"),
                                  labels = c("0-5 cm","5-10 cm","10-20 cm","20-40 cm","40-60 cm","60-80 cm","80-100 cm",
                                             "100-120 cm","120-140 cm","140-160 cm","160-180 cm","180-200 cm"))

print(nmds_dis$stress)#fill in the plot

# nmds_all ----------------------------------------------------------------

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/05_NMDS/i&e_nmds.pdf'),
    width=8,height=6)

ggplot(data = nmds_dis_site2, 
       aes(MDS1, MDS2)) +
  geom_point(aes(color = site,
                          shape = dna_type,
                          alpha = depth_order),
             size = 3) +
  stat_ellipse(data = nmds_dis_site2 %>% 
                 filter(dna_type=="iDNA"),
               aes(col = site),
               linetype = 1,
               level = 0.95, 
               show.legend = FALSE) +
  stat_ellipse(data = nmds_dis_site2 %>% 
                 filter(dna_type=="eDNA"),
               aes(col = site),
               linetype = 2,
               level = 0.95, 
               show.legend = FALSE) +
  scale_shape_manual(values=c(16, 17)) +
  scale_color_manual(values =RColorBrewer::brewer.pal(4, "Set1")) +
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),
        legend.key = element_rect(fill = "transparent")
        ) +
  #legend.position = 'none'+
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  theme(axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12))+
  guides(fill=guide_legend(title="Sites"))+
  ggtitle("i&e DNA community Stress = 0.1847374")

  dev.off()
  
  

# iDNA --------------------------------------------------------------------
  env1 <- env %>% 
    filter(sample_id %in% colnames(asv),dna_typen==0)
filter1=colnames(asv)%in%env1$sample_id

asv0.5=asv[,filter1]
filt2=rowSums(asv0.5)>0

asv0.8=asv0.5[filt2,]

asv1 <- data.frame(t(asv0.8))

bray_dis <- vegdist(asv1, method = 'bray')#computes dissimilarity indices
nmds_dis <- metaMDS(bray_dis, k = 2)#performs Nonmetric Multidimensional Scaling (NMDS)
nmds_dis_site <- data.frame(nmds_dis$points)

sample_id <- row.names(nmds_dis_site)
nmds_dis_site1 <- cbind(sample_id,nmds_dis_site)
nmds_dis_site2 <- left_join(nmds_dis_site1,env,by=c("sample_id"))

nmds_dis_site2$site <- factor(nmds_dis_site2$site,
                              levels = c('PandeAzucar',
                                         'SantaGracia',
                                         'LaCampana',
                                         'Nahuelbuta'))

print(nmds_dis$stress)#fill in plot with the stress

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/05_NMDS/iDNA_nmds.pdf'),
    width=8,height=6)

ggplot(data = nmds_dis_site2, aes(MDS1, MDS2)) +
  geom_point(size=1.5,aes(color = site)) +
  stat_ellipse(aes(fill = interaction(site)),#添加置信椭圆 
               geom = 'polygon', level = 0.95, 
               alpha = 0.1, show.legend = FALSE) +
  scale_shape_manual(values=c(16, 1)) +
  scale_color_manual(values =RColorBrewer::brewer.pal(4, "Set1")) +
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  #legend.position = 'none'+
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  theme(axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12))+
  ggtitle("iDNA community Stress = 0.1836138")
dev.off()

# eDNA --------------------------------------------------------------------

env1 <- env %>% 
  filter(sample_id %in% colnames(asv),dna_typen==1)

filter1=colnames(asv)%in%env1$sample_id

asv0.5=asv[,filter1]
filt2=rowSums(asv0.5)>0

asv0.8=asv0.5[filt2,]

asv1 <- data.frame(t(asv0.8))

bray_dis <- vegdist(asv1, method = 'bray')#computes dissimilarity indices
nmds_dis <- metaMDS(bray_dis, k = 2)#performs Nonmetric Multidimensional Scaling (NMDS)
nmds_dis_site <- data.frame(nmds_dis$points)

sample_id <- row.names(nmds_dis_site)
nmds_dis_site1 <- cbind(sample_id,nmds_dis_site)
nmds_dis_site2 <- left_join(nmds_dis_site1,env,by=c("sample_id"))

nmds_dis_site2$site <- factor(nmds_dis_site2$site,
                              levels = c('PandeAzucar','SantaGracia','LaCampana','Nahuelbuta'))
print(nmds_dis$stress)#fill in the data in the plot

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/05_NMDS/eDNA_nmds.pdf'),
    width=8,height=6)

ggplot(data = nmds_dis_site2, aes(MDS1, MDS2)) +
  geom_point(size=1.5,aes(color = site)) +
  stat_ellipse(aes(fill = interaction(site)),#添加置信椭圆 
               geom = 'polygon', level = 0.95, 
               alpha = 0.1, show.legend = FALSE) +
  scale_shape_manual(values=c(16, 1)) +
  scale_color_manual(values =RColorBrewer::brewer.pal(4, "Set1")) +
  theme(panel.grid.major = element_line(color = 'gray', 
                                        size = 0.2), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  #legend.position = 'none'+
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  theme(axis.text.x= element_text(size=12),
        axis.text.y= element_text(size=12))+
  ggtitle("eDNA community Stress = 0.1561664")

dev.off()
