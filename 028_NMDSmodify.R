library(vegan)   
library(ggplot2)
library(tidyverse)

asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/processed/bac_200_9435rarefild/asv_rare_200.txt",
                row.names = 1, header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/data/raw/env_manage_220812.csv",
                header = TRUE)
# clean -------------------------------------------------------------------
env1 <- env %>% 
  filter(sample_id %in% colnames(asv))
asv1 <- data.frame(t(asv))

bray_dis <- vegdist(asv1, method = 'bray')
nmds_dis <- metaMDS(bray_dis, k = 2)#performs Nonmetric Multidimensional Scaling (NMDS)

nmds_dis_site <- data.frame(nmds_dis$points)
sample_id <- row.names(nmds_dis_site)

nmds_dis_site1 <- cbind(sample_id,nmds_dis_site)
nmds_dis_site2 <- left_join(nmds_dis_site1,env,by=c("sample_id"))

nmds_dis_site2$site <- factor(nmds_dis_site2$site,
                              levels = c('AZ','SG','LC','NB'))

nmds_dis_site2$dna_type <- factor(nmds_dis_site2$dna_type,
                              levels = c('iDNA','eDNA'))
print(nmds_dis$stress)#fill in the plot

# nmds_all ----------------------------------------------------------------

pdf(file.path('/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/👨‍🎓phd_project/phd_es16s/reports/05_NMDS/i&e_nmds.pdf'),
    width=8,height=6)

ggplot(data = nmds_dis_site2, 
       aes(MDS1, MDS2)) +
  geom_point(size=1.5,aes(color = site,
                          shape = dna_type)) +
  stat_ellipse(aes(fill = interaction(dna_type, site)),#添加置信椭圆 
               geom = 'polygon', level = 0.95, 
               alpha = 0.1, show.legend = FALSE) +
  scale_shape_manual(values=c(16, 1)) +
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
############分别为IDNA和eDNA做DBRDA
###### i DNA
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
                              levels = c('PandeAzucar','SantaGracia','LaCampana','Nahuelbuta'))

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


################ eDNA
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
