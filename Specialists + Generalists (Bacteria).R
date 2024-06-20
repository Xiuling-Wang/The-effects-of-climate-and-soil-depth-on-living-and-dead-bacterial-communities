library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(rdacca.hp)
library(RColorBrewer)
library(labdsv)#Indicator Value
# load data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1,
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_2024.csv", 
                row.names = 1,
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  row.names = 1,
                  header = TRUE)
# 01 iDNA specialist-------------------------------------------------------------
#clean env
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv),
         depth < 40,
         dna_type == "iDNA") %>% 
  na.omit() %>% 
  data.frame()

#clean asv
filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]
#clean tax
filt3 <- rownames(tax) %in% rownames(asv2)
tax1 <-tax[filt3,]
# new asv (grouped)
group_info <- data.frame(col_name = rownames(env1),
                         group = env1$site_dna_depth)

asv2.5 <- asv2 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  gather(key = "col_name", value = "value", -id) %>%  # Convert to long format
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value)) %>%  # Compute the mean for each group
  spread(key = "group", value = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env1$site_dna_depth)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]
# trans and relative abundance
# total_reads_i <- colSums(asv2.8)
# relative_abundance_i<- 100*asv2.8 / total_reads_i
# relative_abundance_i$mean_ra <- rowMeans(relative_abundance_i) 
#以上3行为原本运行的错误代码（已在revised version 更正图片）
asv2.8_t <- t(asv2.8)
total_reads <- rowSums(asv2.8_t)
relative_abundance <- 100*asv2.8_t / total_reads
relative_abundance_i <- t(relative_abundance) %>% 
  as.data.frame()
relative_abundance_i$mean_ra <- rowMeans(relative_abundance_i) 
#set threshold
filt_asv_i <- relative_abundance_i %>%
  filter(mean_ra > 0.1) %>% #keep over 0.001 asv
  select(-"mean_ra")
filt_asv_i_t <- t(filt_asv_i)
#create group for analyse
group_info2 <- env1 %>%
  select(site_dna_depth, site) %>%
  distinct(site_dna_depth, .keep_all = TRUE) %>% 
  as.data.frame()
rownames(group_info2) <- group_info2[,1]
group_info2 <- group_info2[,-1,drop = F]

group_info2$site <- as.factor(group_info2$site)
levels_info <- levels(group_info2$site)
#indval analyse
iva <- indval(filt_asv_i_t,group_info2$site)
# relfrq = relative frequency of species in classes
# relabu = relative abundance of species in classes
# indval = the indicator value for each species
# maxcls = the class each species has maximum indicator value for
# indcls = the indicator value for each species to its maximum class
# pval😀 = the probability of obtaining as high an indicator values as observed over the
# specified iterations
gr <- levels_info[iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]]
#gr <- iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]
iv <- iva$indcls[iva$pval <= 0.05 & iva$indcls > 0.8]
pv <- iva$pval[iva$pval <= 0.05 & iva$indcls > 0.8]
fr <- apply(filt_asv_i_t > 0, 2, sum)[iva$pval <= 0.05 & iva$indcls > 0.8]

indvalsummary <- data.frame(group = gr, 
                            indval = iv, 
                            pvalue = pv, 
                            freq = fr)
indvalsummary <- indvalsummary[order(indvalsummary$fr, -indvalsummary$indval, decreasing = TRUE),]
indvalsummary
tab_s1 <- indvalsummary[order(indvalsummary$group), ]

tab_s1

tax_tmp <- tax %>% 
  filter(rownames(.) %in% rownames(tab_s1)) %>% 
  as.data.frame()

TableS1 <- merge(tab_s1,
                 tax_tmp,
                 by = "row.names", 
                 all.x = TRUE) %>% 
  arrange(Phylum)

head(TableS1)
write.table(TableS1, 
            file = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/Figure for papers😊/table/Table_S1_indvalsummary.txt", 
            sep = "\t", 
            row.names = T)

filt4 <- rownames(filt_asv_i) %in% rownames(indvalsummary)
specialist_i <- filt_asv_i[filt4,]
filt5 <- rownames(tax1) %in% rownames(specialist_i)
tax2 <- tax1[filt5,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- specialist_i %>%
  rownames_to_column(var = "asv_id")
tax_table <- tax2 %>%
  rownames_to_column(var = "asv_id")
# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
)
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)
# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)
unique(sorted_data$sample_id)
tmp123 <- c("AZiDNA0_5cm","AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm",
            "SGiDNA0_5cm","SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm",
            "LCiDNA0_5cm","LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm",
            "NBiDNA0_5cm","NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZ 0-5cm","AZ 5-10cm","AZ 10-20cm","AZ 20-40cm",
                                       "SG 0-5cm","SG 5-10cm","SG 10-20cm","SG 20-40cm",
                                       "LC 0-5cm","LC 5-10cm","LC 10-20cm","LC 20-40cm",
                                       "NB 0-5cm","NB 5-10cm","NB 10-20cm","NB 20-40cm"),
                            ordered = TRUE)
         )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
summary(filtered_data)
unique(long_data$Phylum)
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
    "p__Chloroflexi"="#27cdb2",
    "p__Proteobacteria"="#afa523",
    "p__Actinobacteriota"="#b2df8a",
    "p__Verrucomicrobiota"="#238f1c",
    "p__GAL15"="#1f78b4",
    "p__Acidobacteriota"="#f1393b",
    "p__Nitrospirota"="#fb9a99",
    "p__Methylomirabilota"="#ffcd8d",
    "p__Myxococcota"="#ff7f00",
    "p__Gemmatimonadota"="#e4c4f2",
    "Environment"="#2c59ed",
    "p__NA"="#dfd3c7",
    "p__RCP2-54"="#e37bf2",
    "p__Planctomycetota"="#b15928")

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/08 specialist and generalist/Specialist(iDNA pool).pdf",
    width = 13,
    height = 8)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(0.1,1,2,4,8),
                        limits = c(0,10)) +
  labs(x = "Sample", y = "", title = "Specialist (iDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 02 iDNA generalist --------------------------------------------------------------
asv_table <- asv2.8
# 计算 ASV 在样本中出现的频率
asv_freq <- asv_table %>%
  apply(1, function(x) sum(x > 0)) %>%
  as.data.frame() %>%
  rename(frequency = 1) %>%
  rownames_to_column(var = "ASV")
# 计算 ASV 在所有站点中的相对丰度
asv_relab <- asv_table %>%
  apply(1, function(x) sum(x)) %>%
  as.data.frame() %>%
  rename(relative_abundance = 1) %>%
  rownames_to_column(var = "ASV") %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance) * 100)
# 合并两个数据框
asv_freq_relab <- inner_join(asv_freq, asv_relab, by = "ASV")
# 筛选出高度广泛分布的 ASV
generalist_asv <- asv_freq_relab %>%
  filter(frequency / ncol(asv_table) > 0.75 & relative_abundance > 0.1)

filt6 <- rownames(filt_asv_i) %in% generalist_asv$ASV
generalist_tab <- filt_asv_i[filt6,]

filt7 <- rownames(tax1) %in% rownames(generalist_tab)
tax3 <- tax1[filt7,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- generalist_tab %>%
  rownames_to_column(var = "asv_id")
tax_table <- tax3 %>%
  rownames_to_column(var = "asv_id")
# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)

# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)
unique(sorted_data$sample_id)
tmp123 <- c("AZiDNA0_5cm","AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm",
            "SGiDNA0_5cm","SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm",
            "LCiDNA0_5cm","LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm",
            "NBiDNA0_5cm","NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZ 0-5cm","AZ 5-10cm","AZ 10-20cm","AZ 20-40cm",
                                       "SG 0-5cm","SG 5-10cm","SG 10-20cm","SG 20-40cm",
                                       "LC 0-5cm","LC 5-10cm","LC 10-20cm","LC 20-40cm",
                                       "NB 0-5cm","NB 5-10cm","NB 10-20cm","NB 20-40cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
summary(filtered_data)
unique(filtered_data$Phylum)
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")
# 使用自定义颜色调色板
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/08 specialist and generalist/Generalist (iDNA pool).pdf",
    width = 11,
    height = 6)

ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(0.1,1,2,4,8),
                        limits = c(0,8)) +
  labs(x = "Sample", y = "", title = "Generalist (iDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()

rm(list = ls())

# load data ---------------------------------------------------------------
asv <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/processed/bac_200_9435rarefild/asv_rare_200.txt",
                  row.names = 1,
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/env_2024.csv", 
                row.names = 1,
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/01 🌸Chile Bacteria/raw/tax.txt",
                  row.names = 1,
                  header = TRUE)
# 03 eDNA specialist-------------------------------------------------------------
# clean env
env1 <- env %>%
  filter(rownames(env) %in% colnames(asv),
         depth < 40,
         dna_type == "eDNA") %>% 
  na.omit() %>% 
  data.frame()
#clean asv
filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]
filt2 <- rowSums(asv1) > 0
asv2 <- asv1[filt2,]
#clean tax
filt3 <- rownames(tax) %in% rownames(asv2)
tax1 <-tax[filt3,]
# new asv (grouped)
group_info <- data.frame(col_name = rownames(env1),group = env1$site_dna_depth)

asv2.5 <- asv2 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  gather(key = "col_name", value = "value", -id) %>%  # Convert to long format
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value)) %>%  # Compute the mean for each group
  spread(key = "group", value = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env1$site_dna_depth)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]
# trans and relative abundance
# total_reads_i <- colSums(asv2.8)
# relative_abundance_i<- 100*asv2.8 / total_reads_i
# relative_abundance_i$mean_ra <- rowMeans(relative_abundance_i) 
#
asv2.8_t <- t(asv2.8)
total_reads <- rowSums(asv2.8_t)
relative_abundance <- 100*asv2.8_t / total_reads
relative_abundance_i <- t(relative_abundance) %>% 
  as.data.frame()
relative_abundance_i$mean_ra <- rowMeans(relative_abundance_i) 
#set threshold
filt_asv_i <- relative_abundance_i %>%
  filter(mean_ra > 0.1) %>% #keep over 0.001 asv
  select(-"mean_ra")

filt_asv_i_t <- t(filt_asv_i)
#create group for analyse
group_info2 <- env1 %>%
  select(site_dna_depth, site) %>%
  distinct(site_dna_depth, .keep_all = TRUE) %>% 
  as.data.frame()
rownames(group_info2) <- group_info2[,1]
group_info2 <- group_info2[,-1,drop = F]
#indval analyse
iva <- indval(filt_asv_i_t,group_info2$site)
# specified iterations
gr <- iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]
iv <- iva$indcls[iva$pval <= 0.05 & iva$indcls > 0.8]
pv <- iva$pval[iva$pval <= 0.05 & iva$indcls > 0.8]
fr <- apply(filt_asv_i_t > 0, 2, sum)[iva$pval <= 0.05 & iva$indcls > 0.8]

indvalsummary <- data.frame(group = gr, indval = iv, pvalue = pv, freq = fr)
indvalsummary <- indvalsummary[order(indvalsummary$fr, -indvalsummary$indval, decreasing = TRUE),]
indvalsummary

filt4 <- rownames(filt_asv_i) %in% rownames(indvalsummary)
specialist_i <- filt_asv_i[filt4,]

filt5 <- rownames(tax1) %in% rownames(specialist_i)
tax2 <- tax1[filt5,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- specialist_i %>%
  rownames_to_column(var = "asv_id")

tax_table <- tax2 %>%
  rownames_to_column(var = "asv_id")

# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum, Class,Order,Family,Genus,asv_id)

# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)

unique(sorted_data$sample_id)
tmp123 <- c("AZeDNA0_5cm","AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm",
            "SGeDNA0_5cm","SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm",
            "LCeDNA0_5cm","LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm",
            "NBeDNA0_5cm","NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZ 0-5cm","AZ 5-10cm","AZ 10-20cm","AZ 20-40cm",
                                       "SG 0-5cm","SG 5-10cm","SG 10-20cm","SG 20-40cm",
                                       "LC 0-5cm","LC 5-10cm","LC 10-20cm","LC 20-40cm",
                                       "NB 0-5cm","NB 5-10cm","NB 10-20cm","NB 20-40cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
summary(filtered_data)
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")
# 使用自定义颜色调色板

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/08 specialist and generalist/Specialist (eDNA pool).pdf",
    width = 12,
    height = 8)

ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(0.1,1,2.5,5),
                        limits = c(0,5)) +
  labs(x = "Sample", y = "", title = "Specialist (eDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 04 eDNA generalist --------------------------------------------------------------
asv_table <- asv2.8
asv_freq <- asv_table %>%
  apply(1, function(x) sum(x > 0)) %>%
  as.data.frame() %>%
  rename(frequency = 1) %>%
  rownames_to_column(var = "ASV")
# 计算 ASV 在所有站点中的相对丰度
asv_relab <- asv_table %>%
  apply(1, function(x) sum(x)) %>%
  as.data.frame() %>%
  rename(relative_abundance = 1) %>%
  rownames_to_column(var = "ASV") %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance) * 100)
# 合并两个数据框
asv_freq_relab <- inner_join(asv_freq, asv_relab, by = "ASV")
# 筛选出高度广泛分布的 ASV
generalist_asv <- asv_freq_relab %>%
  filter(frequency / ncol(asv_table) > 0.75 & relative_abundance > 0.1)

filt6 <- rownames(filt_asv_i) %in% generalist_asv$ASV
generalist_tab <- filt_asv_i[filt6,]

filt7 <- rownames(tax1) %in% rownames(generalist_tab)
tax3 <- tax1[filt7,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- generalist_tab %>%
  rownames_to_column(var = "asv_id")

tax_table <- tax3 %>%
  rownames_to_column(var = "asv_id")
# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  gather(key = "sample_id", 
         value = "relative_abundance", -asv_id,-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum, Class,Order,Family,Genus,asv_id)

# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)

unique(sorted_data$sample_id)
tmp123 <- c("AZeDNA0_5cm","AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm",
            "SGeDNA0_5cm","SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm",
            "LCeDNA0_5cm","LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm",
            "NBeDNA0_5cm","NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZ 0-5cm","AZ 5-10cm","AZ 10-20cm","AZ 20-40cm",
                                       "SG 0-5cm","SG 5-10cm","SG 10-20cm","SG 20-40cm",
                                       "LC 0-5cm","LC 5-10cm","LC 10-20cm","LC 20-40cm",
                                       "NB 0-5cm","NB 5-10cm","NB 10-20cm","NB 20-40cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
summary(filtered_data)
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
custom_color_palette <- c("p__Firmicutes"="#a6cee3",
                          "p__Chloroflexi"="#27cdb2",
                          "p__Proteobacteria"="#afa523",
                          "p__Actinobacteriota"="#b2df8a",
                          "p__Verrucomicrobiota"="#238f1c",
                          "p__GAL15"="#1f78b4",
                          "p__Acidobacteriota"="#f1393b",
                          "p__Nitrospirota"="#fb9a99",
                          "p__Methylomirabilota"="#ffcd8d",
                          "p__Myxococcota"="#ff7f00",
                          "p__Gemmatimonadota"="#e4c4f2",
                          "Environment"="#2c59ed",
                          "p__NA"="#dfd3c7",
                          "p__RCP2-54"="#e37bf2",
                          "p__Planctomycetota"="#b15928")

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/01 🌸Chile Bacteria/08 specialist and generalist/Generalist (eDNA pool).pdf",
    width = 13,
    height = 5)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(0.1,0.5,1,2.5),
                        limits = c(0,2.5)) +
  labs(x = "Sample", y = "", title = "Generalist (eDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
rm(list = ls())